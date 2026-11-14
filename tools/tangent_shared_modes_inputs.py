"""Verify historical raw inputs against the published score-cohort receipts."""
import csv
import hashlib
from pathlib import Path


def verify_raw_inputs(raw, manifest, score_path):
    """Check all raw bytes and their aggregated scores before report generation.

    Locate relocated raw files by their original basename; the manifest retains
    the producer's absolute command paths. Never replace the producer receipts
    with newly computed hashes from an unrelated raw directory.
    """
    rows, seen = [], set()
    for receipt in manifest['commands']:
        command = receipt['command']
        case = command[command.index('--scenario')+1]
        seed = int(command[command.index('--seed')+1])
        model = ('population' if '--population-prior' in command else
                 'noise' if '--matched-noise' in command else 'native')
        name = f'{case}_{seed}_{model}'
        if Path(command[command.index('--output')+1]).name != name or name in seen:
            raise ValueError(f'invalid or duplicate experiment receipt: {name}')
        seen.add(name)
        for suffix, key in (('_scores.csv', 'scores_sha256'),
                            ('_trace.csv', 'trace_sha256'),
                            ('_model.json', 'model_sha256')):
            path = Path(raw)/(name+suffix)
            expected = receipt.get(key)
            if not expected or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f'raw input differs from scored receipt: {path}')
        with (Path(raw)/(name+'_scores.csv')).open(newline='') as stream:
            experiment_rows = list(csv.DictReader(stream))
        for row in experiment_rows:
            if row['scenario'] != case or int(row['seed']) != seed:
                raise ValueError(f'score identity differs from receipt: {name}')
            rows.append(dict(row, noise_profile=model))
    with Path(score_path).open(newline='') as stream:
        published = list(csv.DictReader(stream))
    # Compare exact values including failure markers, independently of row order
    # and CSV line endings. This binds the aggregate to the receipt score hashes.
    def canonical(values):
        return sorted(tuple(sorted(row.items())) for row in values)
    if canonical(rows) != canonical(published):
        raise ValueError('published scores differ from raw experiment receipts')
    return len(seen)
