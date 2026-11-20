"""Compare the original and transported TFG shared-mode response on matched records."""
import argparse
import random
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from marine_mag_filter_replays import REFERENCE_CSV_SHA256_BY_CASE, ROOT, digest, replay, write

VALIDATION_SEEDS = (180161, 211313, 233021, 257371, 280001, 307009, 330017, 362161)


def replay_score(task):
    case, record, mode, seed = task
    result = replay(("tfg_shared", case, record, mode), seed=seed)
    return dict(seed=seed if seed is not None else "default", **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", choices=("released", "validation", "replication"))
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    # The replication panel is fixed independently of the eight validation draws.
    seeds = {"released": (None,), "validation": VALIDATION_SEEDS,
             "replication": random.Random(7716020).sample(range(400000, 4000000), 32)}[args.panel]
    tasks = []
    for case, entries in REFERENCE_CSV_SHA256_BY_CASE.items():
        for record, checksum in sorted(entries.items()):
            if args.panel != "released" and ("pmstokes" not in record or "H4.000" in record):
                continue
            path = ROOT / ".sim-data/cases" / case / record
            if digest(path) != checksum:
                raise ValueError(f"Unexpected wave input: {path}")
            tasks.extend((case, path, mode, seed) for seed in seeds for mode in ("native", "production"))
    with ThreadPoolExecutor(args.jobs) as pool:
        scores = list(pool.map(replay_score, tasks))
    write(scores, args.out)


if __name__ == "__main__":
    main()
