"""Optional carried-source audit of exact same-cell AG row grouping.

Reuses the existing read-only observer and checks untapped control parity.
Exported float operands are exact rational inputs to the row factorization;
80-digit singular values are non-promoting feasibility diagnostics only.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from .ag_readout_source_diagnostic import HEADER, REPO, instrument
from .lin_path_certificate import inverse
from .matrix_certificates import matmul
from .moving_pivots import same_prediction_cell_groups


def verify_diagnostic(summary):
    """Verify provenance and non-promoting scope, not the unretained trace."""
    source = (REPO/HEADER).read_text()
    expected = {'shipping_header_sha256': hashlib.sha256(source.encode()).hexdigest(),
                'observer_header_sha256': hashlib.sha256(instrument(source).encode()).hexdigest(),
                'driver_sha256': hashlib.sha256((REPO/'tools/stability/ag_readout_source.cpp').read_bytes()).hexdigest()}
    if any(summary.get(k) != v for k, v in expected.items()):
        raise ValueError('same-cell source/observer/driver provenance changed')
    if summary.get('qualification') != 'OU3_CARRIED_SAME_CELL_GROUP_DIAGNOSTIC_V1':
        raise ValueError('same-cell diagnostic qualification changed')
    for flag in ('source_uniform_verified', 'all_time_magnetic_service_certified', 'theorem_closed'):
        if summary.get(flag) is not False:
            raise ValueError('finite same-cell diagnostic cannot promote '+flag)
    cases = summary.get('cases', [])
    if [c.get('input_profile') for c in cases] != ['0', 'wave']:
        raise ValueError('both carried source profiles required')
    for case in cases:
        if (case.get('literal_terminal_parity') is not True or
                case.get('row_defect_exactly_zero') is not True or
                case.get('groups', 0) < 2 or
                case.get('source_uniform_verified') is not False or
                case.get('real_trajectory_enclosed') is not False):
            raise ValueError('same-cell source audit scope/parity failed')
    return True


def analyze(trace):
    import mpmath as mp
    groups = same_prediction_cell_groups(trace['events'])
    if len(groups) < 2:
        raise ValueError('two actual complete row groups required')
    first, last = groups[0], groups[-1]
    between = matmul(last['anchor_transport'], inverse(first['anchor_transport']))
    # Exact E=0 is checked in the ORIGINAL historical-root coordinates.
    rows = first['raw_rows'] + last['raw_rows']
    local_rows = matmul(rows, inverse(first['anchor_transport']))
    expected = [r + [F(0)]*3 for r in first['C']] + matmul(last['C'], between[:3])
    if local_rows != expected:
        raise ArithmeticError('two-group chronological identity failed')
    with mp.workdps(80):
        scalar = lambda x: mp.mpf(x.numerator)/x.denominator
        matrix = lambda x: mp.matrix([[scalar(v) for v in row] for row in x])
        c = min(min(mp.svd(matrix(g['C']), compute_uv=False)) for g in (first, last))
        a = max(mp.svd(matrix(between)[:3, :3], compute_uv=False))
        b = min(mp.svd(matrix(between)[:3, 3:], compute_uv=False))
        s = c/(1+(a+1)/b) if b > 0 else mp.mpf(0)
        actual = min(mp.svd(matrix(local_rows), compute_uv=False))
        if not 0 < s <= actual:
            raise ArithmeticError('supplied source two-group budget infeasible')
        values = {'minimum_group_attitude_singular_value': c,
                  'interanchor_attitude_norm': a, 'interanchor_gyro_singular_value': b,
                  'two_group_six_column_budget': s, 'two_group_actual_singular_value': actual,
                  'all_six_pivot_budget': s/4}
        return {'groups': len(groups), 'selected_acc_events': [first['acc_event'], last['acc_event']],
                'selected_mag_events': [first['mag_event'], last['mag_event']],
                'row_defect_exactly_zero': True,
                **{k: mp.nstr(v, 40) for k, v in values.items()},
                'exported_trace_sha256': hashlib.sha256(json.dumps(trace, sort_keys=True).encode()).hexdigest(),
                'source_uniform_verified': False, 'real_trajectory_enclosed': False}


def run(eigen):
    source = (REPO/HEADER).read_text()
    driver = REPO/'tools/stability/ag_readout_source.cpp'
    cases = []
    with tempfile.TemporaryDirectory(prefix='ou3-group-source-') as directory:
        tmp = Path(directory)
        include = tmp/'kalman_ou_iii'
        include.mkdir()
        observed_source = instrument(source)
        (include/HEADER.name).write_text(observed_source)
        binaries = []
        for name, inc in [('observed', ['-I'+str(tmp)]), ('control', [])]:
            binary = tmp/name
            subprocess.run(['g++', '-O2', '-std=c++20', *inc, '-I'+str(REPO/'src'),
                            '-isystem', str(eigen), str(driver), '-o', str(binary)], check=True)
            binaries.append(binary)
        for profile in ('0', 'wave'):
            observed, control = [json.loads(subprocess.check_output([str(p), profile], text=True))
                                 for p in binaries]
            if any(control[k] != observed[k] for k in control if k != 'events'):
                raise ArithmeticError('observer changed literal source output')
            cases.append({'input_profile': profile, 'literal_terminal_parity': True,
                          'live_step': observed['live_step'], 'refined_step': observed['refined_step'],
                          'active_step': observed['active_step'], **analyze(observed)})
    return {'qualification': 'OU3_CARRIED_SAME_CELL_GROUP_DIAGNOSTIC_V1',
            'shipping_header_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'observer_header_sha256': hashlib.sha256(observed_source.encode()).hexdigest(),
            'driver_sha256': hashlib.sha256(driver.read_bytes()).hexdigest(),
            'decimal_digits': 80, 'cases': cases, 'source_uniform_verified': False,
            'all_time_magnetic_service_certified': False, 'theorem_closed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.eigen), indent=2, sort_keys=True)+'\n')
