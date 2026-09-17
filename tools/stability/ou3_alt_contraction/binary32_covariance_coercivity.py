"""Exact coercivity of finite SPD binary32 covariance, conditional on preservation.

This proves a format-wide matrix implication, not shipping SPD preservation.
No sampled eigenvalues, diagonal loading, covariance floors or fitted metric.
"""
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

QUANTUM = F(1, 2**149)
MAX_FINITE = F((2**24-1)*2**104)
QUALIFICATION = 'OU3_ALT_FINITE_SPD_BINARY32_COERCIVITY_V1'


def bounds(dimension=21):
    if type(dimension) is not int or dimension < 1:
        raise ValueError('positive integer dimension required')
    upper = dimension*MAX_FINITE
    lower = QUANTUM**dimension / upper**(dimension-1)
    return {
        'dimension': dimension,
        'covariance_lambda_min_lower': lower,
        'covariance_lambda_max_upper': upper,
        'joint_storage_lambda_min_lower': min(1/upper, F(1)),
        'joint_storage_lambda_max_upper': max(1/lower, F(1)),
        'conditional_format_coercivity_lemma_proved': True,
        'requires_every_carried_covariance_finite_symmetric_SPD': True,
        'shipping_SPD_preservation_proved': False,
        'source_uniform_rho_certified': False,
        'storage_search_allowed': False,
        'ALT_STARTUP_PASS': False, 'ALT_LIVE_PASS': False, 'ALT_END_TO_END_PASS': False,
    }


def binary32_value(value):
    value = F(value)
    try:
        point = float(value)
        rounded = struct.unpack('!f', struct.pack('!f', point))[0]
    except (OverflowError, ValueError) as exc:
        raise ValueError('finite binary32 covariance entry required') from exc
    if not math.isfinite(point) or not math.isfinite(rounded) or F(rounded) != value:
        raise ValueError('exact finite binary32 covariance entry required')
    return value


def exact_spd(matrix):
    """Exact unpivoted LDL^T/Sylvester test, with no floating tolerance."""
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError('nonempty square covariance required')
    a = [list(map(binary32_value, row)) for row in matrix]
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('exact covariance symmetry required')
    pivots = []
    for k in range(n):
        pivot = a[k][k]
        if pivot <= 0:
            raise ValueError('covariance is not SPD')
        pivots.append(pivot)
        for i in range(k+1, n):
            if a[i][k] == 0:
                continue
            for j in range(i, n):
                if a[j][k] == 0:
                    continue
                a[i][j] -= a[i][k]*a[j][k]/pivot
                a[j][i] = a[i][j]
    determinant = math.prod(pivots)
    lattice_determinant = determinant / QUANTUM**n
    if lattice_determinant.denominator != 1 or lattice_determinant < 1:
        raise ArithmeticError('dyadic determinant separation failed')
    return {'dimension': n, 'exact_SPD': True, 'determinant': determinant,
            'minimum_exact_LDL_pivot': min(pivots),
            'determinant_grid_separation_verified': True}


def audit_native_endpoints(path):
    rows = []
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for line in stream:
            digest.update(line)
            words = line.decode().split()
            step, kind = map(int, words[:2])
            if (step, kind) == (0, 0) or (kind == 100 and step in (600, 1200)):
                if len(words) != 444:
                    raise ValueError('full 21-state endpoint covariance required')
                # 17-digit native serialization roundtrips each binary32 value
                # through binary64 exactly; binary32_value independently checks.
                values = [F(float(x)) for x in words[3:]]
                result = exact_spd([values[i:i+21] for i in range(0, 441, 21)])
                rows.append({'step': step, 'active_bias': bool(int(words[2])), **result})
    if [r['step'] for r in rows] != [0, 600, 1200]:
        raise ValueError('expected reached root and consecutive native endpoints')
    return {'raw_trace_sha256': digest.hexdigest(), 'endpoints': rows,
            'all_intermediate_covariances_audited': False,
            'infinite_SPD_preservation_proved': False}


def audit_recorded_covariances(path):
    """Check each recorded event-boundary covariance, not internal partial writes."""
    counts, digest, minimum, previous, last_result = {}, hashlib.sha256(), None, None, None
    last_step, sample_boundaries, unique_checks = 0, 0, 0
    kinds = {0, 100, 1002, 1003, 1004, 1011, 1021, 1031, 1061}
    with Path(path).open('rb') as stream:
        for line in stream:
            digest.update(line)
            words = line.decode().split()
            step, kind = map(int, words[:2])
            if kind not in kinds:
                continue
            if len(words) != 444:
                raise ValueError('full covariance record required')
            current = tuple(words[3:])
            if current != previous:
                values = [F(float(x)) for x in current]
                last_result = exact_spd([values[i:i+21] for i in range(0, 441, 21)])
                previous = current
                unique_checks += 1
            minimum = (last_result['minimum_exact_LDL_pivot'] if minimum is None else
                       min(minimum, last_result['minimum_exact_LDL_pivot']))
            counts[str(kind)] = counts.get(str(kind), 0)+1
            if kind == 100:
                if step != last_step+1:
                    raise ValueError('native sample ancestry gap')
                last_step = step
                sample_boundaries += 1
    if counts.get('0') != 1 or last_step != 1200 or sample_boundaries != 1200:
        raise ValueError('complete root-to-1200 trace required')
    return {'raw_trace_sha256': digest.hexdigest(), 'records_by_kind': counts,
            'covariance_records_checked': sum(counts.values()),
            'distinct_consecutive_covariances_factored': unique_checks,
            'minimum_exact_LDL_pivot': minimum, 'sample_boundaries': sample_boundaries,
            'all_recorded_covariances_exact_finite_symmetric_SPD': True,
            'internal_partial_matrix_writes_audited': False,
            'infinite_SPD_preservation_proved': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--audit-retained-traces', action='store_true')
    args = parser.parse_args()
    histories = {}
    for mode in ('H', 'A', 'HA'):
        for stride in (0, 8, 40, 200):
            name = f'{mode}-stride{stride}'
            histories[name] = audit_native_endpoints(args.native_directory/(name+'-observed.txt'))
    report = {'qualification': QUALIFICATION, 'format_bounds': bounds(), 'histories': histories,
              'coercivity_on_finite_SPD_binary32_domain': True,
              'shipping_infinite_domain_membership_proved': False}
    if args.audit_retained_traces:
        report['retained_service_finite_trace_audit'] = {
            mode: audit_recorded_covariances(args.native_directory/(mode+'-stride8-observed.txt'))
            for mode in ('H', 'A', 'HA')}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, default=str, sort_keys=True, indent=2)+'\n')
    print('Exact SPD and determinant-grid checks passed at 36 native endpoints.')


if __name__ == '__main__':
    main()
