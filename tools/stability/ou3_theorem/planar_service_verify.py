"""Fail-closed verification of finite planar evidence, never theorem admission.

The numerical limits below are regression tolerances for this native profile,
not all-time arithmetic bounds. A successful check verifies only the supplied
finite word and its source bindings. Missing fields and promoted status fail.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .planar_service_stream import HEADER, PROBE, instrument


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _finite(value: Any, name: str, lower: float, upper: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name}: finite numeric value required')
    if not math.isfinite(value) or not lower <= value <= upper:
        raise ValueError(f'{name}: outside finite diagnostic range [{lower}, {upper}]')
    return float(value)


def verify(stream: dict, audit: dict, frechet: dict) -> dict:
    """Validate current-source evidence and return a non-promoting summary."""
    expected_qualifications = (
        'OU3_PLANAR_SERVICE_OPERATION_STREAM_V1',
        'OU3_PLANAR_SERVICE_OPERATION_AUDIT_V1',
        'OU3_PLANAR_COVARIANCE_FRECHET_DIAGNOSTIC_V1')
    for report, qualification in zip((stream, audit, frechet), expected_qualifications):
        if report.get('qualification') != qualification:
            raise ValueError('wrong diagnostic qualification')
        if report.get('result_type') != 'FINITE DIAGNOSTIC ONLY':
            raise ValueError('finite evidence cannot be promoted')
        for key in ('all_time_magnetic_service_verified', 'theorem_closed'):
            if report.get(key) is not False:
                raise ValueError(f'{key}: must remain explicitly false')
    forbidden_promotions = ((stream, ('shipping_behavior_changed', 'full_shipping_counterexample_admitted')),
        (audit, ('scheduler_phase_cell_forward_invariant', 'local_commutator_is_complete_scheduler_word_defect')),
        (frechet, ('forward_invariant_cell_certified', 'rounding_certified',
                   'scheduler_branch_coverage_certified', 'generated_mean_coefficient_variations_included')))
    for report, keys in forbidden_promotions:
        for key in keys:
            if report.get(key) is not False:
                raise ValueError(f'{key}: unverified promotion or missing qualification')
    if stream.get('untapped_complete_tail_samples_bitwise_equal') is not True:
        raise ValueError('untapped shipping control is required')
    if frechet.get('all_AW_sync_derivatives_retained') is not True:
        raise ValueError('AW derivative missing')
    sha = stream.get('stream_sha256')
    if not isinstance(sha, str) or len(sha) != 64 or any(c not in '0123456789abcdef' for c in sha):
        raise ValueError('invalid stream digest')
    if audit.get('stream_sha256') != sha or frechet.get('stream_sha256') != sha:
        raise ValueError('reports do not belong to the same literal stream')
    bound = {
        'probe_sha256': _digest(PROBE),
        'driver_sha256': _digest(PROBE.with_name('planar_service_stream.py')),
        'shipping_header_sha256': _digest(HEADER),
        'instrumented_header_sha256': hashlib.sha256(instrument(HEADER.read_text()).encode()).hexdigest(),
    }
    for key, value in bound.items():
        if stream.get(key) != value:
            raise ValueError(f'source binding changed: {key}')
    counts = stream['record_counts']
    for kind, count in counts.items():
        _finite(count, f'count {kind}', 1, 10**8)
        if audit['counts'].get(kind) != count:
            raise ValueError('operation coverage mismatch')
    expected_kinds=set(map(str,range(1,10)))
    if set(counts) != expected_kinds:
        raise ValueError('missing operation kind')
    if not (counts['1'] == counts['2'] == counts['6'] == counts['7']):
        raise ValueError('sample/prediction/accelerometer/AW coverage mismatch')
    if '9' in counts and counts['9'] != sum(counts[k] for k in ('2','3','4')):
        raise ValueError('mean tangent/correction coverage mismatch')
    if counts['5'] != counts['8'] or counts['5'] != sum(counts[k] for k in ('2', '3', '4')):
        raise ValueError('correction/reset coverage mismatch')
    first, last = audit['sample_range']
    if last-first+1 != counts['1']:
        raise ValueError('noncontiguous sample range')
    if not first-1 <= frechet['root_sample'] < frechet['end_sample'] <= last:
        raise ValueError('Frechet word lies outside the audited stream')
    if frechet['retained_operations']['1'] != frechet['end_sample']-frechet['root_sample']:
        raise ValueError('incomplete Frechet word')
    native = stream['native']
    if native['parity_max_abs'] != 0:
        raise ValueError('lossless parity compression failed')
    _finite(native['service_min'], 'finite service floor', 1, 1000)
    _finite(native['service_windows'], 'finite service windows', 1, 10**7)
    _finite(audit['smallest_covariance_eigenvalue_observed'], 'observed covariance floor', 1e-12, 1)
    _finite(audit['smallest_innovation_eigenvalue_observed'], 'observed innovation floor', 1e-8, 1e6)
    residuals = audit['maximum_absolute_operation_residuals']
    limits = {'prediction':1e-7, 'innovation':1e-6, 'PHt':1e-7, 'gain_solve':1e-7,
        'joseph':1e-7, 'reset_congruence':1e-8, 'aw_sync':1e-7,
        'next_prediction_input':1e-11, 'correction_input':1e-11, 'sample_output':1e-11,
        'linked_even_identity':1e-12, 'linked_odd_identity':1e-12}
    for key, cap in limits.items():
        _finite(residuals[key], key, 0, cap)
    gains = {}
    for name, dimension in (('even', 78), ('odd', 45)):
        row = frechet['parity_blocks'][name]
        if row['symmetric_dimension'] != dimension or row['same_root_end_covariance_used'] is not False:
            raise ValueError('invalid covariance-map coordinates')
        # A measured gain >=1 fails this feasibility candidate, not the filter.
        gains[name] = _finite(row['relative_Frobenius_gain'], name+' partial gain', 0, 1-1e-10)
    return {'qualification':'OU3_PLANAR_FINITE_EVIDENCE_VERIFICATION_V1',
        'finite_reproduction_pass':True, 'result_type':'FINITE DIAGNOSTIC ONLY',
        'stream_sha256':sha, 'finite_service_floor':native['service_min'],
        'covariance_partial_gains':gains, 'complete_closed_loop_contraction_claimed':False,
        'all_time_magnetic_service_verified':False, 'theorem_closed':False}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stream-report',type=Path,required=True)
    parser.add_argument('--operation-audit',type=Path,required=True)
    parser.add_argument('--frechet',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        result=verify(*(json.loads(p.read_text()) for p in
            (args.stream_report,args.operation_audit,args.frechet)))
    except (OSError, ValueError, KeyError, TypeError) as error:
        result={'finite_reproduction_pass':False,'failure':str(error),
            'all_time_magnetic_service_verified':False,'theorem_closed':False}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')
    return 0 if result['finite_reproduction_pass'] else 1

if __name__=='__main__':
    raise SystemExit(main())
