#!/usr/bin/env python3
"""Non-promoting ALT experiments: exact graph checks and frozen-map diagnostics.

Neither a captured trajectory nor high precision qualifies a source family.
This runner never searches/fits a metric to replay and never promotes a gate.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT / 'tools/stability'),
                str(ROOT / 'tests/kalman_ou_iii')]
from tools.stability.ou3_alt_contraction.core import open_obligations, structural_findings, masked_supply_certificate
import ou3_p4_complete_brmm_word_feasibility as old


def high_precision_ratio(A, P0, P1, digits: int) -> dict:
    """Reevaluate a STORED binary32 product; not high-precision event composition."""
    with mp.workdps(digits):
        a = mp.matrix(A.tolist())
        p0, p1 = mp.matrix(P0.tolist()), mp.matrix(P1.tolist())
        p0, p1 = (p0+p0.T)/2, (p1+p1.T)/2
        L = mp.cholesky(p0)
        T = a*L
        U = mp.matrix(a.rows, a.cols)
        for j in range(a.cols):
            col = mp.lu_solve(p1, T[:, j])
            for i in range(a.rows):
                U[i, j] = col[i]
        B = T.T*U
        B = (B+B.T)/2
        values, vectors = mp.eigsy(B)
        rho, v = values[values.rows-1, 0], vectors[:, vectors.cols-1]
        residual = mp.norm(B*v-rho*v)
        return {'digits': digits, 'rho': mp.nstr(rho, digits-10),
                'distance_to_one': mp.nstr(1-rho, digits-10),
                'eigen_residual_norm': mp.nstr(residual, 12)}


def bias_contracts() -> dict:
    families = {}
    keys = ('phi_true_interval', 'true_bias_norm_upper_mps2',
            'driver_increment_norm_upper_mps2')
    for name in ('BIAS0', 'BIAS1', 'BIAS2'):
        mod = importlib.import_module('ou3_p4_'+name.lower()+'_family')
        definition = mod.build()
        failures = mod.validate(definition)
        families[name] = {k: definition[k] for k in keys}
        families[name].update(definition_validation_pass=not failures,
                              definition_validation_failures=failures,
                              hardware_qualification=False)
    aggregate = importlib.import_module('ou3_p4_bias_family_joint_iss_supply')
    try:
        result = aggregate.build()
    except RuntimeError as exc:
        if not str(exc).startswith('qualified fresh-entry source failed:'):
            raise
        admission = {'builder_completed': False, 'error': str(exc),
                     'classification': 'shared fresh-Live hard-entry qualification',
                     'ALT_entry_certified': False}
    else:
        # Successful construction alone would still not certify ALT entry.
        admission = {'builder_completed': True, 'ALT_entry_certified': False,
                     'qualification': result.get('qualification')}
    return {'families': families, 'aggregate': admission}


def analyze_traces(map_path: Path, cov_path: Path) -> dict:
    ms, maps = old.read_maps(map_path)
    cs, covs = old.read_covariances(cov_path)
    if ms != cs or len(maps) != len(covs):
        raise ValueError('map/covariance traces must have equal stride AND record count')
    baseline = old.analyze(map_path, cov_path)
    modes = {}
    for mode, dim in (('H18', 18), ('A21', 21)):
        row = baseline['modes'][mode]
        if not row['legal_blocks']:
            raise ValueError('no retained '+mode+' diagnostic word')
        worst = row['worst']
        i = worst['record']
        a, p0, p1 = maps[i]['M'][:dim, :dim], covs[i]['P0'][:dim, :dim], covs[i]['P1'][:dim, :dim]
        hp = [high_precision_ratio(a, p0, p1, d) for d in (80, 120)]
        with mp.workdps(120):
            agreement = abs(mp.mpf(hp[0]['rho'])-mp.mpf(hp[1]['rho']))
            agreement_text = mp.nstr(agreement, 12)
        modes[mode] = {'retained_words': row['legal_blocks'], 'worst': worst,
                       'high_precision_stored_map': hp,
                       'precision_comparison_abs_difference': agreement_text}
    # The existing H18 observer's homogeneous map leaves the held coordinate
    # unchanged. Its reconstructed H=P^-1*N is not the full residual Jacobian;
    # do NOT label this a full physical held-bias/true-bias experiment.
    i = modes['H18']['worst']['record']
    A = maps[i]['M']
    modes['H18']['observer_held_block'] = {
        'lower_left_norm': float(np.linalg.norm(A[18:21, :18])),
        'lower_right_minus_identity_norm': float(np.linalg.norm(A[18:21, 18:21]-np.eye(3))),
        'upper_right_norm': float(np.linalg.norm(A[:18, 18:21])),
        'physical_held_bias_forcing_bound': False,
    }
    return {'qualification': 'CAPTURED_FROZEN_COEFFICIENT_DIAGNOSTIC_ONLY',
            'metric_fitted_to_replay': False,
            'full_joint24_word_feasibility_established': False,
            'event_product_recomputed_at_high_precision': False,
            'endogenous_gain_derivatives_included': False,
            'H18_held_error_zero_in_18_state_ratio': True,
            'A21_full_21_state_covariance_kept': True,
            'source_uniform': False, 'nonlinear_basin_certified': False,
            'modes': modes, 'rejected_counts': baseline['records_rejected'],
            'limitation': 'OU3MAP3 uses covariance-reconstructed H and a binary32 accumulated product; actual residual/held-bias forcing, source/tuner dependence, guards and joint24 truth are not certified by this trace.'}


def run(map_path=None, cov_path=None) -> dict:
    if (map_path is None) != (cov_path is None):
        raise ValueError('--maps and --cov must be supplied together')
    files = ['src/kalman_ou_iii/Kalman3D_Wave_OU_III.h',
             'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h',
             'tools/stability/ou3_brmm_physical_wave_source.py',
             'tools/stability/ou3_p4_closure_domain.json',
             'tests/kalman_ou_iii/ou3-operation-ledger-sim.cpp']
    report = {'qualification': 'ALT_ARCHITECTURE_AND_EXPLORATORY_ALGEBRA_ONLY',
              'status': open_obligations(), 'exact_findings': structural_findings(),
              'exact_masked_supply_analogue': masked_supply_certificate(),
              'bias_contracts': bias_contracts(),
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
    if map_path is not None:
        report['trace_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (map_path, cov_path)}
        report['frozen_map_diagnostic'] = analyze_traces(map_path, cov_path)
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--maps', type=Path)
    ap.add_argument('--cov', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    report = run(args.maps, args.cov)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'qualification': report['qualification'], 'status': report['status'],
                      'aggregate_entry': report['bias_contracts']['aggregate']}, indent=2))
    # Successful experiment execution, not mathematical PASS.
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
