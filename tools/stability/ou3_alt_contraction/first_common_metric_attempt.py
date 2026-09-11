#!/usr/bin/env python3
"""Guarded entry to the first ALT common-joint24 storage calculation.

The retained candidate calculation is not allowed to run on the current
Jacobian cocycle. A finite physical descriptor with reference forcing and all
branch/product graphs is a prerequisite. The assembly flags alone do not
establish that prerequisite. See docs/ou3-alt-finite-measurement-proof.md.

If a valid complete finite master is attached, the same source-independent M
must be checked for every H18/A21 x BIAS family and asynchronous event branch.
"""
from __future__ import annotations

import math
import numpy as np
from tools.stability.ou3_alt_contraction import proof_plan as PLAN
from tools.stability.ou3_alt_contraction import physical_word as WORD

QUALIFICATION='OU3_ALT_FIRST_COMMON_JOINT24_METRIC_ATTEMPT_V1'
RHO_CANDIDATES=(0.9,0.99,0.999,0.9999,0.99999)


def normalized_diagonal_metric():
    PLAN.assert_finite_storage_master(WORD.finite_storage_readiness())
    # Imports remain behind the finite-word guard in build(); a hard-entry
    # builder and a covariance outer cannot repair a missing finite identity.
    import ou3_p4_hard_entry_set as ENTRY
    from tools.stability.ou3_alt_contraction import bias_families as BIAS
    e=ENTRY.build();f=ENTRY.validate(e)
    if f:raise RuntimeError('hard regional scale prerequisite failed: '+repr(f))
    r=e['coordinate_radii']
    scales=[]
    for name in ('attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2'):
        scales.extend([float(r[name])]*3)
    contracts=BIAS.contracts()
    beta=max(c.true_bias_norm_bound for c in contracts)
    eba=max(float(r['accelerometer_bias_error_norm_mps2']), beta+float(r['accelerometer_bias_error_norm_mps2']))
    scales.extend([eba]*3);scales.extend([beta]*3)
    if len(scales)!=24 or any(not math.isfinite(x) or x<=0 for x in scales):raise RuntimeError('invalid normalization scales')
    diag=[1.0/(x*x) for x in scales]
    M=np.diag(diag)
    return M,scales

def build():
    PLAN.assert_finite_storage_master(WORD.finite_storage_readiness())
    from tools.stability.ou3_alt_contraction import coarse_endpoint_outer_attempt as ENDPOINT
    from tools.stability.ou3_alt_contraction import projected_storage_certificate as CERT
    ep=ENDPOINT.build();ef=ENDPOINT.validate(ep)
    if ef:raise RuntimeError('endpoint outer prerequisite failed: '+repr(ef))
    M,scales=normalized_diagonal_metric();attempts=[];winner=None
    endpoints_finite=bool(ep['all_600_step_IMU_endpoint_outers_finite'])
    if endpoints_finite:
        for rho in RHO_CANDIDATES:
            rows={};all_pass=True
            for mode,modes in ep['reports'].items():
                rows[mode]={}
                for family,row in modes.items():
                    cert=CERT.certify_joint24_bias_supply_family(row['endpoint_outer'],M,rho)
                    rows[mode][family]=cert
                    all_pass=all_pass and bool(cert['metric_spd'] and cert['projected_strict'])
            attempts.append({'rho':rho,'all_families_projected_strict':all_pass,'families':rows})
            if all_pass and winner is None:winner=rho
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'candidate_metric_source':'declared regional coordinate normalization only',
      'trajectory_or_replay_fit_used_for_M':False,'per_mode_metric_used':False,'per_bias_family_metric_used':False,
      'metric_dimension':24,'metric_coercive_point_candidate':True,'normalization_scales':scales,
      'same_M_used_for_all_H_A_BIAS_families':True,'rho_candidates':list(RHO_CANDIDATES),
      'endpoint_outer_families_finite':endpoints_finite,'attempts':attempts,'first_passing_rho':winner,
      'common_M_projected_IMU_endpoint_LDLT_closed':winner is not None,
      'async_magnetometer_nonexpansive_same_M_closed':False,
      'common_joint24_storage_complete':False,'ALT_LIVE_PASS':False,
      'next_obligation':('certify every admitted asynchronous magnetometer event nonexpansive with this same M, then add bounded physical source/bias supplies' if winner is not None else ('optimize one common coercive M against the certified endpoint outers; every candidate must still pass outward LDLT on all six families' if endpoints_finite else 'refine the coarse universal endpoint representation before metric search'))
    }
def summary(d):
    out={'endpoint_outer_families_finite':d['endpoint_outer_families_finite'],'first_passing_rho':d['first_passing_rho'],'attempts':[]}
    for a in d['attempts']:
        piv={m:{f:r.get('worst_projected_ldlt_pivot_lower') for f,r in rows.items()} for m,rows in a['families'].items()}
        out['attempts'].append({'rho':a['rho'],'all_pass':a['all_families_projected_strict'],'worst_pivots':piv})
    return out
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('metric_dimension')!=24:f.append('qualification/dimension mismatch')
    for k in ('metric_coercive_point_candidate','same_M_used_for_all_H_A_BIAS_families'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_or_replay_fit_used_for_M','per_mode_metric_used','per_bias_family_metric_used','async_magnetometer_nonexpansive_same_M_closed','common_joint24_storage_complete','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('common_M_projected_IMU_endpoint_LDLT_closed') is not (d.get('first_passing_rho') is not None):f.append('winner flag mismatch')
    return f
