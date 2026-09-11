#!/usr/bin/env python3
"""First 600-step numerical outer for the ALT Normal-Live IMU core.

For each H/A x BIAS0/1/2 family, ``coarse_imu_core_outer`` supplies two outward
24x24 per-sample maps: S not-due and S due.  Their entrywise interval hull H is
an outer set containing either literal map.  Therefore every length-N product
A_N...A_1 with A_k chosen from those maps belongs to the interval matrix power
H^N computed with outward interval multiplication.

This intentionally forgets scheduler correlations after proving set inclusion;
it does NOT assert the hull is a realizable source transition.  Fast binary
exponentiation is only an exact arithmetic optimization of that outer product
set.  If dependency blow-up makes H^600 useless/nonfinite, the result says the
single-hull representation is too coarse and must be partitioned; it does not
shrink the physical theorem family.

Asynchronous magnetometer events remain excluded and are governed separately by
``asynchronous_event_star``.
"""
from __future__ import annotations

import math
from ou3_interval import matrix_identity,matrix_mul
from tools.stability.ou3_alt_contraction import coarse_imu_core_outer as CORE
from tools.stability.ou3_alt_contraction import endpoint_family_induction as OUTER
from tools.stability.ou3_alt_contraction import asynchronous_event_star as MAGSTAR

QUALIFICATION='OU3_ALT_COARSE_600_STEP_IMU_ENDPOINT_OUTER_V1'
TRANSITIONS=600


def interval_matrix_power(A,n:int):
    n=int(n)
    if n<0:raise ValueError('nonnegative power required')
    dim=len(A)
    if dim==0 or any(len(r)!=dim for r in A):raise ValueError('square matrix required')
    result=matrix_identity(dim);base=A
    while n:
        if n&1:result=matrix_mul(base,result)
        n//=2
        if n:base=matrix_mul(base,base)
    return result

def finite_matrix(A):return all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in A for x in row)
def max_width(A):return max((x.hi-x.lo for row in A for x in row),default=0.0)
def max_abs(A):return max((max(abs(x.lo),abs(x.hi)) for row in A for x in row),default=0.0)

def build():
    core=CORE.build();cf=CORE.validate(core);star=MAGSTAR.build();sf=MAGSTAR.validate(star)
    if cf or sf:raise RuntimeError(f'endpoint outer prerequisites failed core={cf} star={sf}')
    reports={}
    for mode,rows in core['reports'].items():
        reports[mode]={}
        for family,row in rows.items():
            r={'imu_core_family_finite':row['imu_core_family_finite'],'endpoint_outer_finite':False,'endpoint_outer':None,'max_entry_width':None,'max_entry_abs':None,'failure':None}
            if row['imu_core_family_finite']:
                try:
                    H=OUTER.hull_matrices(row['imu_core_family']);E=interval_matrix_power(H,TRANSITIONS)
                    r.update(endpoint_outer=E,endpoint_outer_finite=finite_matrix(E),max_entry_width=max_width(E),max_entry_abs=max_abs(E))
                    if not r['endpoint_outer_finite']:r['failure']='600-step interval hull power became nonfinite'
                except Exception as exc:r['failure']=type(exc).__name__+': '+str(exc)
            else:r['failure']='local IMU-core outer not finite; endpoint power not attempted'
            reports[mode][family]=r
    all_finite=all(r['endpoint_outer_finite'] for rows in reports.values() for r in rows.values())
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','word_transitions':TRANSITIONS,
      'S_due_two_pattern_hull_is_outer_only_not_history_generator':True,
      'binary_interval_power_encloses_all_length_600_IMU_core_products':True,
      'scheduler_correlations_relaxed_only_after_local_same_cell_outer':True,
      'async_magnetometer_excluded_and_star_obligation_retained':True,
      'reports':reports,'all_600_step_IMU_endpoint_outers_finite':all_finite,
      'single_hull_600_step_representation_closed':all_finite,
      'async_magnetometer_nonexpansive_same_M_closed':False,
      'common_storage_closed':False,'ALT_LIVE_PASS':False,
      'next_obligation':('use these endpoint outers only as the first common-M candidate family and certify projected storage; separately certify source-uniform magnetometer nonexpansiveness with the same M' if all_finite else 'partition the local source/scheduler/state/covariance outer before multiplication; do not shrink COMPLETE-BRMM')
    }
def summary(d):
    return {m:{f:{k:r[k] for k in ('imu_core_family_finite','endpoint_outer_finite','max_entry_width','max_entry_abs','failure')} for f,r in rows.items()} for m,rows in d['reports'].items()}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('word_transitions')!=TRANSITIONS:f.append('qualification/window mismatch')
    for k in ('S_due_two_pattern_hull_is_outer_only_not_history_generator','binary_interval_power_encloses_all_length_600_IMU_core_products','scheduler_correlations_relaxed_only_after_local_same_cell_outer','async_magnetometer_excluded_and_star_obligation_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('async_magnetometer_nonexpansive_same_M_closed','common_storage_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
