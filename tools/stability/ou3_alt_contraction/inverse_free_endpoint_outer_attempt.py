#!/usr/bin/env python3
"""First 600-step endpoint outer built from the valid inverse-free local family.

For each H/A x BIAS family, the two IMU-core maps are first hulled and then
raised to the 600-step product set with outward interval multiplication.  This
is a conservative representation of all S-due/not-due schedules.  It is not a
source generator.  Async magnetometer events remain separate.
"""
from __future__ import annotations
import math
from ou3_interval import matrix_identity,matrix_mul
from tools.stability.ou3_alt_contraction import inverse_free_universal_local_outer as LOCAL
from tools.stability.ou3_alt_contraction import endpoint_family_induction as OUTER
QUALIFICATION='OU3_ALT_INVERSE_FREE_600_STEP_ENDPOINT_OUTER_V1';N=600

def power(A,n=N):
    r=matrix_identity(len(A));b=A;k=int(n)
    while k:
        if k&1:r=matrix_mul(b,r)
        k//=2
        if k:b=matrix_mul(b,b)
    return r
def finite(A):return all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in A for x in row)
def width(A):return max((x.hi-x.lo for row in A for x in row),default=0.0)
def magnitude(A):return max((max(abs(x.lo),abs(x.hi)) for row in A for x in row),default=0.0)
def build():
    local=LOCAL.build();lf=LOCAL.validate(local)
    if lf:raise RuntimeError('local inverse-free outer invalid: '+repr(lf))
    reports={}
    for mode,rows in local['reports'].items():
        reports[mode]={}
        for fam,row in rows.items():
            out={'local_core_finite':row['imu_core_family_finite'],'endpoint':None,'endpoint_finite':False,'max_width':None,'max_abs':None,'failure':None}
            if row['imu_core_family_finite']:
                try:
                    H=OUTER.hull_matrices(row['imu_core_family']);E=power(H);out.update(endpoint=E,endpoint_finite=finite(E),max_width=width(E),max_abs=magnitude(E))
                    if not out['endpoint_finite']:out['failure']='interval product became nonfinite'
                except Exception as exc:out['failure']=type(exc).__name__+': '+str(exc)
            else:out['failure']='valid inverse-free local IMU core was not finite'
            reports[mode][fam]=out
    closed=all(r['endpoint_finite'] for rows in reports.values() for r in rows.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','inverse_free_local_family_consumed':True,'endpoint_Pbar_route_consumed':False,'word_transitions':N,'outward_hull_power_is_outer_only':True,'reports':reports,'all_endpoint_outers_finite':closed,'numeric_600_step_IMU_outer_materialized':closed,'async_magnetometer_star_still_separate':True,'common_storage_closed':False,'ALT_LIVE_PASS':False,'next_obligation':('test one common coercive M with outward projected LDLT on every endpoint outer and separately on async magnetic nonexpansiveness' if closed else 'partition the valid inverse-free local outer before multiplication; do not return to endpoint Pbar or replay')}
def summary(d):return {m:{f:{k:r[k] for k in ('local_core_finite','endpoint_finite','max_width','max_abs','failure')} for f,r in rows.items()} for m,rows in d['reports'].items()}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('word_transitions')!=N:f.append('qualification/window mismatch')
    for k in ('inverse_free_local_family_consumed','outward_hull_power_is_outer_only','async_magnetometer_star_still_separate'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('endpoint_Pbar_route_consumed','common_storage_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
