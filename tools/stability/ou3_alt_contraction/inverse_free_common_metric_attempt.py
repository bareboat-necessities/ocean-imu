#!/usr/bin/env python3
"""First common-M outward certificate on the corrected inverse-free endpoint family.

The same source-independent coercive 24x24 diagonal candidate is tested against
all H/A x BIAS0/1/2 600-step endpoint outers.  Only outward projected LDLT can
close an endpoint candidate.  Async magnetometer nonexpansiveness remains a
separate same-M obligation, so this module never promotes ALT_LIVE_PASS alone.
"""
from __future__ import annotations
from tools.stability.ou3_alt_contraction import inverse_free_endpoint_outer_attempt as END
from tools.stability.ou3_alt_contraction import first_common_metric_attempt as BASE
from tools.stability.ou3_alt_contraction import projected_storage_certificate as CERT
QUALIFICATION='OU3_ALT_INVERSE_FREE_COMMON_JOINT24_METRIC_ATTEMPT_V1'
RHO=(0.9,0.99,0.999,0.9999,0.99999)
def build():
    ep=END.build();ef=END.validate(ep)
    if ef:raise RuntimeError('inverse-free endpoint prerequisite invalid: '+repr(ef))
    M,scales=BASE.normalized_diagonal_metric();attempts=[];winner=None;finite=ep['all_endpoint_outers_finite']
    if finite:
        for rho in RHO:
            rows={};ok=True
            for mode,modes in ep['reports'].items():
                rows[mode]={}
                for fam,row in modes.items():
                    c=CERT.certify_joint24_bias_supply_family(row['endpoint'],M,rho);rows[mode][fam]=c;ok=ok and c['metric_spd'] and c['projected_strict']
            attempts.append({'rho':rho,'all_pass':bool(ok),'families':rows})
            if ok and winner is None:winner=rho
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','corrected_inverse_free_endpoint_outer_consumed':True,'superseded_endpoint_Pbar_local_route_consumed':False,'same_M_all_modes_bias_families':True,'trajectory_replay_fit_used':False,'metric_scales':scales,'endpoint_outers_finite':finite,'attempts':attempts,'passing_rho':winner,'common_M_IMU_endpoint_projected_LDLT_closed':winner is not None,'async_magnetometer_nonexpansive_same_M_closed':False,'bounded_source_supply_completion_closed':False,'ALT_LIVE_PASS':False,'next_obligation':('prove every admitted async magnetometer event nonexpansive with this exact M and complete bounded physical supplies' if winner is not None else ('search/optimize one common M against this certified endpoint family' if finite else 'refine the inverse-free endpoint partition before common-M search'))}
def summary(d):
    return {'endpoint_outers_finite':d['endpoint_outers_finite'],'passing_rho':d['passing_rho'],'attempts':[{'rho':a['rho'],'all_pass':a['all_pass'],'pivots':{m:{f:r.get('worst_projected_ldlt_pivot_lower') for f,r in rows.items()} for m,rows in a['families'].items()}} for a in d['attempts']]}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('corrected_inverse_free_endpoint_outer_consumed','same_M_all_modes_bias_families'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('superseded_endpoint_Pbar_local_route_consumed','trajectory_replay_fit_used','async_magnetometer_nonexpansive_same_M_closed','bounded_source_supply_completion_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('common_M_IMU_endpoint_projected_LDLT_closed') is not (d.get('passing_rho') is not None):f.append('winner mismatch')
    return f
