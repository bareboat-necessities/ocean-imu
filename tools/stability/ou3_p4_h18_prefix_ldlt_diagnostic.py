#!/usr/bin/env python3
"""Actual H18 moving-information every-prefix LDLT diagnostic.

This is deliberately NOT the production P4 certificate.  It uses the exact
synchronized H18 Phi transport and actual event-local covariance history, but
collapses the entire transported finite defect d_l to one conservative source
coordinate w.  That loses COMPLETE-BRMM moment/radial correlation and therefore
cannot promote P4.  Its purpose is to answer a narrower question with the real
matrices: does the homogeneous moving-information prefix margin survive, and
what source-supply coefficient would make the augmented interval LDLT close?

For each prefix l,

    D_l = J_0 - M_l^T J_l M_l,
    d_l = B_l w,   B_l := interval enclosure of direct_defect_l,

and the exact quadratic master is constructed by the canonical joint-sector
master.  A one-dimensional ISS source selector is then penalized by gamma_s and
tested through the same outward LDLT backend used by production.
"""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_complete_brmm_joint_sector_master as MASTER
import ou3_p4_joint_iss_augmented_master as ISS
import ou3_p4_h18_source_indexed_prefix_transport as HTR

SCHEMA=1
QUALIFICATION='OU3_P4_H18_ACTUAL_PREFIX_LDLT_DIAGNOSTIC_V1'
P3_DELTA=1.0e-18

def I(x):return Interval.point(float(x))
def zeros(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def col(v):return [[x] for x in v]
def covariance_nodes(samples):
    nodes=[samples[0].cells[0].P]
    for si,s in enumerate(samples):
        if si:nodes.append(nodes[-1]) # proof-coordinate rebase: covariance unchanged
        trusted={c.event_index_in_sample:c for c in s.selector.H_event_cells}
        for cell in s.cells:nodes.append(trusted[cell.event_ordinal].P_after)
    return nodes
def defect_matrix(J0,JN,M):
    return matrix_symmetric_hull(matrix_sub(J0,matrix_mul(matrix_mul(matrix_transpose(M),JN),M)))
def augmented_for_prefix(J0,JN,tr):
    M=tr['M_word'];D=defect_matrix(J0,JN,M);B=col(tr['direct_defect'])
    L=MASTER.master_quadratic_matrix(D,M,JN,B)
    return D,L
def selector(n,index):
    A=[[I(0) for _ in range(n)]];A[0][index]=I(1);return A
def zero_sector(n):return zeros(n,n)
def sweep_gamma(L):
    n=len(L);S=selector(n,n-1);Z=zero_sector(n)
    # Broad deterministic sweep.  This is a diagnostic coefficient, not a
    # theorem constant and not used by production.
    candidates=[0.0]+[10.0**p for p in range(-18,37)]
    last=[]
    for g in candidates:
        ok,p=ISS.certify_iss(L,[Z],[0.0],S,g,None,0.0);last=p
        if ok:return True,g,p
    return False,None,last
def audit(samples):
    t=HTR.build_transport(samples);P=covariance_nodes(samples)
    if len(P)!=len(t['events'])+1:raise RuntimeError('covariance/transport boundary count mismatch')
    J0=matrix_inverse_gauss_jordan(P[0]);rows=[]
    for i,r in enumerate(t['decompositions']):
        JN=matrix_inverse_gauss_jordan(P[i+1]);D,L=augmented_for_prefix(J0,JN,r['transport']);ok,g,p=sweep_gamma(L)
        rows.append({'prefix':i+1,'kind':t['events'][i]['kind'],'ldlt_closes_with_scalar_source_supply':ok,'first_gamma':g,'pivot_lowers':p,'D_diagonal_lower_min':min(float(D[k][k].lo) for k in range(18)),'identity_residual_contains_zero':r['transport']['identity_residual_contains_zero']})
    return rows
def build():
    h=HTR.build();hf=HTR.validate(h);m=MASTER.build();mf=MASTER.validate(m);iss=ISS.build();ifail=ISS.validate(iss)
    if hf or mf or ifail:raise RuntimeError(f'diagnostic prerequisites failed H18={hf} master={mf} ISS={ifail}')
    rows=audit(HTR.smoke_objects());closed=bool(rows) and all(r['ldlt_closes_with_scalar_source_supply'] for r in rows)
    fail=next((r for r in rows if not r['ldlt_closes_with_scalar_source_supply']),None)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,'actual_synchronized_H18_prefix_transport_consumed':True,'actual_event_covariance_information_matrices_consumed':True,'full_augmented_interval_LDLT_backend_consumed':True,'one_conservative_independent_source_column_used_for_diagnostic':True,'source_moment_radial_correlation_preserved_in_diagnostic':False,'binary32_channel_included_in_diagnostic':False,'diagnostic_all_prefixes_close':closed,'first_failed_prefix':fail,'prefixes':rows,'production_endpoint_LDLT_promoted_here':False,'production_every_prefix_LDLT_promoted_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'replace the diagnostic defect column by the common COMPLETE-BRMM moment/radial/source graph map and add conditional binary32; retain the same actual D_l,M_l,J_l matrices and rerun every-prefix LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('actual_synchronized_H18_prefix_transport_consumed','actual_event_covariance_information_matrices_consumed','full_augmented_interval_LDLT_backend_consumed','one_conservative_independent_source_column_used_for_diagnostic'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_moment_radial_correlation_preserved_in_diagnostic','binary32_channel_included_in_diagnostic','production_endpoint_LDLT_promoted_here','production_every_prefix_LDLT_promoted_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not d.get('prefixes'):f.append('empty diagnostic prefix set')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'all_close':d['diagnostic_all_prefixes_close'],'first_failed':d['first_failed_prefix'],'prefix_count':len(d['prefixes']),'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
