#!/usr/bin/env python3
"""H18 post-prediction homogeneous endpoint LDLT diagnostic.

This diagnostic uses the *actual* synchronized source-indexed H18 structural
transport on the canonical post-prediction word and the actual shipping
covariance at its two boundaries.  It asks whether the homogeneous map alone
satisfies the P3-scale information-storage inequality

    J0 - M^T JN M - delta J0 >= 0,
    delta = 1e-18,

where J=P^{-1}.  The reset tangent G(I-KH), source-indexed Phi rebases, and
prediction finite-map linear part are therefore present in M.  No source defect,
nonlinear graph sector, bias input or binary32 supply is included.

A pass is only a feasibility signal for the homogeneous block.  A failure is not
an instability result: the production augmented master may recover the exact
signed cancellations using its nonlinear sectors.  This module never promotes
P4.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_h18_post_prediction_word as WORD
import ou3_p4_h18_source_indexed_prefix_transport as HTR

SCHEMA=1
QUALIFICATION='OU3_P4_H18_POST_PREDICTION_HOMOGENEOUS_ENDPOINT_LDLT_DIAGNOSTIC_V1'
DELTA=Interval.point(1e-18)

def scale(A,s):return [[x*s for x in row] for row in A]
def endpoint_matrix(w):
    tr=w['decompositions'][-1]['transport'];M=tr['M_word'];J0=matrix_inverse_gauss_jordan(w['entrance_P']);JN=matrix_inverse_gauss_jordan(w['terminal_P'])
    raw=matrix_symmetric_hull(matrix_sub(J0,matrix_mul(matrix_mul(matrix_transpose(M),JN),M)))
    target=matrix_symmetric_hull(matrix_sub(raw,scale(J0,DELTA)))
    ok,piv=symmetric_positive_definite_ldlt(target)
    return {'closed':ok,'pivots':piv,'raw_diagonal_lower_min':min(x[i].lo for i,x in enumerate(raw)),'target_diagonal_lower_min':min(x[i].lo for i,x in enumerate(target)),'endpoint_kind':w['events'][-1]['kind'],'event_count':len(w['events'])}
def build():
    base=WORD.build();bf=WORD.validate(base)
    if bf:raise RuntimeError('post-prediction word prerequisite failed: '+repr(bf))
    w=WORD.build_word(HTR.smoke_objects());r=endpoint_matrix(w)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':1e-18,
      'actual_post_prediction_storage_boundaries_consumed':True,'actual_source_indexed_H18_transport_consumed':True,
      'finite_reset_tangent_present_in_homogeneous_map':True,'source_coordinate_rebases_present_in_homogeneous_map':True,
      'endpoint_is_following_prediction':r['endpoint_kind']=='prediction','homogeneous_endpoint_LDLT_closes_on_smoke':r['closed'],
      'endpoint_result':r,'source_defect_included':False,'nonlinear_graph_sectors_included':False,'binary32_supply_included':False,
      'production_endpoint_LDLT_promoted_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'retain this actual endpoint D/M/J block, add the common same-history nonlinear graph and correlated BRMM/binary32 supplies, then run the production augmented LDLT; do not infer source-uniform closure from this smoke diagnostic'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=1e-18:f.append('P3 delta changed')
    for k in ('actual_post_prediction_storage_boundaries_consumed','actual_source_indexed_H18_transport_consumed','finite_reset_tangent_present_in_homogeneous_map','source_coordinate_rebases_present_in_homogeneous_map','endpoint_is_following_prediction'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_defect_included','nonlinear_graph_sectors_included','binary32_supply_included','production_endpoint_LDLT_promoted_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not d.get('endpoint_result'):f.append('endpoint result missing')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['homogeneous_endpoint_LDLT_closes_on_smoke'],'result':d['endpoint_result'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
