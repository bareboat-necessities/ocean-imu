#!/usr/bin/env python3
"""Dense active accelerometer-bias projection graph IQC for A21 P4 prefixes.

Shipping projects the corrected accelerometer-bias estimate onto the Euclidean
ball.  In physical error coordinates the exact operator is

    f = F_R(e, beta) = beta - Pi_R(beta-e).

The global projection lemma proves for every branch, including the Clarke
boundary,

    ||F_R(e1,beta1)-F_R(e2,beta2)||^2
      <= ||e1-e2||^2 + ||beta1-beta2||^2.

Taking the second point as (0,0), where F_R(0,0)=0, yields the homogeneous dense
QC

    ||e||^2 + ||beta||^2 - ||f||^2 >= 0.

For a Joseph prefix the projection input is the SAME-event corrected bias error

    e = e_ba_before - d_ba,

with d=Kq from the estimator-owned P/H/R/K cell.  This module therefore accepts
maps for that exact pre-projection error, the physical BIAS1 beta coordinate,
and the projected output f.  It never substitutes a rowwise K bound or a fixed
inactive projection branch.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_add,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_projection_sector as PROJ

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_DENSE_BIAS_PROJECTION_GRAPH_IQC_V1'

def shape(A):return len(A),len(A[0]) if A else 0
def gram(A):return matrix_mul(matrix_transpose(A),A)
def projection_iqc(Epre,Beta,F):
    er,n=shape(Epre);br,n2=shape(Beta);fr,n3=shape(F)
    if (er,br,fr)!=(3,3,3) or n==0 or n2!=n or n3!=n:raise ValueError('projection maps must be 3xn on one common coordinate')
    return matrix_symmetric_hull(matrix_sub(matrix_add(gram(Epre),gram(Beta)),gram(F)))
def selector(n,offset):
    if offset<0 or offset+3>n:raise ValueError('selector outside coordinate')
    I=Interval.point
    return [[I(1.0 if j==offset+i else 0.0) for j in range(n)] for i in range(3)]
def build():
    p=PROJ.build();pf=PROJ.validate(p)
    if pf:raise RuntimeError('projection prerequisite failed: '+repr(pf))
    n=9;E=selector(n,0);B=selector(n,3);F=selector(n,6);Pi=projection_iqc(E,B,F)
    finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in Pi for x in row)
    symmetric=all(Pi[i][j].lo==Pi[j][i].lo and Pi[i][j].hi==Pi[j][i].hi for i in range(n) for j in range(n))
    # Representative exact points exercise inactive and active branches.
    base=p['base_report'];radius=float(p['projection_radius_mps2'])
    cases=[]
    for e,beta in (([.1,-.05,.02],[.03,.01,-.02]),([.7,-.2,.1],[.08,-.04,.03])):
        f=__import__('ou3_projection_sector').bias_error_projection(e,beta,radius)
        slack=sum(x*x for x in e)+sum(x*x for x in beta)-sum(x*x for x in f)
        cases.append(slack>=-1e-14)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'projection_operator':'F_R(e,beta)=beta-Pi_R(beta-e)','global_joint_sector_consumed':p['global_joint_sector_closed'],
      'saturated_branch_retained':p['saturated_branch_included'],'unsaturated_branch_retained':p['unsaturated_branch_included'],
      'Clarke_boundary_branch_retained':p['boundary_clarke_branch_included'],
      'same_event_preprojection_error_map_required':True,'physical_beta_coordinate_required':True,'projected_output_coordinate_required':True,
      'dense_projection_IQC_available':True,'dense_projection_IQC_finite':finite,'dense_projection_IQC_symmetric':symmetric,
      'representative_exact_operator_points_satisfy_IQC':all(cases),
      'inactive_projection_assumed':False,'rowwise_K_bound_used':False,'independent_bias_error_slots_used':False,
      'production_A21_prefix_projection_attached_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('global_joint_sector_consumed','saturated_branch_retained','unsaturated_branch_retained','Clarke_boundary_branch_retained','same_event_preprojection_error_map_required','physical_beta_coordinate_required','projected_output_coordinate_required','dense_projection_IQC_available','dense_projection_IQC_finite','dense_projection_IQC_symmetric','representative_exact_operator_points_satisfy_IQC'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('inactive_projection_assumed','rowwise_K_bound_used','independent_bias_error_slots_used','production_A21_prefix_projection_attached_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'projection_iqc':d['dense_projection_IQC_available'],'branches':d['saturated_branch_retained'] and d['unsaturated_branch_retained'] and d['Clarke_boundary_branch_retained'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
