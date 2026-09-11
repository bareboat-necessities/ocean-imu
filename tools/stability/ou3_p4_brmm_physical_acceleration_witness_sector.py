#!/usr/bin/env python3
"""Joint quadratic outer relation for the 15D BRMM acceleration witness.

One prediction forcing block is generated from

  q = [a0_xyz, a1_xyz, J0_xyz, J1_xyz, J2_xyz]

by ``ou3_p4_brmm_physical_prediction_forcing``.  The source relation never
replaces q by independent defect boxes.  For an admitted physical acceleration
history with ||a(t)||<=A, the same q necessarily satisfies

  ||a0||^2 <= A^2,
  ||a1||^2 <= A^2,
  sum_axis x_axis^T G^-1 x_axis <= A^2,

where x=(J0/h,J1/h^2,J2/h^3).  These three inequalities form a rigorous
quadratic outer relaxation of the SAME 15D witness.  They do not assert a
converse and do not claim that endpoint values are independent of the moments;
the exact physical-history ancestry remains retained by the source cover.

All inequalities are homogenized with a physical-source coordinate h_s=1,
which is explicitly distinct from the P4 initial-error radial coordinate r_e.
Thus zero initial error may coexist with a full-amplitude admissible sea.
"""
from __future__ import annotations
import argparse,json
from typing import Sequence
from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_brmm_acceleration_moment_iqc as MOM
import ou3_p4_acceleration_moment_iqc_sector as MOMSECTOR
import ou3_p4_brmm_physical_prediction_forcing as FORCING

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_PHYSICAL_ACCELERATION_WITNESS_SECTOR_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))
def shape(A):return len(A),len(A[0]) if A else 0
def gram(A):return matrix_mul(matrix_transpose(A),A)
def scale_gram(smap,A):
    return matrix_mul(matrix_mul(matrix_transpose(smap),[[I(A*A)]]),smap)
def ball_sector(vmap,smap,A):
    vr,n=shape(vmap);sr,m=shape(smap)
    if vr!=3 or sr!=1 or n==0 or n!=m:raise ValueError('need 3xn vector and 1xn source-scale maps')
    return matrix_symmetric_hull(matrix_sub(scale_gram(smap,A),gram(vmap)))
def physical_witness_sectors(a0_map,a1_map,normalized_moment_map,source_scale_map,A):
    if shape(normalized_moment_map)[0]!=9:raise ValueError('normalized moment map must have 9 rows')
    n=shape(normalized_moment_map)[1]
    if shape(a0_map)!=(3,n) or shape(a1_map)!=(3,n) or shape(source_scale_map)!=(1,n):raise ValueError('all witness maps must share one coordinate')
    return (
      ('physical_a0_ball',ball_sector(a0_map,source_scale_map,A)),
      ('physical_a1_ball',ball_sector(a1_map,source_scale_map,A)),
      ('physical_J012_moment_iqc',MOMSECTOR.moment_sector(normalized_moment_map,source_scale_map,A)),
    )
def selector(rows,n,indices):
    if len(indices)!=rows:raise ValueError('selector mismatch')
    z=I(0);o=I(1);A=[[z for _ in range(n)] for _ in range(rows)]
    for r,c in enumerate(indices):A[r][c]=o
    return A

def canonical_maps(h:float):
    """Coordinate z=[h_s,q15] with q order from physical forcing module."""
    if h<=0:raise ValueError('positive sample required')
    n=16;Cs=selector(1,n,[0]);A0=selector(3,n,[1,2,3]);A1=selector(3,n,[4,5,6])
    z=I(0);M=[[z for _ in range(n)] for _ in range(9)]
    # MOM uses axis-major [J0/h,J1/h^2,J2/h^3] for x,y,z.
    for axis in range(3):
        M[3*axis+0][7+axis]=I(1/h)
        M[3*axis+1][10+axis]=I(1/(h*h))
        M[3*axis+2][13+axis]=I(1/(h*h*h))
    return A0,A1,M,Cs
def build():
    mom=MOM.build();mf=MOM.validate(mom);forcing=FORCING.build();ff=FORCING.validate(forcing);ms=MOMSECTOR.build();sf=MOMSECTOR.validate(ms)
    if mf or ff or sf:raise RuntimeError(f'witness sector prerequisites failed moment={mf} forcing={ff} sector={sf}')
    h=.005;A=float(mom['A_max_mps2']);maps=canonical_maps(h);sectors=physical_witness_sectors(*maps,A)
    if len(sectors)!=3 or any(shape(P)!=(16,16) for _,P in sectors):raise RuntimeError('witness sector shape failure')
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'physical_witness_coordinate_order':forcing['source_coordinate_order'],'physical_witness_dimension':15,'homogeneous_physical_source_scale_dimension':1,
      'same_15D_witness_drives_all_four_prediction_forcing_rows':True,'a0_and_a1_endpoint_balls_attached_to_same_witness':True,'J0_J1_J2_coupled_three_axis_IQC_attached_to_same_witness':True,
      'physical_source_scale_is_homogeneous_one':True,'physical_source_scale_independent_of_hard_entry_radial':True,'hard_entry_radial_used_in_source_sector':False,
      'zero_initial_error_nonzero_source_allowed':True,'independent_prediction_defect_boxes_used':False,'independent_J_moment_boxes_used':False,'converse_outer_relation_claimed':False,
      'quadratic_sector_names':[x[0] for x in sectors],'joint_15D_source_quadratic_outer_relation_closed':True,
      'prefix_response_columns_attached_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'append one q15 physical witness block at each prediction, apply the exact forcing injection, suffix-propagate it through all later literal event maps, and embed these three sectors per transition in the same augmented prefix coordinate'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('same_15D_witness_drives_all_four_prediction_forcing_rows','a0_and_a1_endpoint_balls_attached_to_same_witness','J0_J1_J2_coupled_three_axis_IQC_attached_to_same_witness','physical_source_scale_is_homogeneous_one','physical_source_scale_independent_of_hard_entry_radial','zero_initial_error_nonzero_source_allowed','joint_15D_source_quadratic_outer_relation_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('hard_entry_radial_used_in_source_sector','independent_prediction_defect_boxes_used','independent_J_moment_boxes_used','converse_outer_relation_claimed','prefix_response_columns_attached_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('physical_witness_dimension')!=15:f.append('physical witness dimension changed')
    if d.get('quadratic_sector_names')!=['physical_a0_ball','physical_a1_ball','physical_J012_moment_iqc']:f.append('physical sector set changed')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    from pathlib import Path
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'source_sector':d['joint_15D_source_quadratic_outer_relation_closed'],'separate_scale':d['physical_source_scale_independent_of_hard_entry_radial'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
