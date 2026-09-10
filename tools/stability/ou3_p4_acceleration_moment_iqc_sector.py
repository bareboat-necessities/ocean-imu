#!/usr/bin/env python3
"""Lift the same-history BRMM acceleration-moment IQC into a P4 graph sector.

For a radialized source cell with common scale r in [0,1], normalized moment
coordinates x=(J0/h,J1/h^2,J2/h^3) satisfy

    sum_axis x_axis^T G^-1 x_axis <= A_max^2 r^2.

If r=C_r z and the nine interleaved-by-axis moment coordinates m=C_m z in the
COMMON endpoint/prefix augmented coordinate, the admissible source graph obeys

    z^T Pi_mom z >= 0,
    Pi_mom = A_max^2 C_r^T C_r - C_m^T (I_3 kron G^-1) C_m.

This is directly consumable by the existing nonnegative-multiplier joint
S-procedure.  It is not a detached supply port, and the same radial coordinate
must also parameterize the source/entry history used by the rest of the word.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_brmm_acceleration_moment_iqc as MOM
import ou3_p4_complete_brmm_joint_sector_master as JOINT

SCHEMA=1
QUALIFICATION='OU3_P4_ACCELERATION_MOMENT_JOINT_SECTOR_V1'


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _block_ginv():
    z=Interval.point(0.0)
    W=[[z for _ in range(9)] for _ in range(9)]
    for axis in range(3):
        for i in range(3):
            for j in range(3):
                W[3*axis+i][3*axis+j]=Interval.point(float(MOM.GINV[i][j]))
    return W


def moment_sector(
    normalized_moment_map: Sequence[Sequence[Interval]],
    radial_scale_map: Sequence[Sequence[Interval]],
    A_max: float,
):
    """Return Pi with convention z^T Pi z>=0 on admitted radial source cells."""
    mr,mz=_shape(normalized_moment_map)
    rr,rz=_shape(radial_scale_map)
    if mr!=9 or rr!=1 or mz==0 or rz!=mz:
        raise ValueError('need 9xn moment map and 1xn common radial map')
    if not (float(A_max)>0):
        raise ValueError('positive acceleration cap required')
    positive=matrix_mul(matrix_mul(matrix_transpose(radial_scale_map),
                                   [[Interval.point(float(A_max)**2)]]),
                        radial_scale_map)
    negative=matrix_mul(matrix_mul(matrix_transpose(normalized_moment_map),
                                   _block_ginv()),
                        normalized_moment_map)
    return matrix_symmetric_hull(matrix_sub(positive,negative))


def _selector(rows,n,indices):
    if len(indices)!=rows:
        raise ValueError('selector index count mismatch')
    z=Interval.point(0.0); o=Interval.point(1.0)
    A=[[z for _ in range(n)] for _ in range(rows)]
    for r,c in enumerate(indices): A[r][c]=o
    return A


def build():
    mom=MOM.build(); mf=MOM.validate(mom)
    joint=JOINT.build(); jf=JOINT.validate(joint)
    if mf or jf:
        raise RuntimeError(f'prerequisite failure moment={mf} joint={jf}')
    A=float(mom['A_max_mps2'])

    # Common coordinate z=[r, x0x,x1x,x2x, x0y,...,x2z].
    n=10
    Cr=_selector(1,n,[0])
    Cm=_selector(9,n,list(range(1,10)))
    Pi=moment_sector(Cm,Cr,A)

    def q(v): return JOINT.quadratic_value(Pi,[Interval.point(float(x)) for x in v])
    aligned=[1.0,A,A/2,A/6,0,0,0,0,0,0]
    detached=[1.0,A,-A/2,A/6,0,0,0,0,0,0]
    zero=[0.0]*n
    qa=q(aligned); qd=q(detached); qz=q(zero)
    aligned_contains_zero=qa.lo<=0<=qa.hi
    detached_strict_negative=qd.hi<0
    zero_contains_zero=qz.lo<=0<=qz.hi
    if not (aligned_contains_zero and detached_strict_negative and zero_contains_zero):
        raise RuntimeError('moment sector smoke failed')

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'sector_convention':'z^T Pi_mom z >= 0',
      'common_augmented_coordinate_required':True,
      'same_radial_coordinate_must_parameterize_entire_source_history':True,
      'radial_domain':'closed [0,1]',
      'normalized_moment_coordinate_dimension':9,
      'joint_three_axis_Ginv_block_materialized':True,
      'nonnegative_S_procedure_multiplier_compatible':True,
      'independent_source_supply_port_used':False,
      'independent_moment_or_axis_box_used':False,
      'constant_affine_one_coordinate_used':False,
      'moment_sector_matrix_available':True,
      'aligned_constant_acceleration_boundary_contains_zero':aligned_contains_zero,
      'detached_box_corner_strictly_excluded':detached_strict_negative,
      'zero_radial_zero_moment_boundary_contains_zero':zero_contains_zero,
      'P3_delta':1e-18,'P4_PASS':False,'P5_MAY_START':False,
      'production_source_map_bound_here':False,
      'next_obligation':(
        'embed the 9 normalized moments and the SAME source radial coordinate in each '
        'reachable provider/joint24 endpoint and literal-prefix map; then include this Pi '
        'with the reset/projection/residual sectors in outward augmented LDLT')}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema mismatch')
    for k in ('common_augmented_coordinate_required','same_radial_coordinate_must_parameterize_entire_source_history','joint_three_axis_Ginv_block_materialized','nonnegative_S_procedure_multiplier_compatible','moment_sector_matrix_available','aligned_constant_acceleration_boundary_contains_zero','detached_box_corner_strictly_excluded','zero_radial_zero_moment_boundary_contains_zero'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_source_supply_port_used','independent_moment_or_axis_box_used','constant_affine_one_coordinate_used','P4_PASS','P5_MAY_START','production_source_map_bound_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('normalized_moment_coordinate_dimension')!=9:f.append('moment dimension changed')
    if d.get('P3_delta')!=1e-18:f.append('P3 changed')
    return f


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'sector':d['moment_sector_matrix_available'],'detached_excluded':d['detached_box_corner_strictly_excluded'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
