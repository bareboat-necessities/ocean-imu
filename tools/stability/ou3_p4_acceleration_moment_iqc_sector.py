#!/usr/bin/env python3
"""Lift the same-history BRMM acceleration-moment IQC into a P4 graph sector.

The physical sea forcing and the initial P4 error are independent theorem
quantifiers.  Therefore the acceleration-moment bound MUST NOT use the hard-entry
radial coordinate.  In particular, zero initial error does not imply zero sea.

For one physical transition, normalized moment coordinates
x=(J0/h,J1/h^2,J2/h^3) satisfy

    sum_axis x_axis^T G^-1 x_axis <= A_max^2 h_s^2,

where h_s is the homogeneous physical-source scale.  Production uses h_s=1.
It is a source-coordinate homogenizer, distinct from the hard-entry radial r_e
used to cover the error set.

If h_s=C_s z and the nine interleaved-by-axis moment coordinates m=C_m z in the
COMMON endpoint/prefix augmented coordinate, the admissible source graph obeys

    z^T Pi_mom z >= 0,
    Pi_mom = A_max^2 C_s^T C_s - C_m^T (I_3 kron G^-1) C_m.

The same physical-source homogeneous coordinate may be shared by all transition
moment sectors in a prefix; it must never be identified with SourceCoverCell's
hard-entry radial_scale.  This corrects the earlier radial/source conflation.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from typing import Sequence

from ou3_interval import Interval, matrix_mul, matrix_sub, matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_brmm_acceleration_moment_iqc as MOM
import ou3_p4_complete_brmm_joint_sector_master as JOINT

SCHEMA=2
QUALIFICATION='OU3_P4_ACCELERATION_MOMENT_JOINT_SECTOR_V2'


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _exact_number(x):
    """Accept decimal or exact rational JSON scalars without changing semantics."""
    return float(Fraction(str(x)))


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
    physical_source_scale_map: Sequence[Sequence[Interval]],
    A_max: float,
):
    """Return Pi with z^T Pi z>=0; source scale is NOT hard-entry radial."""
    mr,mz=_shape(normalized_moment_map)
    sr,sz=_shape(physical_source_scale_map)
    if mr!=9 or sr!=1 or mz==0 or sz!=mz:
        raise ValueError('need 9xn moment map and 1xn physical-source scale map')
    if not (float(A_max)>0):
        raise ValueError('positive acceleration cap required')
    positive=matrix_mul(matrix_mul(matrix_transpose(physical_source_scale_map),
                                   [[Interval.point(float(A_max)**2)]]),
                        physical_source_scale_map)
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
    A=_exact_number(mom['A_max_mps2'])

    # Common coordinate z=[h_source, x0x,x1x,x2x, x0y,...,x2z].
    # h_source=1 is a homogeneous physical-source coordinate, NOT entry radial.
    n=10
    Cs=_selector(1,n,[0])
    Cm=_selector(9,n,list(range(1,10)))
    Pi=moment_sector(Cm,Cs,A)

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
      'physical_source_homogeneous_scale_used':True,
      'physical_source_scale_value_for_production':1.0,
      'physical_source_scale_independent_of_hard_entry_radial':True,
      'hard_entry_radial_used_for_physical_source_bound':False,
      'same_radial_coordinate_must_parameterize_entire_source_history':False,
      'zero_initial_error_may_have_nonzero_physical_source':True,
      'normalized_moment_coordinate_dimension':9,
      'joint_three_axis_Ginv_block_materialized':True,
      'nonnegative_S_procedure_multiplier_compatible':True,
      'independent_source_supply_port_used':False,
      'independent_moment_or_axis_box_used':False,
      'moment_sector_matrix_available':True,
      'aligned_constant_acceleration_boundary_contains_zero':aligned_contains_zero,
      'detached_box_corner_strictly_excluded':detached_strict_negative,
      'zero_homogeneous_source_zero_moment_boundary_contains_zero':zero_contains_zero,
      'P3_delta':1e-18,'P4_PASS':False,'P5_MAY_START':False,
      'production_source_map_bound_here':False,
      'next_obligation':(
        'embed the 9 normalized moments and a homogeneous physical-source scale h_s=1 in each '
        'reachable endpoint/prefix map while retaining the independent hard-entry radial r_e; '
        'include every transition Pi_mom with reset/projection/residual sectors in outward LDLT')}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema mismatch')
    for k in ('common_augmented_coordinate_required','physical_source_homogeneous_scale_used','physical_source_scale_independent_of_hard_entry_radial','zero_initial_error_may_have_nonzero_physical_source','joint_three_axis_Ginv_block_materialized','nonnegative_S_procedure_multiplier_compatible','moment_sector_matrix_available','aligned_constant_acceleration_boundary_contains_zero','detached_box_corner_strictly_excluded','zero_homogeneous_source_zero_moment_boundary_contains_zero'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('hard_entry_radial_used_for_physical_source_bound','same_radial_coordinate_must_parameterize_entire_source_history','independent_source_supply_port_used','independent_moment_or_axis_box_used','P4_PASS','P5_MAY_START','production_source_map_bound_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('physical_source_scale_value_for_production')!=1.0:f.append('physical source homogeneous scale changed')
    if d.get('normalized_moment_coordinate_dimension')!=9:f.append('moment dimension changed')
    if d.get('P3_delta')!=1e-18:f.append('P3 changed')
    return f


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    from pathlib import Path
    p=Path(a.output);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'sector':d['moment_sector_matrix_available'],'source_scale_separate':d['physical_source_scale_independent_of_hard_entry_radial'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
