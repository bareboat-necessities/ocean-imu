#!/usr/bin/env python3
"""Affine-homogeneous hard-entry/correction/prefix IQCs for P4.

A finite regional theorem is not homogeneous in the physical state because its
entry and retention sets have fixed radii. Introduce one lift coordinate h and
work with the cone over the hard set. For a selector E_g and radius r_g >= 0,

    ||E_g x|| <= r_g h

is represented by

    z^T Pi_g z >= 0,
    Pi_g = r_g^2 e_h e_h^T - E_g^T E_g.

The r_g=0 case is intentional and exact: -||E_g x||^2 >= 0 is equivalent to
E_g x=0. This is how the shipping-reachable handoff fiber represents position
and integral-displacement error at its local origin. It must not be replaced by
an epsilon-radius ball.

Selector-ball IQCs are constructed directly on their known diagonal structure.
Using generic outward interval matrix multiplication for E^T E needlessly
widens structural zeros to tiny signed intervals; at r=0 that would turn an
exact equality fiber into a numerical fuzzy cone. Direct construction retains
exact zero and exact -1 coefficients, while a nonzero r^2 coefficient alone is
outward rounded.
"""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction
from typing import Sequence

from ou3_interval import Interval,down,up,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1

QUALIFICATION='OU3_P4_AFFINE_HOMOGENEOUS_HARD_TUBE_IQC_V3'
GROUPS={
 'attitude_cayley_norm':(0,1,2),
 'gyro_bias_norm_rad_s':(3,4,5),
 'velocity_norm_mps':(6,7,8),
 'position_norm_m':(9,10,11),
 'integral_displacement_norm_m_s':(12,13,14),
 'latent_acceleration_norm_mps2':(15,16,17),
 'accelerometer_bias_error_norm_mps2':(18,19,20),
}

def _shape(A):return len(A),len(A[0]) if A else 0

def _scale(A,a):
    c=Interval.point(float(a));return [[c*x for x in row] for row in A]
def _gram(A):return matrix_mul(matrix_transpose(A),A)

def selector(n:int,indices:Sequence[int]):
    if n<=0 or not indices or any(i<0 or i>=n for i in indices):raise ValueError('selector indices outside coordinate')
    if len(set(indices)) != len(tuple(indices)):raise ValueError('selector indices must be unique')
    A=[[Interval.point(0.0) for _ in range(n)] for _ in indices]
    for r,i in enumerate(indices):A[r][i]=Interval.point(1.0)
    return A

def ball_iqc(n:int,h_index:int,indices:Sequence[int],radius:float):
    """Exact selector IQC r^2 h^2-sum(x_i^2), including exact r=0 fiber.

    ``E`` is a coordinate selector, so E^T E is known exactly and no interval
    matrix product is needed. Structural zeros and -1 coefficients are point
    intervals. Only a positive binary64 radius-square is outward rounded.
    """
    r=float(radius);idx=tuple(indices)
    if not (math.isfinite(r) and r>=0 and 0<=h_index<n):raise ValueError('invalid hard-ball/fiber radius or lift index')
    if not idx or any(i<0 or i>=n for i in idx) or len(set(idx))!=len(idx):raise ValueError('invalid/duplicate hard-ball selector')
    if h_index in idx:raise ValueError('lift coordinate may not be one of the selected physical coordinates')
    Pi=[[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
    for i in idx:Pi[i][i]=Interval.point(-1.0)
    if r==0.0:
        Pi[h_index][h_index]=Interval.point(0.0)
    else:
        rr=r*r
        Pi[h_index][h_index]=Interval(down(rr),up(rr))
    return Pi

def mapped_ball_target(Dmap:Sequence[Sequence[Interval]],h_index:int,radius:float):
    """T with z'Tz = r^2 h^2-||D z||^2."""
    rows,n=_shape(Dmap);r=float(radius)
    if rows==0 or n==0 or not (0<=h_index<n) or not(math.isfinite(r) and r>0):raise ValueError('invalid target map/radius')
    T=_scale(_gram(Dmap),-1.0);T[h_index][h_index]=T[h_index][h_index]+Interval.point(r*r)
    return matrix_symmetric_hull(T)

def sprocedure_remainder(target,premises,multipliers):
    n,m=_shape(target)
    if n==0 or n!=m or len(premises)!=len(multipliers):raise ValueError('target/premise dimensions')
    out=matrix_symmetric_hull(target)
    for Pi,lam in zip(premises,multipliers):
        if _shape(Pi)!=(n,n):raise ValueError('premise dimension mismatch')
        x=float(lam)
        if not(math.isfinite(x) and x>=0):raise ValueError('nonnegative finite multiplier required')
        out=matrix_symmetric_hull(matrix_sub(out,_scale(Pi,x)))
    return out

def certify_strict_target(target,premises,multipliers):
    R=sprocedure_remainder(target,premises,multipliers);ok,p=symmetric_positive_definite_ldlt(R)
    return bool(ok),[float(x.lo) for x in p]

def hard_entry_iqcs(mode:str,n:int,h_index:int=0,state_offset:int=1):
    e=ENTRY.build();bad=ENTRY.validate(e)
    if bad:raise RuntimeError('hard entry prerequisite failed: '+repr(bad))
    if mode not in ('H18','A21'):raise ValueError('mode must be H18/A21')
    dim=18 if mode=='H18' else 21
    if n<state_offset+dim:raise ValueError('augmented coordinate too small for state')
    names=list(e['H18_groups'])+([e['A21_additional_group']] if mode=='A21' else [])
    out={}
    for name in names:
        idx=tuple(state_offset+i for i in GROUPS[name]);out[name]=ball_iqc(n,h_index,idx,float(e['coordinate_radii'][name]))
    return out

def build():
    e=ENTRY.build();b=BIAS1.build();bad={'entry':ENTRY.validate(e),'bias1':BIAS1.validate(b)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('affine hard-tube prerequisites failed: '+repr(bad))
    n=4;Pi=ball_iqc(n,0,(1,2,3),2.0)
    D=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
    for i in range(3):D[i][1+i]=Interval.point(0.4)
    target=mapped_ball_target(D,0,0.9)
    rem=sprocedure_remainder(target,[Pi],[0.16])
    interval_consistent=(rem[0][0].lo>0.0 and all(rem[i][i].lo<=0.0<=rem[i][i].hi for i in range(1,n))
                         and all(rem[i][j].lo<=0.0<=rem[i][j].hi for i in range(n) for j in range(n) if i!=j))
    lam=Fraction(4,25);h_slack=Fraction(9,10)**2-lam*Fraction(2)**2;x_slack=-Fraction(2,5)**2+lam
    smoke_semidefinite=bool(h_slack>0 and x_slack==0 and interval_consistent)

    zero_pi=ball_iqc(2,0,(1,),0.0)
    zero_fiber_exact=(zero_pi[0][0]==Interval.point(0.0)
                      and zero_pi[1][1]==Interval.point(-1.0)
                      and zero_pi[0][1]==Interval.point(0.0)
                      and zero_pi[1][0]==Interval.point(0.0))
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'lift_coordinate':'z=[h; physical_error; chord/reset/projection/source/roundoff auxiliaries]',
      'hard_entry_full_declared_scale_consumed':e['full_declared_scale_enforced'],
      'hard_entry_correlated_fiber_consumed':e['entry_is_correlated_fiber_not_cartesian_box'],
      'selector_gram_built_directly_not_by_interval_product':True,
      'exact_zero_radius_fiber_IQC_supported':True,
      'exact_zero_radius_fiber_smoke_closed':zero_fiber_exact,
      'covariance_membership_used':False,'trajectory_fit_used':False,'domain_shrunk':False,
      'hard_entry_ball_IQC_available':True,'same_graph_correction_domain_target_available':True,
      'same_graph_every_prefix_ball_target_available':True,'nonnegative_multiplier_Sprocedure_available':True,
      'strict_outward_LDLT_checker_available':True,'epsilon_added_to_force_semidefinite_pass':False,
      'BIAS1_true_bias_bound_available_for_source_lift':True,
      'BIAS1_true_bias_norm_upper_mps2':b['true_bias_norm_upper_mps2'],
      'smoke_semidefinite_implication_closed_exact_diagonal':smoke_semidefinite,
      'production_same_cell_correction_domain_closed_here':False,
      'production_every_prefix_hard_domain_closed_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('hard_entry_full_declared_scale_consumed','hard_entry_correlated_fiber_consumed',
              'selector_gram_built_directly_not_by_interval_product','exact_zero_radius_fiber_IQC_supported','exact_zero_radius_fiber_smoke_closed',
              'hard_entry_ball_IQC_available','same_graph_correction_domain_target_available',
              'same_graph_every_prefix_ball_target_available','nonnegative_multiplier_Sprocedure_available',
              'strict_outward_LDLT_checker_available','BIAS1_true_bias_bound_available_for_source_lift',
              'smoke_semidefinite_implication_closed_exact_diagonal'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_membership_used','trajectory_fit_used','domain_shrunk','epsilon_added_to_force_semidefinite_pass','production_same_cell_correction_domain_closed_here','production_every_prefix_hard_domain_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if not(math.isfinite(float(d.get('BIAS1_true_bias_norm_upper_mps2',math.nan))) and float(d['BIAS1_true_bias_norm_upper_mps2'])>0):f.append('BIAS1 bound invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
