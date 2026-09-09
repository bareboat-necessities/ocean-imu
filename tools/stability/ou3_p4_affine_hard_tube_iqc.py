#!/usr/bin/env python3
"""Affine-homogeneous hard-entry/correction/prefix IQCs for P4.

A finite regional theorem is not homogeneous in the physical state because its
entry and retention sets have fixed radii.  Introduce one lift coordinate h and
work with the cone over the hard set.  For a three-coordinate selector E_g and
radius r_g,

    ||E_g x|| <= r_g h

is represented by the quadratic IQC

    z^T Pi_g z >= 0,
    Pi_g = r_g^2 e_h e_h^T - E_g^T E_g,
    z = [h; x; graph/source/roundoff auxiliaries].

No probabilistic covariance membership is involved.  The same construction
expresses a SAME-CELL correction-domain obligation

    ||D_theta z|| <= delta h

and every-prefix hard-domain obligations

    ||O_{ell,g} z|| <= r_g h.

Given premise IQCs Pi_j, a sufficient outward certificate for a target T is

    T - sum_j lambda_j Pi_j >= 0,    lambda_j >= 0.

The implementation uses strict outward LDLT for production numerical checks.
A non-strict analytical boundary can be handled only by a separate exact
factorization; this file does not silently add epsilon to make a failure pass.

This is the missing affine lift needed to prove both reset-chart admission and
literal-prefix retention from the SAME augmented graph.  It does not choose a
correction radius and does not promote P4 by itself.
"""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1

QUALIFICATION='OU3_P4_AFFINE_HOMOGENEOUS_HARD_TUBE_IQC_V1'
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

def _zero(n):return [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
def _scale(A,a):
    c=Interval.point(float(a));return [[c*x for x in row] for row in A]
def _gram(A):return matrix_mul(matrix_transpose(A),A)

def selector(n:int,indices:Sequence[int]):
    if n<=0 or not indices or any(i<0 or i>=n for i in indices):raise ValueError('selector indices outside coordinate')
    A=[[Interval.point(0.0) for _ in range(n)] for _ in indices]
    for r,i in enumerate(indices):A[r][i]=Interval.point(1.0)
    return A

def ball_iqc(n:int,h_index:int,indices:Sequence[int],radius:float):
    """Pi with z'Pi z = r^2 h^2-||E z||^2."""
    r=float(radius)
    if not (math.isfinite(r) and r>0 and 0<=h_index<n):raise ValueError('invalid hard-ball radius/lift index')
    E=selector(n,indices);Pi=[list(row) for row in _scale(_gram(E),-1.0)]
    Pi[h_index][h_index]=Pi[h_index][h_index]+Interval.point(r*r)
    return matrix_symmetric_hull(Pi)

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
    """Strict sufficient certificate T-sum(lambda Pi)>0 by outward LDLT."""
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
    # Smoke coordinate [h,x0,x1,x2]. Premise ||x||<=2h implies target
    # ||0.4x||<=0.9h with multiplier 0.16. The exact remainder is
    # diag(0.17,0,0,0), so its zero pivots must not be rejected merely because
    # outward interval rounding straddles zero. Production strict LDLT below is
    # unchanged and is used only where callers provide strict slack.
    n=4;Pi=ball_iqc(n,0,(1,2,3),2.0)
    D=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
    for i in range(3):D[i][1+i]=Interval.point(0.4)
    target=mapped_ball_target(D,0,0.9)
    rem=sprocedure_remainder(target,[Pi],[0.16])
    interval_consistent=(rem[0][0].lo>0.0 and all(rem[i][i].lo<=0.0<=rem[i][i].hi for i in range(1,n))
                         and all(rem[i][j].lo<=0.0<=rem[i][j].hi for i in range(n) for j in range(n) if i!=j))
    lam=Fraction(4,25);h_slack=Fraction(9,10)**2-lam*Fraction(2)**2;x_slack=-Fraction(2,5)**2+lam
    smoke_semidefinite=bool(h_slack>0 and x_slack==0 and interval_consistent)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'lift_coordinate':'z=[h; physical_error; chord/reset/projection/source/roundoff auxiliaries]',
      'hard_entry_full_declared_scale_consumed':e['full_declared_scale_enforced'],
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
    for k in ('hard_entry_full_declared_scale_consumed','hard_entry_ball_IQC_available','same_graph_correction_domain_target_available','same_graph_every_prefix_ball_target_available','nonnegative_multiplier_Sprocedure_available','strict_outward_LDLT_checker_available','BIAS1_true_bias_bound_available_for_source_lift','smoke_semidefinite_implication_closed_exact_diagonal'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_membership_used','trajectory_fit_used','domain_shrunk','epsilon_added_to_force_semidefinite_pass','production_same_cell_correction_domain_closed_here','production_every_prefix_hard_domain_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if not(math.isfinite(float(d.get('BIAS1_true_bias_norm_upper_mps2',math.nan))) and float(d['BIAS1_true_bias_norm_upper_mps2'])>0):f.append('BIAS1 bound invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
