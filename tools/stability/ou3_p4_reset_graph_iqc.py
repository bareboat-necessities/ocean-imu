#!/usr/bin/env python3
"""Parameterized dense finite-reset graph IQC tied to SAME Joseph correction.

For one accepted shipping update let d=K y, d_theta=E_theta d and
b=G(d_theta)^-1 rho.  Given a SAME-GRAPH certified correction domain
||d_theta||<=delta h, the reset primitive supplies a uniform gain

    ||b|| <= ||rho|| <= mu_R(delta) ||d_theta||.

If Y maps the common augmented coordinate z to the physical residual y and B
maps z to b, then D=E_theta*K*Y uses the actual same P/H/R/K cell and

    z^T [mu_R^2 D^T D - B^T B] z >= 0.

The correction-domain proof is separate but must use this same D map and the
affine hard-tube IQCs.  No rowwise/global K bound is accepted as delta.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_reset_signed_absorption_diagnostic as RESET_SECTOR
import ou3_p4_affine_hard_tube_iqc as HARD

QUALIFICATION='OU3_P4_PARAMETERIZED_SAME_CELL_FINITE_RESET_GRAPH_IQC_V3'

def _shape(A):return len(A),len(A[0]) if A else 0
def _scale(A,a):
    x=Interval.point(float(a));return [[x*v for v in row] for row in A]
def _gram(A):return matrix_mul(matrix_transpose(A),A)

def reset_sector_matrix(Dmap:Sequence[Sequence[Interval]],Bmap:Sequence[Sequence[Interval]],gamma_reset:float):
    dr,n=_shape(Dmap);br,m=_shape(Bmap)
    if dr!=3 or br!=3 or n==0 or m!=n:raise ValueError('Dmap/Bmap must be 3xn')
    g=float(gamma_reset)
    if not(math.isfinite(g) and g>=0):raise ValueError('finite nonnegative reset gain required')
    return matrix_symmetric_hull(matrix_sub(_scale(_gram(Dmap),g*g),_gram(Bmap)))

def same_cell_correction_map(K:Sequence[Sequence[Interval]],Ymap:Sequence[Sequence[Interval]]):
    kr,kc=_shape(K);yr,yc=_shape(Ymap)
    if kr not in (18,21) or kc!=3 or yr!=3 or yc==0:raise ValueError('K must be H18/A21 x3 and Ymap 3xn')
    return matrix_mul([list(row) for row in K[:3]],Ymap)

def correction_domain_target(Dmap,h_index:int,delta:float):
    """Affine target delta^2 h^2-||D z||^2>=0 for same-cell reset validity."""
    return HARD.mapped_ball_target(Dmap,h_index,float(delta))

def parameterized_reset_iqc(Dmap,Bmap,q:float,delta:float):
    sector=RESET_SECTOR.homogeneous_sector(float(q),float(delta))
    Pi=reset_sector_matrix(Dmap,Bmap,float(sector['reset_defect_over_correction_norm_upper']))
    return Pi,sector

def kalman_correction_information_identity(PJinv,K,Rinv,Sinv):
    n,n2=_shape(PJinv);kr,kc=_shape(K)
    if n not in (18,21) or n2!=n or (kr,kc)!=(n,3) or _shape(Rinv)!=(3,3) or _shape(Sinv)!=(3,3):raise ValueError('identity dimensions')
    lhs=matrix_mul(matrix_mul(matrix_transpose(K),PJinv),K);rhs=matrix_sub(Rinv,Sinv)
    return matrix_symmetric_hull(matrix_sub(lhs,rhs))
def _zero_contained(A):return all(v.lo<=0<=v.hi for row in A for v in row)

def build():
    e=ENTRY.build();rs=RESET_SECTOR.build();hard=HARD.build();bad={'entry':ENTRY.validate(e),'reset':RESET_SECTOR.validate(rs),'hard':HARD.validate(hard)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset IQC prerequisites failed: '+repr(bad))
    q=float(e['coordinate_radii']['attitude_cayley_norm'])
    # Algebra-only smoke at delta=.1. It is not a source correction claim.
    delta=.1;n=9;Y=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
    for i in range(3):Y[i][1+i]=Interval.point(1.0)
    K=[[Interval.point(0.0) for _ in range(3)] for _ in range(18)]
    for i in range(3):K[i][i]=Interval.point(.25+.05*i)
    D=same_cell_correction_map(K,Y);sector=RESET_SECTOR.homogeneous_sector(q,delta);gamma=float(sector['reset_defect_over_correction_norm_upper'])
    B=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
    for i in range(3):B[i][1+i]=Interval.point(gamma*(.25+.05*i))
    Pi=reset_sector_matrix(D,B,gamma)
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'exact_K_PJinv_K_identity':'K^T*P_J^-1*K=R^-1-S^-1','same_cell_K_required':True,
      'independent_correction_port_forbidden':True,'parameterized_reset_defect_dense_IQC_available':True,
      'reset_IQC_keeps_residual_direction_via_Etheta_K_Y':True,'reset_inverse_operator_norm_upper':1.0,
      'production_delta_must_be_certified_from_same_Dmap':True,'affine_correction_domain_target_available':True,
      'rowwise_K_correction_domain_forbidden':True,'parameterized_homogeneous_reset_sector_consumed':True,
      'smoke_delta_not_theorem_bound':delta,'smoke_reset_gain':gamma,'smoke_reset_sector_zero_contained':_zero_contained(Pi),
      'production_same_graph_correction_domain_closed_here':False,'source_uniform_augmented_LDLT_closed_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_K_required','independent_correction_port_forbidden','parameterized_reset_defect_dense_IQC_available',
              'reset_IQC_keeps_residual_direction_via_Etheta_K_Y','production_delta_must_be_certified_from_same_Dmap',
              'affine_correction_domain_target_available','rowwise_K_correction_domain_forbidden',
              'parameterized_homogeneous_reset_sector_consumed','smoke_reset_sector_zero_contained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_same_graph_correction_domain_closed_here','source_uniform_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('reset_inverse_operator_norm_upper')!=1.0:f.append('reset inverse norm changed')
    if not(0<float(d.get('smoke_delta_not_theorem_bound',0))<3):f.append('smoke delta invalid')
    if not(math.isfinite(float(d.get('smoke_reset_gain',math.nan))) and float(d['smoke_reset_gain'])>=0):f.append('smoke reset gain invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
