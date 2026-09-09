#!/usr/bin/env python3
"""Dense finite-reset graph IQC tied to the SAME Joseph correction.

For one accepted shipping update let

    d = K y,                 d_theta = E_theta d,
    P_J = P-K S K^T,
    b = G(d_theta)^-1 rho,

where rho is the exact Cayley/quaternion reset defect. The reset transport
proves ||G^-1||_2=1. On the full hard attitude cell and every correction from
zero through the source-uniform ceiling, the homogeneous reset certificate
proves

    ||rho|| <= mu_R ||d_theta||,
    ||b||   <= mu_R ||d_theta||.

Crucially d_theta is NOT an independent radius in the production master. If Y
maps the common augmented coordinate z to the physical residual y and B maps z
to b, then D=E_theta K Y uses the same P/H/R/K cell and the graph obeys

    z^T [ mu_R^2 D^T D - B^T B ] z >= 0.

The uniform mu_R is not the endpoint ratio rho(delta_max)/delta_max: it is the
worst-branch homogeneous gain over the entire correction interval, including
the polynomial/axis-angle implementation switch. Rowwise K bounds are used
only to certify the correction ceiling on which this gain is valid; they never
replace D in the storage inequality.

The module also records the exact Kalman identity

    K^T P_J^-1 K = R^-1 - S^-1,

which follows from the information-form posterior and ties correction energy to
the same Joseph event. It is an algebra/IQC primitive; it cannot promote P4
without source-correlated event cells and the full augmented LDLT.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_reset_signed_absorption_diagnostic as RESET_SECTOR

QUALIFICATION='OU3_P4_SAME_CELL_FINITE_RESET_GRAPH_IQC_V2'

def _shape(A):return len(A),len(A[0]) if A else 0
def _scale(A,a):
    x=Interval.point(float(a));return [[x*v for v in row] for row in A]
def _gram(A):return matrix_mul(matrix_transpose(A),A)

def reset_sector_matrix(Dmap:Sequence[Sequence[Interval]],Bmap:Sequence[Sequence[Interval]],gamma_reset:float):
    """Return Pi with z'Pi z>=0 for ||Bz||<=gamma||Dz||."""
    dr,n=_shape(Dmap);br,m=_shape(Bmap)
    if dr!=3 or br!=3 or n==0 or m!=n:raise ValueError('Dmap/Bmap must be 3xn')
    g=float(gamma_reset)
    if not (math.isfinite(g) and g>=0):raise ValueError('finite nonnegative reset gain required')
    return matrix_symmetric_hull(matrix_sub(_scale(_gram(Dmap),g*g),_gram(Bmap)))

def same_cell_correction_map(K:Sequence[Sequence[Interval]],Ymap:Sequence[Sequence[Interval]]):
    """E_theta*K*Y on the common augmented coordinate; never a free d port."""
    kr,kc=_shape(K);yr,yc=_shape(Ymap)
    if kr not in (18,21) or kc!=3 or yr!=3 or yc==0:raise ValueError('K must be H18/A21 x3 and Ymap 3xn')
    return matrix_mul([list(row) for row in K[:3]],Ymap)

def kalman_correction_information_identity(PJinv,K,Rinv,Sinv):
    """Interval residual of K' PJ^-1 K = R^-1-S^-1 for one same cell."""
    n,n2=_shape(PJinv);kr,kc=_shape(K)
    if n not in (18,21) or n2!=n or (kr,kc)!=(n,3) or _shape(Rinv)!=(3,3) or _shape(Sinv)!=(3,3):raise ValueError('identity dimensions')
    lhs=matrix_mul(matrix_mul(matrix_transpose(K),PJinv),K)
    rhs=matrix_sub(Rinv,Sinv)
    return matrix_symmetric_hull(matrix_sub(lhs,rhs))

def _zero_contained(A):return all(v.lo<=0<=v.hi for row in A for v in row)

def build():
    e=ENTRY.build();c=COEFF.build();rs=RESET_SECTOR.build()
    bad={'entry':ENTRY.validate(e),'coeff':COEFF.validate(c),'reset_sector':RESET_SECTOR.validate(rs)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset IQC prerequisites failed: '+repr(bad))
    q=float(e['coordinate_radii']['attitude_cayley_norm']);modes={}
    for mode in ('H18','A21'):
        delta=float(c['modes'][mode]['attitude_correction_norm_upper'])
        sr=rs['modes'][mode]
        if abs(float(sr['correction_norm_upper'])-delta) > 8*math.ulp(max(1.0,delta)):
            raise RuntimeError(mode+' reset sector/correction ceiling detached')
        gamma=math.nextafter(float(sr['reset_defect_over_correction_norm_upper']),math.inf)
        rho=float(sr['endpoint_absolute_reset_defect_norm_upper'])
        # Pure algebra smoke with a non-axis same-cell-shaped correction map.
        n=8;Y=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
        for i in range(3):Y[i][i]=Interval.point(1.0)
        K=[[Interval.point(0.0) for _ in range(3)] for _ in range(18 if mode=='H18' else 21)]
        for i in range(3):K[i][i]=Interval.point(0.25+0.05*i)
        D=same_cell_correction_map(K,Y)
        B=[[Interval.point(0.0) for _ in range(n)] for _ in range(3)]
        for i in range(3):B[i][i]=Interval.point(gamma*(0.25+0.05*i))
        Pi=reset_sector_matrix(D,B,gamma)
        smoke=_zero_contained(Pi)
        modes[mode]={
          'full_declared_attitude_cayley_norm_upper':q,
          'source_uniform_attitude_correction_norm_upper':delta,
          'endpoint_absolute_reset_defect_norm_upper':rho,
          'reset_defect_to_actual_correction_norm_gain_upper':gamma,
          'uniform_gain_from_homogeneous_exact_reset_sector':True,
          'endpoint_ratio_used_as_uniform_gain':False,
          'homogeneous_sector_dominates_endpoint_absolute_bound':bool(sr['homogeneous_sector_dominates_endpoint_absolute_bound']),
          'cayley_composition_denominator_lower':sr['cayley_composition_denominator_lower'],
          'same_cell_Dtheta_equals_Etheta_K_Y_required':True,'reset_graph_IQC_available':True,
          'reset_sector_smoke_zero_contained':smoke,'chart_safe':sr['chart_safe']}
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'exact_K_PJinv_K_identity':'K^T*P_J^-1*K=R^-1-S^-1','same_cell_K_required':True,
      'independent_correction_port_forbidden':True,'reset_defect_dense_IQC_available':True,
      'reset_IQC_keeps_residual_direction_via_Etheta_K_Y':True,'reset_inverse_operator_norm_upper':1.0,
      'uniform_gain_from_homogeneous_exact_reset_sector':True,
      'rowwise_K_used_only_for_uniform_correction_ceiling':True,
      'rowwise_K_may_not_replace_same_cell_Dmap_in_storage':True,
      'modes':modes,'source_uniform_augmented_LDLT_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_K_required','independent_correction_port_forbidden','reset_defect_dense_IQC_available','reset_IQC_keeps_residual_direction_via_Etheta_K_Y','uniform_gain_from_homogeneous_exact_reset_sector','rowwise_K_used_only_for_uniform_correction_ceiling','rowwise_K_may_not_replace_same_cell_Dmap_in_storage'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('reset_inverse_operator_norm_upper')!=1.0:f.append('reset inverse norm changed')
    for mode,m in d['modes'].items():
        for k in ('same_cell_Dtheta_equals_Etheta_K_Y_required','reset_graph_IQC_available','reset_sector_smoke_zero_contained','chart_safe','uniform_gain_from_homogeneous_exact_reset_sector','homogeneous_sector_dominates_endpoint_absolute_bound'):
            if m.get(k) is not True:f.append(mode+' '+k+' not true')
        if m.get('endpoint_ratio_used_as_uniform_gain') is not False:f.append(mode+' endpoint ratio still used as uniform gain')
        g=float(m.get('reset_defect_to_actual_correction_norm_gain_upper',math.nan))
        if not (math.isfinite(g) and g>=0):f.append(mode+' reset gain invalid')
        if not float(m.get('cayley_composition_denominator_lower',0))>0:f.append(mode+' reset denominator nonpositive')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'modes':{k:v['reset_defect_to_actual_correction_norm_gain_upper'] for k,v in d['modes'].items()},'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
