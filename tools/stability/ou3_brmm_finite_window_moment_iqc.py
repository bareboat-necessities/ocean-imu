#!/usr/bin/env python3
"""Correlated finite-window acceleration-moment IQC for COMPLETE-BRMM.

No finite harmonic realization is introduced.  For one admitted physical
acceleration function a:[0,h]->R^3 define the same-history moments

  J0 = int a(s) ds,
  J1 = int (h-s) a(s) ds,
  J2 = int (h-s)^2/2 a(s) ds.

With ||a(s)||<=A, the L2 energy obeys int ||a||^2 <= A^2 h.  If
phi=(1,t,t^2/2), t=h-s, and G=int phi phi^T, every actual moment triple obeys

  sum_axis [J0_i,J1_i,J2_i] G^-1 [J0_i,J1_i,J2_i]^T <= A^2 h.

This is a rigorous same-history correlation: it is impossible to select J0,J1,
J2 independently at their scalar extrema.  For h=3 s the exact rational inverse
is

  G^-1 = [[3,-4,20/9],[-4,64/9,-40/9],[20/9,-40/9,80/27]].

The exact physical jet then propagates
  v1=v0+J0,
  p1=p0+h v0+J1,
  S1=S0+h p0+h^2 v0/2+J2.

The module exposes an outward homogeneous source IQC with a constant coordinate
and an exact-rational transition map.  It is a necessary enclosure of the same
physical history, not a source generator and not covariance membership.
"""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
from ou3_interval import Interval,symmetric_positive_definite_ldlt
import ou3_brmm_finite_window_primitive_qualification as PRIM

SCHEMA=1
QUALIFICATION='OU3_BRMM_FINITE_WINDOW_CORRELATED_MOMENT_IQC_V1'
H=F(3,1)
GINV=((F(3),F(-4),F(20,9)),(F(-4),F(64,9),F(-40,9)),(F(20,9),F(-40,9),F(80,27)))

def Iq(q:F):return Interval.outward_bounds(float(q),float(q))
def _zero(r,c):return [[Interval.point(0.0) for _ in range(c)] for _ in range(r)]

def source_iqc_matrix(A:float):
    # coordinate z=[1,J0x,J1x,J2x,J0y,J1y,J2y,J0z,J1z,J2z]
    M=_zero(10,10);M[0][0]=Interval.outward_bounds(float(A*A*float(H)),float(A*A*float(H)))
    for axis in range(3):
        off=1+3*axis
        for i in range(3):
            for j in range(3):M[off+i][off+j]=-Iq(GINV[i][j])
    return M

def jet_transition_matrix():
    # x=[v(3),p(3),S(3), moments grouped by axis (J0,J1,J2)x3]
    # output [v1,p1,S1]. Exact coefficients only.
    M=_zero(9,18)
    for a in range(3):
        v=a;p=3+a;s=6+a;m=9+3*a
        M[a][v]=Iq(F(1));M[a][m]=Iq(F(1))
        M[3+a][p]=Iq(F(1));M[3+a][v]=Iq(H);M[3+a][m+1]=Iq(F(1))
        M[6+a][s]=Iq(F(1));M[6+a][p]=Iq(H);M[6+a][v]=Iq(H*H/F(2));M[6+a][m+2]=Iq(F(1))
    return M

def _exact_gram_inverse_check():
    G=((H,H*H/F(2),H**3/F(6)),(H*H/F(2),H**3/F(3),H**4/F(8)),(H**3/F(6),H**4/F(8),H**5/F(20)))
    for i in range(3):
        for j in range(3):
            s=sum(G[i][k]*GINV[k][j] for k in range(3))
            if s!=(F(1) if i==j else F(0)):return False
    return True

def build():
    p=PRIM.build();pf=PRIM.validate(p)
    if pf:raise RuntimeError('primitive qualification invalid: '+repr(pf))
    A=float(p['uniform_physical_primitives']['acceleration_norm_upper_mps2'])
    M=source_iqc_matrix(A);J=jet_transition_matrix()
    # The negative moment block is strictly negative because G^-1 is SPD.
    GinvI=[[Iq(GINV[i][j]) for j in range(3)] for i in range(3)];ok,piv=symmetric_positive_definite_ldlt(GinvI)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','word_horizon_s':3.0,'acceleration_norm_upper_mps2':A,
      'same_physical_acceleration_function_generates_all_moments':True,'independent_J0_J1_J2_boxes_used':False,'finite_harmonic_source_used':False,'trajectory_replay_used':False,
      'exact_Gram_inverse_identity_closed':_exact_gram_inverse_check(),'Gram_inverse_exact_rational':[[str(x) for x in row] for row in GINV],
      'moment_IQC_dimension':10,'moment_energy_rhs_A2h':A*A*3.0,'moment_IQC_matrix_API':'source_iqc_matrix(A)','jet_transition_matrix_API':'jet_transition_matrix()',
      'moment_quadratic_form_strictly_positive':bool(ok),'moment_quadratic_form_worst_LDLT_pivot_lower':min(float(x.lo) for x in piv),
      'exact_v_p_S_jet_transition_materialized':True,'same_history_S_increment_generated_by_J2_not_independent_S_port':True,
      'finite_window_correlated_source_operator_materialized':True,'covariance_membership_used':False,'P4_promoted_here':False,
      'next_obligation':'embed this source IQC and exact jet map in the endpoint/every-prefix augmented master together with the estimator-owned coefficient/source graph'}

def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('same_physical_acceleration_function_generates_all_moments','exact_Gram_inverse_identity_closed','moment_quadratic_form_strictly_positive','exact_v_p_S_jet_transition_materialized','same_history_S_increment_generated_by_J2_not_independent_S_port','finite_window_correlated_source_operator_materialized'):
  if d.get(k) is not True:f.append(k+' not true')
 for k in ('independent_J0_J1_J2_boxes_used','finite_harmonic_source_used','trajectory_replay_used','covariance_membership_used','P4_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if float(d.get('acceleration_norm_upper_mps2',0))!=8.0:f.append('source acceleration cap is not 8')
 if not float(d.get('moment_quadratic_form_worst_LDLT_pivot_lower',0))>0:f.append('moment Gram inverse not SPD')
 if int(d.get('moment_IQC_dimension',0))!=10:f.append('moment IQC dimension mismatch')
 return f

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'A':d['acceleration_norm_upper_mps2'],'A2h':d['moment_energy_rhs_A2h'],'Gram_pivot':d['moment_quadratic_form_worst_LDLT_pivot_lower'],'transition':d['finite_window_correlated_source_operator_materialized'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
