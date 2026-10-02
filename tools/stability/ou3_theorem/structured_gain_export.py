"""Structured gain export for AW rows without entrywise K intervals.

All blocks are expressed in the common force-axis INJ algebra.  This mirrors the
literal shipping PC^T construction: attitude + AW + optional active BA +
optional lever-arm gyro-bias columns, followed by the same innovation inverse.
"""
from __future__ import annotations
from .inj_block_algebra import INJ
from .structured_riccati_closure import inverse_coefficients

ZERO=INJ(0.,0.,0.)

def selected_aw_gain(terms:list[tuple[INJ,INJ]],S:INJ):
 """Return (sum_j P_aw,j H_j^T) S^-1 in the shared INJ basis."""
 num=ZERO
 for P,H_T in terms:num=num+P.mul(H_T)
 return num.mul(inverse_coefficients(S))

def acc_aw_gain(P_aw_th:INJ,P_aw_aw:INJ,P_aw_ba:INJ|None,
                Jatt_T:INJ,Jaw_T:INJ,S:INJ,use_ba:bool,
                P_aw_bg:INJ|None=None,Jbg_T:INJ|None=None,
                use_lever_bg:bool=False):
 """Literal AW row of accelerometer gain.

 K_aw=(P_aw,th Jatt' + P_aw,aw Jaw' + [active BA]P_aw,ba
       + [lever]P_aw,bg Jbg') S_acc^-1.
 """
 terms=[(P_aw_th,Jatt_T),(P_aw_aw,Jaw_T)]
 if use_ba and P_aw_ba is not None:terms.append((P_aw_ba,INJ(1.,0.,0.)))
 if use_lever_bg:
  if P_aw_bg is None or Jbg_T is None:
   raise ValueError("lever-arm gain requires P_aw_bg and Jbg_T")
  terms.append((P_aw_bg,Jbg_T))
 return selected_aw_gain(terms,S)

def inj_spectral_norm(x:INJ):
 """Exact 2-norm: parallel scalar i+n; transverse complex magnitude sqrt(i^2+j^2)."""
 return max(abs(x.i+x.n),(x.i*x.i+x.j*x.j)**.5)

def s_aw_gain(P_aw_s:INJ,P_ss:INJ,R_s:INJ):
 """Literal K_aw,S=P_aw,S(P_SS+R_S)^-1; no isotropy assumption."""
 return P_aw_s.mul(inverse_coefficients(P_ss+R_s))

def s_aw_gain_scalar(p_aw_s,p_ss,r_s):
 """Compatibility helper for the special scalar*I case."""
 return s_aw_gain(INJ(p_aw_s,0,0),INJ(p_ss,0,0),INJ(r_s,0,0)).i

def mag_aw_gain(P_aw_th:INJ,Jmag_T:INJ,S_mag:INJ):
 """Literal AW row of an accepted magnetic correction."""
 return selected_aw_gain([(P_aw_th,Jmag_T)],S_mag)

def export(P_aw_th,P_aw_aw,P_aw_ba,Jatt_T,Jaw_T,S,use_ba,
           P_aw_s,P_ss,R_s,P_aw_bg=None,Jbg_T=None,use_lever_bg=False):
 ka=acc_aw_gain(P_aw_th,P_aw_aw,P_aw_ba,Jatt_T,Jaw_T,S,use_ba,
                P_aw_bg,Jbg_T,use_lever_bg)
 ks=s_aw_gain(P_aw_s,P_ss,R_s)
 return {"Kaw":ka,"Kaw_norm":inj_spectral_norm(ka),
         "KawS":ks,"KawS_norm":inj_spectral_norm(ks),
         "gain_intervalized":False,
         "lever_bg_column_retained":bool(use_lever_bg),
         "S_gain_isotropy_assumed":False}
