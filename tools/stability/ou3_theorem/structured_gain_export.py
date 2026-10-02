"""Structured gain export for AW rows without entrywise K intervals."""
from __future__ import annotations
from dataclasses import dataclass
from .inj_block_algebra import INJ
from .structured_riccati_closure import inverse_coefficients

def acc_aw_gain(P_aw_th:INJ,P_aw_aw:INJ,P_aw_ba:INJ|None,
                Jatt_T:INJ,Jaw_T:INJ,S:INJ,use_ba:bool):
 """K_aw=(P_aw,th Jatt' + P_aw,aw Jaw' + P_aw,ba) S^-1."""
 num=P_aw_th.mul(Jatt_T)+P_aw_aw.mul(Jaw_T)
 if use_ba and P_aw_ba is not None:num=num+P_aw_ba
 return num.mul(inverse_coefficients(S))

def inj_spectral_norm(x:INJ):
 """Exact 2-norm: parallel scalar i+n; transverse complex magnitude sqrt(i^2+j^2)."""
 return max(abs(x.i+x.n),(x.i*x.i+x.j*x.j)**.5)

def s_aw_gain_scalar(p_aw_s,p_ss,r_s):
 """When LIN blocks are scalar*I: K_aw,S = k I exactly."""
 den=p_ss+r_s
 if not den>0:raise ArithmeticError("S pseudo innovation SPD")
 return p_aw_s/den

def export(P_aw_th,P_aw_aw,P_aw_ba,Jatt_T,Jaw_T,S,use_ba,p_aw_s,p_ss,r_s):
 ka=acc_aw_gain(P_aw_th,P_aw_aw,P_aw_ba,Jatt_T,Jaw_T,S,use_ba)
 ks=s_aw_gain_scalar(p_aw_s,p_ss,r_s)
 return {"Kaw":ka,"Kaw_norm":inj_spectral_norm(ka),
         "KawS_scalar":ks,"KawS_norm":abs(ks),
         "gain_intervalized":False}
