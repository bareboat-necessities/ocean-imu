"""Kernel-aware linked quotient over one temporal source domain.

Uses one augmented quadratic form in z=[e,u].  Linear source constraints remain
attached to the same u. Nonlinear MARINE/magnetic constraints are mandatory
side certificates; absence fails closed.
"""
from __future__ import annotations
import math
from .interval_riccati import symmetric_interval_gershgorin
from .rank_loss_interval_factor import IMat,verified_inverse,matmul,transpose,add,scale

def _sub(a,b):return add(a,scale(b,-1))

def schur_homogeneous_lower(iq,source_metric:IMat):
 """Verified lower bound after linked source elimination.

 source_metric is a SAME-history quadratic certificate R for u'R u<=1. This
 function does not manufacture R from independent maxima.
 """
 if source_metric.shape!=iq.C.shape:raise ValueError("source metric shape")
 # For lambda>=0, q + lambda(1-u'Ru) has homogeneous Schur
 # A-B(lambda R-C)^-1 B'. Optimize lambda outside this primitive.
 return {"A":iq.A,"B":iq.B,"C":iq.C,"R":source_metric,"verified":False}

def evaluate_lambda(iq,R:IMat,lam:float):
 if not(lam>0 and math.isfinite(lam)):raise ValueError("positive lambda")
 G=_sub(scale(R,lam),iq.C);Gi,cert=verified_inverse(G)
 if not cert.get("verified"):return {"verified":False,"reason":"source Schur inverse"}
 L=_sub(iq.A,matmul(matmul(iq.B,Gi),transpose(iq.B)))
 lo,hi=symmetric_interval_gershgorin(L.mid,L.rad)
 return {"verified":lo>0,"dissipation_lower":max(0.,lo),"dissipation_upper":hi,
         "lambda":lam,"source_constant":lam,"inverse":cert}

def optimize_lambda_grid(iq,R:IMat,lams):
 best=None
 for x in lams:
  z=evaluate_lambda(iq,R,float(x))
  if z.get("verified") and (best is None or z["source_constant"]/z["dissipation_lower"]<best["quotient_upper"]):
   z["quotient_upper"]=z["source_constant"]/z["dissipation_lower"];best=z
 return best or {"verified":False,"reason":"no verified linked Schur point"}
