"""Prospective latent-AW mean and weighted INJ-axis defect recurrence.

The operation iterator preserves shipping order.  Vector mode is preferred:
gain/residual products are formed before taking a norm, retaining directional
cancellation that the older scalar triangle diagnostic discarded.
"""
from __future__ import annotations
import math
import numpy as np
G=9.80665

def axis_delta_from_aw(A,phi):
 """Axis change for w=aw-g -> phi aw-g from current ||aw|| upper bound."""
 df=(1-phi)*A;r=max(0.,G-A)
 if r<=0 or df>=2*r:return math.pi
 return 2*math.asin(min(1.,df/(2*r)))

def aw_mean_step(A,phi,Kaw_norm,residual_norm,KawS_norm=0.,Smean_norm=0.,
                 KawMag_norm=0.,mag_residual_norm=0.):
 """Scalar outer recurrence retained only as an explicit relaxation."""
 pred=phi*A
 acc=Kaw_norm*residual_norm
 sint=KawS_norm*Smean_norm
 mag=KawMag_norm*mag_residual_norm
 return {"A_pred":pred,"acc_increment":acc,"S_increment":sint,
         "mag_increment":mag,"A_next":pred+acc+sint+mag,
         "relaxation":"triangle_norm"}

def vector_chronology_step(aw,phi,operations):
 """Iterate literal AW mean through ordered accepted corrections.

 operations is an ordered iterable of (kind,K,residual).  K and residual must
 come from the same history leaf and same pre-correction boundary.
 """
 x=phi*np.asarray(aw,float)
 trace=[{"kind":"prediction","aw":x.copy()}]
 for kind,K,r in operations:
  if kind not in ("acc","S","mag"):
   raise ValueError("unknown AW correction kind "+str(kind))
  x=x+np.asarray(K,float)@np.asarray(r,float)
  trace.append({"kind":kind,"aw":x.copy()})
 return {"aw_next":x,"A_next":float(np.linalg.norm(x)),"trace":trace,
         "products_formed_before_norm":True}

def weighted_inj_defect(b,c,A,phi):
 d=axis_delta_from_aw(A,phi)
 return {"delta_rad":d,"delta_deg":math.degrees(d),
  "D":abs(b)*(math.sin(d) if d<math.pi else 1.)+2*abs(c)*math.sin(d/2)}

def coefficient_step(A,phi,b,c,operations=None,**scalar):
 """One (A_k,b_k,c_k,D_k) release step.

 b,c are the carried INJ coefficients of the block being rebased.  Covariance
 chronology supplies their next values separately after the correction/reset;
 this function prices only the exact axis rebase caused by OU prediction.
 """
 z=weighted_inj_defect(b,c,A,phi)
 if operations is not None:
  if "aw_vector" not in scalar:raise ValueError("vector chronology needs aw_vector")
  mean=vector_chronology_step(scalar["aw_vector"],phi,operations)
 else:
  mean=aw_mean_step(A,phi,scalar.get("Kaw_norm",0.),scalar.get("residual_norm",0.),
                    scalar.get("KawS_norm",0.),scalar.get("Smean_norm",0.),
                    scalar.get("KawMag_norm",0.),scalar.get("mag_residual_norm",0.))
 return {**mean,**z,"A_k":A,"b_k":b,"c_k":c,"D_k":z["D"]}

def step(A,phi,b,c,Kaw_norm,residual_norm,KawS_norm=0.,Smean_norm=0.):
 return coefficient_step(A,phi,b,c,Kaw_norm=Kaw_norm,residual_norm=residual_norm,
                         KawS_norm=KawS_norm,Smean_norm=Smean_norm)
