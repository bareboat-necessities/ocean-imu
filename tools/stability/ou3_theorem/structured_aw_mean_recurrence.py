"""Prospective latent-AW mean and weighted INJ-axis defect recurrence."""
from __future__ import annotations
import math
G=9.80665

def axis_delta_from_aw(A,phi):
 """Axis change for w=aw-g -> phi aw-g from current ||aw|| upper bound."""
 df=(1-phi)*A;r=max(0.,G-A)
 if r<=0 or df>=2*r:return math.pi
 return 2*math.asin(min(1.,df/(2*r)))

def aw_mean_step(A,phi,Kaw_norm,residual_norm,KawS_norm=0.,Smean_norm=0.):
 """Literal prospective norm recurrence: prediction, acc correction, optional S correction."""
 pred=phi*A
 acc=Kaw_norm*residual_norm
 sint=KawS_norm*Smean_norm
 return {"A_pred":pred,"acc_increment":acc,"S_increment":sint,
         "A_next":pred+acc+sint}

def weighted_inj_defect(b,c,A,phi):
 d=axis_delta_from_aw(A,phi)
 return {"delta_rad":d,"delta_deg":math.degrees(d),
  "D":abs(b)*(math.sin(d) if d<math.pi else 1.)+2*abs(c)*math.sin(d/2)}

def step(A,phi,b,c,Kaw_norm,residual_norm,KawS_norm=0.,Smean_norm=0.):
 z=weighted_inj_defect(b,c,A,phi)
 return {**aw_mean_step(A,phi,Kaw_norm,residual_norm,KawS_norm,Smean_norm),**z,
  "same_history_inputs_required":["Kaw","acc residual","KawS","S mean","INJ b","INJ c"]}
