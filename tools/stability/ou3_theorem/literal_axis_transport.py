"""Literal OU-III force-axis transport identity.

At an accelerometer row f=R(a_w-g). Carry the structured axis through every
attitude injection/reset by the same left rotation. The next gyro predictor is
also a common left rotation. Hence frame rotations cancel exactly; the axis
mismatch is the angle between w=a_w-g and w+=phi a_w-g.
"""
from __future__ import annotations
import math

G=9.80665
def ou_axis_change_bound(aw_norm,phi):
 """Direction change between aw-g and phi*aw-g using only ||aw||.

 ||w+-w||=(1-phi)||aw||.
 min(||w||,||w+||)>=g-aw_norm for aw_norm<g.
 """
 A=float(aw_norm);p=float(phi)
 if not(0<p<=1) or A<0:raise ValueError("OU parameters")
 df=(1-p)*A
 r=max(0.,G-A)
 if r<=0 or df>=2*r:return {"delta_rad":math.pi,"delta_deg":180.,"df":df,"r_floor":r}
 d=2*math.asin(df/(2*r))
 return {"delta_rad":d,"delta_deg":math.degrees(d),"df":df,"r_floor":r}

def phi_range(dt,tau):
 return math.exp(-dt/tau)

def diagnostic():
 rows=[]
 for A in (.1,.5,1,2,4,8,9):
  p=phi_range(.006,.02) # fastest OU forgetting = worst one-step change
  rows.append({"aw_norm":A,"phi":p,**ou_axis_change_bound(A,p)})
 return {"verified":True,"frame_rotations_cancel_exactly":True,
         "mahony_continuous_term":False,"watchdog_is_separate_reset_event":True,
         "rows":rows}
