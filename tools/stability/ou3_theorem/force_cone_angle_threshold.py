"""Solve required force-direction cone angle for PSD-factor inverse residual."""
from __future__ import annotations
import math,numpy as np
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .accel_geometry_cell import rotation_cell_from_ball
from .force_direction_cone import ForceCone,accel_h_force_cone
from .psd_factor_innovation import factor_innovation_probe

RHO=18.60665
AXIS=np.array([0.,0.,1.])

def residual_for(alpha_deg,rho_lo=RHO,rho_hi=RHO):
 _,_,_,Racc,*_=one_sample_boxes();Rm,Rr=rotation_cell_from_ball(math.radians(6.9))
 c=ForceCone(rho_lo,rho_hi,AXIS,math.radians(alpha_deg))
 H=accel_h_force_cone(c,Rm,Rr)
 z=factor_innovation_probe(H,covariance_factor_interval(),Racc)
 return {"alpha_deg":alpha_deg,**z}

def sweep():
 angles=[60,45,30,20,15,10,7.5,5,3,2,1,.5,.25,.1,.05,.02,.01,0.]
 return [residual_for(a) for a in angles]

def solve_alpha(tol_deg=1e-5):
 rows=sweep()
 good=[x for x in rows if x["residual"]<1]
 if not good:return {"verified":False,"rows":rows,"reason":"no passing angle"}
 # largest passing angle, then bracket against nearest larger failing sample.
 g=max(good,key=lambda x:x["alpha_deg"]);lo=g["alpha_deg"]
 fails=[x for x in rows if x["alpha_deg"]>lo and x["residual"]>=1]
 if not fails:return {"verified":True,"alpha_star_lower_deg":lo,"alpha_star_upper_deg":None,"rows":rows}
 hi=min(fails,key=lambda x:x["alpha_deg"])["alpha_deg"]
 for _ in range(80):
  if hi-lo<=tol_deg:break
  m=(lo+hi)/2;z=residual_for(m)
  if z["residual"]<1:lo=m
  else:hi=m
 return {"verified":True,"alpha_star_lower_deg":lo,"alpha_star_upper_deg":hi,
         "residual_at_lower":residual_for(lo)["residual"],
         "residual_at_upper":residual_for(hi)["residual"],"rows":rows}
