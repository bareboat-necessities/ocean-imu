"""Joint attitude-cover x force-cone factor residual threshold."""
from __future__ import annotations
import math,numpy as np
from .attitude_ball_cover import lattice_cover
from .accel_geometry_cell import rotation_cell_centered
from .force_direction_cone import ForceCone,accel_h_force_cone
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .psd_factor_innovation import factor_innovation_probe
RHO=18.60665;AX=np.array([0.,0.,1.])

def worst(alpha_deg,drho):
 cells,_=lattice_cover();_,_,_,Racc,*_=one_sample_boxes();L=covariance_factor_interval()
 # Test outer radial shell centered at RHO-drho/2; this maximizes angular term
 # and is the decisive radial cell. Inner shells have smaller rho_hi.
 lo=max(0.,RHO-drho);c=ForceCone(lo,RHO,AX,math.radians(alpha_deg))
 w=None
 for i,a in enumerate(cells):
  Rm,Rr=rotation_cell_centered(a.center,a.radius)
  z=factor_innovation_probe(accel_h_force_cone(c,Rm,Rr),L,Racc)
  if w is None or z["residual"]>w["residual"]:w={"cell":i,**z}
 return {"alpha_deg":alpha_deg,"drho":drho,**w}

def angle_threshold(drho=0.,tol=.001):
 # alpha=0 must pass for this radial width.
 z0=worst(0.,drho)
 if z0["residual"]>=1:return {"verified":False,"reason":"radial width alone fails","zero_angle":z0}
 lo=0.;hi=10.
 while worst(hi,drho)["residual"]<1 and hi<90:lo=hi;hi*=2
 for _ in range(40):
  if hi-lo<=tol:break
  m=(lo+hi)/2
  if worst(m,drho)["residual"]<1:lo=m
  else:hi=m
 return {"verified":True,"drho":drho,"alpha_lower_deg":lo,"alpha_upper_deg":hi,
  "lower":worst(lo,drho),"upper":worst(hi,drho)}

def radial_threshold(tol=.001):
 lo=0.;hi=RHO
 if worst(0.,hi)["residual"]<1:return {"verified":True,"drho_lower":hi,"whole_radius":True}
 for _ in range(50):
  if hi-lo<=tol:break
  m=(lo+hi)/2
  if worst(0.,m)["residual"]<1:lo=m
  else:hi=m
 return {"verified":True,"drho_lower":lo,"drho_upper":hi,
  "lower":worst(0.,lo),"upper":worst(0.,hi)}

def diagnostic():
 return {"radial":radial_threshold(),"angle_at_zero_radial":angle_threshold(0.)}
