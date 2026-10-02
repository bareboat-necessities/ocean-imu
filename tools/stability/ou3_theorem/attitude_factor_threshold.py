"""Solve attitude-ball radius threshold for factor innovation with exact force."""
from __future__ import annotations
import math,numpy as np
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .accel_geometry_cell import rotation_cell_from_ball
from .force_direction_cone import ForceCone,accel_h_force_cone
from .psd_factor_innovation import factor_innovation_probe
RHO=18.60665
def residual_attitude(beta_deg):
 _,_,_,Racc,*_=one_sample_boxes();Rm,Rr=rotation_cell_from_ball(math.radians(beta_deg))
 c=ForceCone(RHO,RHO,np.array([0.,0.,1.]),0.)
 return {"beta_deg":beta_deg,**factor_innovation_probe(accel_h_force_cone(c,Rm,Rr),covariance_factor_interval(),Racc)}
def sweep():
 return [residual_attitude(x) for x in [6.9,5,3,2,1,.5,.25,.1,.05,.02,.01,0]]
