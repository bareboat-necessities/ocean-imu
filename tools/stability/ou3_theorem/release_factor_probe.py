"""Probe factor-structured goLive accelerometer innovation on captured domain."""
import math,numpy as np
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .accel_geometry_cell import accel_h_from_force_rotation,rotation_cell_from_ball
from .psd_factor_innovation import factor_innovation_probe

def probe():
 _,_,_,Racc,*_=one_sample_boxes()
 Rm,Rr=rotation_cell_from_ball(math.radians(6.9))
 H=accel_h_from_force_rotation(np.zeros(3),np.full(3,18.60665),Rm,Rr)
 return factor_innovation_probe(H,covariance_factor_interval(),Racc)
