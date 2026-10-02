# ruff: noqa: F401, F811
"""Probe PSD-factor inverse in rotated accelerometer innovation frame."""
from __future__ import annotations
import math,numpy as np
from .attitude_ball_cover import lattice_cover
from .accel_geometry_cell import rotation_cell_centered
from .rotated_accel_geometry import transformed_h_cell,rotate_measurement_covariance_isotropic
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .psd_factor_innovation import factor_innovation_probe
RHO=18.60665

def exact_force_probe():
 cells,_=lattice_cover();_,_,_,Racc,*_=one_sample_boxes();Racc=rotate_measurement_covariance_isotropic(Racc)
 if Racc is None:return {"verified":False,"reason":"Racc not rotation invariant"}
 L=covariance_factor_interval();bad=[];worst=None
 # Exact local-body force at worst magnitude; direction z. This isolates whether
 # the AW cancellation restores useful attitude margin.
 fb=np.array([0.,0.,RHO])
 for i,c in enumerate(cells):
  Rm,Rr=rotation_cell_centered(c.center,c.radius)
  H=transformed_h_cell(Rm,Rr,fb,np.zeros(3))
  z=factor_innovation_probe(H,L,Racc);row={"cell":i,**z}
  if worst is None or row["residual"]>worst["residual"]:worst=row
  if not z["verified"]:bad.append(row)
 return {"verified":not bad,"worst":worst,"unverified_count":len(bad),
         "aw_block_exact_identity":True}

def full_force_ball_probe():
 cells,_=lattice_cover();_,_,_,Racc,*_=one_sample_boxes();Racc=rotate_measurement_covariance_isotropic(Racc)
 L=covariance_factor_interval();bad=[];worst=None
 # Full body-frame force ball outer component interval, but now independent of
 # attitude rotation because force and innovation use the same body frame.
 for i,c in enumerate(cells):
  Rm,Rr=rotation_cell_centered(c.center,c.radius)
  H=transformed_h_cell(Rm,Rr,np.zeros(3),np.full(3,RHO))
  z=factor_innovation_probe(H,L,Racc);row={"cell":i,**z}
  if worst is None or row["residual"]>worst["residual"]:worst=row
  if not z["verified"]:bad.append(row)
 return {"verified":not bad,"worst":worst,"unverified_count":len(bad),
         "aw_block_exact_identity":True}
def diagnostic():return {"exact_force":exact_force_probe(),"full_force_ball":full_force_ball_probe()}
