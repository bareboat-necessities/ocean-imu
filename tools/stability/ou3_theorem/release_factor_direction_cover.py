# ruff: noqa: F401, F811
"""PSD-factor innovation cover using force magnitude/direction cones."""
from __future__ import annotations
import math,heapq,numpy as np
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .accel_geometry_cell import rotation_cell_from_ball
from .force_direction_cone import octahedral_direction_cones,split_radial,accel_h_force_cone
from .psd_factor_innovation import factor_innovation_probe

def attempt(c):
 _,_,_,Racc,*_=one_sample_boxes();Rm,Rr=rotation_cell_from_ball(math.radians(6.9))
 H=accel_h_force_cone(c,Rm,Rr)
 return factor_innovation_probe(H,covariance_factor_interval(),Racc)

def cover(max_radial_depth=12):
 cones=octahedral_direction_cones();good=[];bad=[];heap=[];counter=0
 for c in cones:heapq.heappush(heap,(-(c.rho_hi-c.rho_lo),counter,0,c));counter+=1
 while heap:
  _,_,d,c=heapq.heappop(heap);z=attempt(c)
  if z["verified"]:good.append((c,z));continue
  if d>=max_radial_depth:bad.append((c,z));continue
  a,b=split_radial(c)
  for x in (a,b):heapq.heappush(heap,(-(x.rho_hi-x.rho_lo),counter,d+1,x));counter+=1
 return {"verified":not bad,"leaf_count":len(good),"unresolved_count":len(bad),
  "worst_verified_residual":max((z["residual"] for _,z in good),default=None),
  "first_unresolved":bad[0][1] if bad else None,
  "direction_cones":6,"direction_half_angle_deg":math.degrees(cones[0].half_angle)}
