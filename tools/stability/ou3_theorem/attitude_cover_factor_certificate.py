"""Run PSD-factor inverse certificate over the complete attitude cover."""
from __future__ import annotations
import numpy as np
from .attitude_ball_cover import lattice_cover,verify_cover_metadata
from .accel_geometry_cell import rotation_cell_centered
from .force_direction_cone import ForceCone,accel_h_force_cone
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .psd_factor_innovation import factor_innovation_probe
RHO=18.60665

def run():
 cells,meta=lattice_cover();_,_,_,Racc,*_=one_sample_boxes();L=covariance_factor_interval()
 # First decisive test requested: exact worst-magnitude force along a fixed
 # direction, across every attitude cover cell. Force cones are added only
 # after the attitude cover itself verifies.
 f=ForceCone(RHO,RHO,np.array([0.,0.,1.]),0.)
 worst=None;bad=[]
 for i,c in enumerate(cells):
  Rm,Rr=rotation_cell_centered(c.center,c.radius)
  z=factor_innovation_probe(accel_h_force_cone(f,Rm,Rr),L,Racc)
  row={"cell":i,"center_norm_deg":float(np.degrees(np.linalg.norm(c.center))),**z}
  if worst is None or row["residual"]>worst["residual"]:worst=row
  if not z["verified"]:bad.append(row)
 return {"verified":not bad,"cover":verify_cover_metadata(),"worst":worst,
         "unverified_count":len(bad),"first_unverified":bad[0] if bad else None}
