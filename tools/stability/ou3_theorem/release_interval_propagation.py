# ruff: noqa: F401, F811
"""Constructive goLive->A21 covariance release enclosure over shipping ranges."""
from __future__ import annotations
import json,math
from pathlib import Path
from .golive_release_seed import covariance_interval
from .interval_riccati_21 import (shipping_prediction_intervals,shipping_acc_update_intervals,
 shipping_integral_update_intervals,shipping_mag_update_intervals,predict_covariance,
 verified_joseph_update,verified_joseph_update_psd,linked_covariance_update,aw_covariance_floor_event,spectral_box)
ROOT=Path(__file__).resolve().parents[3]
C=json.loads((ROOT/"tools/stability/ou3_theorem/constants.json").read_text())

def one_sample_boxes():
 sm=C["sensor_model"];m=C["marine_motion"];mag=C["magnetic_service"]
 F,Q=shipping_prediction_intervals(dt_min=sm["sample_period_min_s"],dt_max=sm["sample_period_max_s"],
  tau_min=.02,tau_max=12.,omega_max=m["Omega_max_rad_s"],tau_bacc=120.,
  gyro_white_density=.00135,gyro_bias_rw_density=math.sqrt(1e-11),
  aw_sigma_max=4.,accel_bias_drive_density=5e-4)
 Hacc,Racc=shipping_acc_update_intervals(9.80665+m["A_max_mps2"],sm["accel_nominal_std_mps2"],sm["accel_effective_std_max_mps2"])
 HS,RS=shipping_integral_update_intervals(math.sqrt(.15),10.)
 Hmag,Rmag=shipping_mag_update_intervals(mag["field_norm_max_uT"],.1,2.)
 return F,Q,Hacc,Racc,HS,RS,Hmag,Rmag

def covariance_release_image(max_steps=None):
 F,Q,Hacc,Racc,HS,RS,Hmag,Rmag=one_sample_boxes();P=covariance_interval();certs=[]
 # Certified release horizon 469 s. Use dt_min for maximum number of predictions.
 if max_steps is None:max_steps=math.ceil(469./C["sensor_model"]["sample_period_min_s"])
 for k in range(max_steps):
  P=predict_covariance(P,F,Q)
  # Literal accelerometer correction every qualified regular sample.
  P,c=verified_joseph_update(P,Hacc,Racc);certs.append(("acc",k,c))
  # Optional S/mag corrections are covariance reducing in pre-reset coordinates.
  # They are not inserted at an invented cadence here; aggregate magnetic
  # information is consumed by the homogeneous action certificate separately.
  # AW sync is an increase and must be retained. Reachable target <=sigma^2<=16.
  P=aw_covariance_floor_event(P,16.)
 return {"verified":True,"P":P,"spectral_box":spectral_box(P),"steps":max_steps,
         "innovation_certificates":certs[-3:],"S_schedule_enumerated":False,
         "mag_schedule_enumerated":False,
         "note":"S/mag covariance reductions omitted; acc literal every sample; reset transport still required for final promotion"}

def probe():
 try:
  z=covariance_release_image(1)
  return {"one_step_verified":z["verified"],"spectral_box":z["spectral_box"]}
 except Exception as e:return {"one_step_verified":False,"reason":str(e)}

def covariance_release_psd_steps(steps):
 from .accel_geometry_cell import accel_h_from_force_rotation,rotation_cell_from_ball
 import numpy as np
 F,Q,_,Racc,*_=one_sample_boxes();P=covariance_interval()
 Rm,Rr=rotation_cell_from_ball(math.radians(6.9));fm=np.zeros(3);fr=np.full(3,9.80665+C["marine_motion"]["A_max_mps2"])
 H=accel_h_from_force_rotation(fm,fr,Rm,Rr)
 checkpoints=[];targets={1,10,100,1000,10000,50000,100000,math.ceil(469./C["sensor_model"]["sample_period_min_s"])}
 for k in range(1,steps+1):
  P=predict_covariance(P,F,Q)
  P,cert=linked_covariance_update(P,H,Racc)
  P=aw_covariance_floor_event(P,16.)
  if k in targets:
   lo,hi=spectral_box(P);checkpoints.append({"step":k,"spectral_lower":lo,"spectral_upper":hi,
      "inverse_certificate":cert.get("certificate"),"finite":math.isfinite(hi)})
   if not math.isfinite(hi):break
 return {"verified":all(x["finite"] for x in checkpoints),"steps_requested":steps,
         "steps_completed":checkpoints[-1]["step"] if checkpoints else steps,"checkpoints":checkpoints}

def staged_psd_probe():
 n=math.ceil(469./C["sensor_model"]["sample_period_min_s"])
 return covariance_release_psd_steps(n)
