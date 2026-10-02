"""Explicit shipping goLive seed for constructive release propagation."""
from __future__ import annotations
import math
from .interval_riccati_21 import diagonal_interval

def certificate():
 # initialize_from_attitude rebuilds attitude covariance. A later tilt watchdog
 # can use accel-only yaw=1.5708, so use that literal worst seed.
 p_att=1.5708**2
 diag=[p_att]*3+[1e-6]*3+[1.0]*3+[400.0]*3+[2500.0]*3+[16.48]*3+[0.004**2]*3
 return {"qualification":"OU3_GOLIVE_CONSTRUCTIVE_SEED_V1",
  "mean_state":[0.0]*21,
  "covariance_diagonal_upper":diag,
  "cross_covariance_zero_at_handoff":True,
  "attitude_seed_variance_upper":p_att,
  "BG_variance":1e-6,"v_variance":1.0,"p_variance":400.0,"S_variance":2500.0,
  "AW_variance_upper":16.48,"BA_variance":0.004**2,
  "frontend_drove_MEKF_before_handoff":False,
  "constructive_seed_available":True}

def covariance_interval():
 d=certificate()["covariance_diagonal_upper"]
 return diagonal_interval(tuple(d),tuple(d))

def covariance_factor_interval():
 from .psd_factor_innovation import point_cholesky_from_diagonal
 return point_cholesky_from_diagonal(certificate()["covariance_diagonal_upper"])
