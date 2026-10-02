"""Explicit numerical release bounds derivable from literal H18/refinement mechanics."""
from __future__ import annotations
from .h18_release import default_captured_refinement
from .nuisance_upper_certificate import certificate as nuisance
from .gyro_bias_projection import source_radius

def bounds():
 h=default_captured_refinement();T=float(h["conservative_release_bound_s"])
 ptheta0=1.5708**2;pb0=1e-6
 qtheta=.00135**2;qb=1e-11
 ptheta=ptheta0+T*T*pb0+qtheta*T+qb*T**3/3
 pbg=pb0+qb*T
 ag_spectral=3*(ptheta+pbg)
 bg_error=float(source_radius())+.02
 ba_graph_norm=9.80665+8.8
 n=nuisance()
 return {"release_horizon_s":T,"AG_covariance_spectral_upper":ag_spectral,
   "attitude_covariance_trace_component_upper":ptheta,"BG_covariance_component_upper":pbg,
   "BG_mean_error_norm_upper_rad_s":bg_error,
   "BA_graph_operator_norm_upper":ba_graph_norm,
   "BA_graph_entry_interval":[-ba_graph_norm,ba_graph_norm],
   "LIN_BA_covariance_upper_diagonal":n["upper_diagonal"],
   "LIN_mean_numeric_radius_available":False,
   "reason_LIN":"held-H18 BIBO theorem proves a finite radius but its rho_L and affine word forcing D_H are existential, not numerical",
   "shipping_counterexample":False}

def certificate():
 b=bounds();return {"qualification":"OU3_NUMERIC_RELEASE_PARTIAL_V1",**b,
  "AG_BG_BA_graph_numeric":True,"full_release_seed_numeric":False,"theorem_closed":False}
