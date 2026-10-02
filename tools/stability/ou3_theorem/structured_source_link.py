# ruff: noqa: F401, F811
"""Same-history source identities for structured AW release propagation.

These are identities/coordinate maps, not independent source boxes.  A promoted
leaf must construct every argument from one persistent physical/filter history.
"""
from __future__ import annotations
import math
import numpy as np

def acc_residual_components():
 return {
  "identity":"r_acc = f_meas - (R_wb*(a_w-g) + lever + b_a)",
  "error_identity":"r_acc = H_acc e + d_acc on the same pre-correction boundary",
  "sources":["physical_specific_force","slow_accel_bias","fast_accel_error",
             "measurement_model_defect","lever_arm","BA_mean","latent_aw","attitude"],
  "independent_suprema_forbidden":True}

def s_pseudo_residual_components():
 return {"identity":"r_S = -S_hat = e_S - S_phys",
         "sources":["physical_S","integrated_OU_chain","previous_acc_corrections",
                    "previous_S_corrections"],
         "independent_suprema_forbidden":True}

def exact_acc_residual(f_meas,R_wb,aw,g_world,lever,ba):
 """Literal shipping accelerometer innovation from one history boundary."""
 return np.asarray(f_meas,float)-(
  np.asarray(R_wb,float)@(np.asarray(aw,float)-np.asarray(g_world,float))
  +np.asarray(lever,float)+np.asarray(ba,float))

def residual_local_defect(r,H,e):
 """Exact additive measurement defect d=r-H e at the same boundary."""
 return np.asarray(r,float)-np.asarray(H,float)@np.asarray(e,float)

def exact_s_residual(S_hat):
 return -np.asarray(S_hat,float)

def physical_s_source(e_S,S_phys):
 """Identity r_S=e_S-S_phys when e_S=S_phys-S_hat."""
 return np.asarray(e_S,float)-np.asarray(S_phys,float)

def linked_norm(v):
 return float(np.linalg.norm(np.asarray(v,float)))

def required_leaf_exports():
 return ["acc_residual_vector","acc_local_defect_vector","S_residual_vector",
         "physical_S_vector","Kaw","KawS","KawMag","INJ_b","INJ_c",
         "aw_mean_vector","phi"]

def audit():
 return {"source_uniform_verified":False,
         "same_history_required":True,
         "independent_suprema_forbidden":True,
         "acc_link":"r_acc=H_acc e+d_acc",
         "S_link":"r_S=e_S-S_phys",
         "note":"Existing local-defect bounds constrain d_acc; they do not by themselves bound r_acc without the linked state e."}
