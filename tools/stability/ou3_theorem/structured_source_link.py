"""Map existing same-history source symbols to literal acc residual and S mean.

This module defines identities, not independent bounds.  A source-uniform leaf
must supply all terms jointly.
"""
from __future__ import annotations

def acc_residual_components():
 return {
  "identity":"r = f_meas - (R_wb*(a_w-g) + lever + b_a(temp))",
  "sources":["physical_specific_force","slow_accel_bias","fast_accel_error",
             "measurement_model_defect","lever_arm","BA_mean","latent_aw","attitude"],
  "independent_suprema_forbidden":True}

def s_pseudo_residual_components():
 return {"identity":"r_S = -S_mean",
         "sources":["integrated_OU_chain","previous_acc_corrections","previous_S_corrections"],
         "independent_suprema_forbidden":True}

def required_leaf_exports():
 return ["acc_residual_norm","S_mean_norm","Kaw_norm","KawS_norm","INJ_b","INJ_c","aw_mean_norm","phi"]
