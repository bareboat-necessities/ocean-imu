"""Held-H18 LIN covariance compactness and conditional BIBO certificate."""
from fractions import Fraction as F
from .nuisance_upper_certificate import bounds as nuisance_bounds

def four_s_zero_action_kernel(times):
    t=tuple(F(x) for x in times)
    if len(t)!=4 or any(b<=a for a,b in zip(t,t[1:])): raise ValueError("four increasing S times")
    a,b,c=t[:3]; det=(b-a)*(c-a)*(c-b)/2
    return {"neutral_three_row_det":det,"full_12D_kernel_zero":det!=0}

def h18_lin_covariance_upper():
    aw,fresh,w,sd,diag=nuisance_bounds()
    return {"coordinate_order":["v","p","S","a_w"],"upper_diagonal":[str(x) for x in diag[:4]],
      "AW_unconditioned_variance_ceiling":str(aw),"fresh_neutral_noise_coefficient_squared":str(fresh),
      "window_s":"17","applied_S_gap_ceiling_s":"0.156",
      "held_BA_activity_used_in_first_four_blocks":False,
      "accelerometer_and_magnetic_corrections_omitted_in_dominating_comparison":True,
      "P_LL_uniform_upper_after_17s_held_H18":True,"cross_covariance_inside_LL_retained":True}

def coefficient_compactness():
    return {"dt_s":["0.004","0.006"],"tau_s":["0.02","12"],"sigma_aw_upper_mps2":"4",
      "R_S_sigma_upper":"100","pseudo_S_gap_upper_s":"0.156","measurement_noise_positive":True,
      "gain_continuity_on_bounded_covariance_and_positive_innovation_noise":True,
      "held_H18_LIN_coefficient_family_compact_after_17s":True}

def certificate():
    return {"qualification":"OU3_HELD_BA_LIN_BIBO_V2","homogeneous_metric_nonexpansion":True,
      "four_S_structural_detectability":True,"held_H18_LIN_covariance_upper":h18_lin_covariance_upper(),
      "held_H18_coefficient_compactness":coefficient_compactness(),
      "source_uniform_strict_factor_exists_by_compactness":True,"numeric_rho_LIN":None,
      "uniform_homogeneous_rho_LIN_lt_1_exists":True,"bounded_affine_input_class_proved":False,
      "all_time_affine_BIBO_closed":False,"qTv_boundary_action_closed":False,
      "release_LIN_mean_compactness_closed":False,
      "reason_remaining":"need source-uniform affine AG/BG/held-BA/physical SLOW+FAST input bound on the fixed strict word",
      "theorem_closed":False}
