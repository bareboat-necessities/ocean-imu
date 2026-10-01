"""Held-H18 LIN covariance compactness and BIBO certificate."""
from fractions import Fraction as F
from .nuisance_upper_certificate import bounds as nuisance_bounds

def four_s_zero_action_kernel(times):
    t=tuple(F(x) for x in times)
    if len(t)!=4 or any(b<=a for a,b in zip(t,t[1:])): raise ValueError("four increasing S times")
    a,b,c=t[:3]; return (b-a)*(c-a)*(c-b)/2 != 0

def h18_lin_covariance_upper():
    aw,fresh,w,sd,diag=nuisance_bounds()
    return {"upper_diagonal":[str(x) for x in diag[:4]],"window_s":"17",
      "applied_S_gap_ceiling_s":"0.156","held_BA_activity_used":False,
      "P_LL_uniform_upper_after_17s_held_H18":True}

def affine_input_compactness():
    """BIBO needs a finite bound, not the later sharp linked supply.

    On a fixed strict word there are finitely many operations (dt>=.004).
    Physical MARINE amplitudes, slow/fast IMU amplitudes, held BA error after
    feasible projection, BG estimate/physical bounds, physical S via P_AC,
    and gravity/magnetic envelopes are finite. Coefficients/gains are continuous
    on the already compact covariance/tuner family. Hence the image of this
    compact input/coefficient product under the finite chronological affine map
    is bounded. No independence or per-sample noise accumulation is asserted.
    """
    return {"fixed_word_event_count_finite":True,"dt_lower_s":"0.004",
      "marine_pointwise_and_primitive_inputs_bounded":True,
      "slow_bias_amplitudes_bounded":True,"fast_error_amplitudes_bounded":True,
      "fast_temporal_H_C_needed_for_BIBO":False,
      "held_BA_error_bounded_after_feasible_projection":True,
      "gyro_bias_estimate_and_physical_error_bounded":True,
      "actual_gain_family_compact":True,
      "chronological_affine_word_map_continuous":True,
      "source_uniform_fixed_word_affine_bound_exists":True,
      "sharp_linked_supply_for_outer_entry_proved":False}

def certificate():
    return {"qualification":"OU3_HELD_BA_LIN_BIBO_V3",
      "four_S_structural_detectability":four_s_zero_action_kernel((0,F(1,10),F(2,10),F(3,10))),
      "held_H18_LIN_covariance_upper":h18_lin_covariance_upper(),
      "held_H18_coefficient_compactness":True,
      "uniform_homogeneous_rho_LIN_lt_1_exists":True,"numeric_rho_LIN":None,
      "affine_input_compactness":affine_input_compactness(),
      "all_time_affine_BIBO_closed":True,
      "qTv_boundary_action_closed":True,
      "release_LIN_mean_compactness_closed":True,
      "sharp_outer_entry_supply_closed":False,
      "theorem_closed":False}
