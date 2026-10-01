"""Status of the corrected outer physical-to-nominal bridge."""
def candidate():
 reserve=.003056118
 return {"profile":{"T_s":60,"theta_deg":1,"P_E_m":.02,"C_a_mps":1.2,"C_g_rad":.0042},
 "old_fast_gyro_amplitude_witness_admissible":False,
 "reason_old_witness":"60-s carried primitive exceeds candidate C_g",
 "physical_reserve_rad_approx":reserve,
 "interior_base_innovations_eliminated_by_signed_adjoint":True,
 "interior_pointwise_AW_error_required":False,
 "E_state_reduced_to_external_endpoint_reader_plus_nonlinear_defects":True,
 "E_state_source_uniform_finite_from_compact_outer_class":True,
 "E_state_numeric_upper":None,
 "strict_outer_bridge_closed":False,
 "next_target":"bound augmented endpoint-reader operator on compact A21 release/outer class",
 "theorem_closed":False}
