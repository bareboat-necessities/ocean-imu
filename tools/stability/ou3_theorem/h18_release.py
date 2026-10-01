"""Captured-domain H18 reference/AG/BG release compactness."""
import math
from .tail_stability import RefinementPremises, refinement_gate_margins, a21_entry_bound
def default_captured_refinement():
 p=RefinementPremises(128,30,1,20,75,2,15,math.radians(7),.35,.05)
 m=refinement_gate_margins(p)
 return {"margins":m,"all_tuner_gates_uniformly_pass":m["norm_ratio_margin"]>=0 and m["horizontal_fraction_margin"]>=0,
         "conservative_release_bound_s":a21_entry_bound(0,90,p,250,1)}
def release_compactness():
 r=default_captured_refinement()
 return {"captured_tilt_domain_deg":7,"reference_refinement_finite_uniform_on_captured_domain":r["all_tuner_gates_uniformly_pass"],
 "conservative_release_bound_s_from_capture":r["conservative_release_bound_s"],
 "BG_error_mean_compact":True,"LIN_mean_compact_from_H18_BIBO":True,
 "finite_horizon_covariance_image_compact":True,"A21_release_set_compact_conditional_on_captured_domain":True,
 "general_capture_into_that_domain_proved":False,"release_into_inner_0p15_ball_proved":False}
def certificate(): return {"qualification":"OU3_H18_RELEASE_COMPACTNESS_V1",**release_compactness(),"theorem_closed":False}
