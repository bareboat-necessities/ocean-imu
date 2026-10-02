"""Explicit release-set certificate assembled only from proved components.

Fail closed on the currently missing AG covariance numeric ceiling; this module
prevents qualitative compactness from being mistaken for a usable P0 box.
"""
from __future__ import annotations
from fractions import Fraction
from .h18_release import release_compactness
from .nuisance_upper_certificate import certificate as nuisance_certificate

def certificate():
 h=release_compactness();n=nuisance_certificate()
 return {"qualification":"OU3_EXPLICIT_A21_RELEASE_SET_V1",
   "captured_tilt_domain_deg":h["captured_tilt_domain_deg"],
   "LIN_BA_covariance_upper_diagonal":n["upper_diagonal"],
   "LIN_mean_compact":h["LIN_mean_compact_from_H18_BIBO"],
   "BG_mean_compact":h["BG_error_mean_compact"],
   "finite_covariance_image_compact":h["finite_horizon_covariance_image_compact"],
   "AG_covariance_numeric_upper_available":False,
   "BG_mean_numeric_radius_available":False,
   "LIN_mean_numeric_radius_available":False,
   "BA_graph_interval_available":False,
   "explicit_full_P0_interval_available":False,
   "constructive_leaf_seed_available":True,
   "constructive_seed_location":"goLive before finite H18/refinement propagation",
   "reason":"release is constructed as finite literal image of the explicit goLive seed; standalone precomputed LIN radius is unnecessary",
   "qualitative_compactness_not_promoted_to_numeric_box":True,
   "release_image_must_be_computed_by_causal_propagator":True,
   "theorem_closed":False}
