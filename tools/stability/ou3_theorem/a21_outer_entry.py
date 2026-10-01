"""Outer-to-inner A21 entry reduction with temporal SLOW+FAST inputs."""
from __future__ import annotations
def entry_reduction():
 return {
  "captured_domain_release_set_compact":True,
  "zero_translation_moving_alias_removed_by_displacement_excitation":True,
  "inner_radius":"0.15",
  "inner_zero_action_invariant_set_empty":True,
  "compact_annulus_positive_homogeneous_loss_exists":True,
  "temporal_fast_profiles_numerically_qualified":False,
  "marine_excitation_constants_numerically_qualified":False,
  "linked_finite_supply_strictly_below_annulus_loss":False,
  "outer_retention_closed":False,"finite_inner_entry_closed":False,
  "reason":"strict supply-vs-loss margin needs qualified H_a,C_a,H_g,C_g,T_E,theta_E,T_P,P_E on the same history; no invented values"}
def conditional_entry(*,loss_floor,supply_ceiling):
 if loss_floor<=0 or supply_ceiling<0: raise ValueError("positive loss and nonnegative supply required")
 return {"strict_margin":loss_floor-supply_ceiling,
         "finite_entry":supply_ceiling<loss_floor}
