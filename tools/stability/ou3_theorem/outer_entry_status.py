"""Existence-level outer A21 entry/absorption theorem status."""
from fractions import Fraction as F
LOCAL=F("0.15")**2

def certificate():
    return {
      "qualification":"OU3_OUTER_ENTRY_EXISTENCE_V1",
      "release_compactness_used":True,
      "recurring_root_covariance_lower_used":True,
      "compact_outer_release_storage_bound_exists":True,
      "finite_superword_homogeneous_strictness_exists":True,
      "linked_chi_supremum_finite_on_compact_outer_class":True,
      "finite_practical_absorbing_storage_radius_exists":True,
      "local_target_storage":str(LOCAL),
      "required_strict_entry_inequality":"sup_same_history chi_gamma/gamma < 0.0225",
      "required_prefix_condition":"linked prefix retention for every intermediate operation",
      "strict_entry_inequality_certified":False,
      "prefix_retention_certified":False,
      "entry_into_local_ball_proved":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
