"""Exact status of the signed correction-jump/twist remainder."""
from fractions import Fraction as F
from .corrected_word import shipping_reset_remainder_bound

R=F("0.20943951023931954") # pairwise error-difference <= 12 deg on captured 6-deg pair
D=F("0.007")               # conservative per-prediction/source angle cap, NOT correction injection cap
HEAD=F("0.52439539501604595")

def certificate():
    # This is intentionally only a scale diagnostic: D is not a proved bound
    # on every measurement injection. It shows why eventwise reset summation
    # cannot be promoted even after the linear signed jumps telescope.
    one=shipping_reset_remainder_bound(D,R,R)
    return {
      "qualification":"OU3_JUMP_TWIST_REMAINDER_STATUS_V1",
      "captured_pair_chart_radius_rad":str(R),
      "illustrative_small_injection_rad":str(D),
      "one_event_exact_reset_remainder_scale_rad":str(one),
      "available_AW_headroom_mps2":str(HEAD),
      "linear_jump_must_telescope":True,
      "eventwise_remainder_sum_allowed":False,
      "complete_causal_reader_action_ceiling":"16",
      "normalized_reader_norm_ceiling":"4",
      "higher_order_remainder_can_use_reader_action":True,
      "injection_squared_times_error_is_source_only":False,
      "injection_squared_times_error_destination":"linked nonlinear finite-error/prefix retention",
      "zero_error_twist_term":"literal quaternion polynomial/arithmetic defect only",
      "measurement_injection_uniform_bound_certified":False,
      "signed_twist_source_uniform_bound_certified":True,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
