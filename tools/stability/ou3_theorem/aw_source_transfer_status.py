"""Controlling source-only AW transfer status after exact telescopes."""
from fractions import Fraction as F
G=F("9.80665")
PHYS=F("0.8375")
SLOW_A=F("0.22516660498395405")
FAST_A=F("0.003125")
GYRO=F("0.0136432352941176")
THRESH=G/5

def certificate():
    charged=PHYS+SLOW_A+FAST_A+GYRO
    # Gyro is a conservative cross-channel charge retained here until the
    # final augmented field-axis normalization proves it belongs solely in
    # injection transport. It is intentionally not used to promote G0.
    return {
      "qualification":"OU3_AW_SOURCE_TRANSFER_STATUS_V1",
      "corollary_Astar_threshold_mps2":str(THRESH),
      "physical_mean_charge_mps2":str(PHYS),
      "slow_accel_charge_mps2":str(SLOW_A),
      "fast_accel_charge_mps2":str(FAST_A),
      "continuous_slow_fast_gyro_cross_charge_mps2":str(GYRO),
      "conservative_total_before_polynomial_float32_mps2":str(charged),
      "remaining_Astar_margin_before_float32_mps2":str(THRESH-charged),
      "real_arithmetic_reset_error_dependent_twist_in_source":False,
      "quaternion_polynomial_is_separate_AW_source":False,
      "quaternion_polynomial_accounting":"already included in qualified gyro prediction transport angle/singular-floor certificate; do not double count",
      "explicit_G0_point4_point12_premises_proved":False,
      "Astar_real_arithmetic_source_budget_closed":True,
      "Astar_source_uniform_strictness_closed":True,
      "reason_open":"float32 transfer and downstream injection/nonlinear retention remain separate obligations",
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
