"""Signed field-axis frame-rotation budget under candidate gyro qualification."""
from fractions import Fraction as F
T=F(17); THETA=F("0.10471975511965977") # retained 6-deg captured domain
BGS=F("0.02"); DGS=F("0.00001"); CG=F("0.002")

def certificate():
    # Lemma I*: injection-product angle <= endpoint attitude errors +
    # integrated gyro residual. Candidate FAST contributes <=Cg on T<60.
    # Slow residual integral is conservatively min(B_s*T, 2 B_s/D_s style
    # amplitude); no zero-mean assumption. This bound is intentionally exposed:
    # it is likely too loose unless the slow component is paired by Abel rate.
    raw_slow=BGS*T
    endpoint=2*THETA
    raw=endpoint+raw_slow+CG
    # Rate-paired slow charge already available in the velocity telescope.
    return {"qualification":"OU3_FIELD_AXIS_ROTATION_BUDGET_V1",
      "endpoint_attitude_angle_rad":str(endpoint),
      "raw_slow_integral_rad":str(raw_slow),
      "fast_integral_rad":str(CG),
      "raw_signed_injection_angle_ceiling_rad":str(raw),
      "raw_angle_route_small":raw<F("0.2"),
      "slow_must_use_rate_paired_Abel":True,
      "per_reset_norm_sum_used":False,
      "source_uniform_twist_sum_closed":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
