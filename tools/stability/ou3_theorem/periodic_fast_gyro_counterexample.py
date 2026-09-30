"""Constructive periodic fast-gyro witness certificate.

Non-promoting until every interval check is true. The forcing period is 6 s:
1200 literal 5 ms IMU samples and 150 deployed 25 Hz magnetic callbacks.
"""
from fractions import Fraction as F
DT=F(1,200); PERIOD_S=F(6); SAMPLES=1200; MAG_STRIDE=8
VAMP=F(7,2); NGAMP=F(1,50)

def certificate():
    pi_lo,pi_hi=F(333,106),F(355,113)
    om_lo,om_hi=pi_lo/3,pi_hi/3
    p_hi=VAMP/om_lo
    a_hi=VAMP*om_hi
    j_hi=VAMP*om_hi*om_hi
    theta_hi=NGAMP/om_lo
    supply=VAMP*NGAMP/2
    return {
      "qualification":"OU3_PERIODIC_FAST_GYRO_WITNESS_V1",
      "period_s":str(PERIOD_S),"samples_per_period":SAMPLES,
      "magnetic_callbacks_per_period":SAMPLES//MAG_STRIDE,
      "omega_rad_s_lower":str(om_lo),"omega_rad_s_upper":str(om_hi),
      "p_amp_upper_m":str(p_hi),"v_amp_mps":str(VAMP),
      "a_amp_upper_mps2":str(a_hi),"jerk_amp_upper_mps3":str(j_hi),
      "field_axis_error_amp_upper_rad":str(theta_hi),
      "signed_fast_gyro_supply_mps2":str(supply),
      "old_margin_mps2":"0.0338084",
      "supply_exceeds_old_margin":supply>F("0.0338084"),
      "physical_envelopes_pre_gravity_compensation":
          p_hi<F("8.1") and VAMP<=F("5.5") and a_hi<F("8.8") and j_hi<F(100),
      "native_periodic_orbit_exported":False,
      "periodic_tuner_word_interval_enclosed":False,
      "periodic_covariance_orbit_interval_enclosed":False,
      "physical_compatibility_operator_interval_enclosed":False,
      "physical_compatibility_operator_nonsingular":False,
      "physical_solution_interval_enclosed":False,
      "all_shipping_gates_verified":False,
      "magnetic_service_verified":False,
      "weighted_signed_functional_interval_lower":None,
      "shipping_counterexample_certified":False,
      "theorem_closed":False}
if __name__=="__main__":
    import json; print(json.dumps(certificate(),indent=2,sort_keys=True))
