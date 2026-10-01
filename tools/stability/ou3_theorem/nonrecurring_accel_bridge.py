"""Non-recurrent physical acceleration -> aggregate G0 bridge.

Uses the existing bounded-velocity/jerk physical-transfer theorem and the
candidate SLOW+FAST acceleration qualification. This module deliberately
separates what is proved from the still-open nominal AW-loop transfer.
"""
from fractions import Fraction as F

G=F("9.80665"); V=F("5.5"); J=F(100); H=F(6,1000)
L=F(16); BS=F("0.22516660498395405"); DS=F("0.001")
BF=F("0.3"); HF=F(60); CF=F(1,20)

def fast_accum(T):
    n=T//HF; r=T-n*HF
    return min(BF*T,n*CF+min(BF*r,CF))

def certificate():
    physical=2*V/L+J*H/4
    # For a normalized signed mean over a complete L-window, the fast mean
    # contribution is bounded by K*(L)/L. A single carried slow component has
    # mean amplitude <=BS; using only amplitude is conservative. This is a
    # sensor-record bound, NOT yet a nominal-AW bound because the literal
    # time-varying Kalman gain can rectify the input.
    fast_mean=fast_accum(L)/L
    sensor_mean=BS+fast_mean
    allowance=G/F(5)-physical
    return {
      "qualification":"OU3_NONRECURRING_ACCEL_BRIDGE_V1",
      "window_s":str(L),
      "physical_signed_mean_supply_mps2":str(physical),
      "corollary_Astar_total_threshold_mps2":str(G/F(5)),
      "signed_error_allowance_after_physical_mps2":str(allowance),
      "candidate_fast_mean_cap_mps2":str(fast_mean),
      "slow_amplitude_outer_mps2":str(BS),
      "sensor_error_mean_outer_mps2":str(sensor_mean),
      "sensor_record_margin_mps2":str(allowance-sensor_mean),
      "sensor_record_margin_positive":sensor_mean<allowance,
      "nominal_AW_mean_follows_sensor_mean":False,
      "reason_open":"literal time-varying acc/S gains can rectify; use signed AW-loop/two-Abel identity",
      "G0_nominal_window_premise_closed":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
