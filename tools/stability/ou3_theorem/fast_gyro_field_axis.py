"""Candidate FAST-gyro charge in the signed field-axis telescope."""
from fractions import Fraction as F
import math
T=F(17); P=F("5.5"); Bs=F("0.02"); Ds=F("0.00001")
Hg=F(60); Cg=F("0.002"); Bf=F("0.02")

def fast_integral(T):
    n=T//Hg; r=T-n*Hg
    return min(Bf*T,n*Cg+min(Bf*r,Cg))

def certificate():
    # For a coefficient normalized by T, the FAST signed primitive contributes
    # at most P*K*(T)/T before coefficient-rotation variation. Slow term is the
    # existing 2 P B_s/T + P D_s bound.
    slow=2*P*Bs/T+P*Ds
    fast=P*fast_integral(T)/T
    return {
      "qualification":"OU3_FAST_GYRO_FIELD_AXIS_CANDIDATE_V1",
      "window_s":str(T),"slow_charge_mps2":str(slow),
      "fast_signed_charge_mps2":str(fast),
      "slow_plus_fast_charge_mps2":str(slow+fast),
      "old_fast_amplitude_charge_mps2":str(P*Bf),
      "old_0p055_sinusoid_requires_halfcycle_integral_rad":"0.04",
      "candidate_Cg_rad":str(Cg),
      "old_sinusoid_excluded_by_candidate_fast_cap":Cg<F("0.04"),
      "coefficient_rotation_variation_included":False,
      "literal_jump_cancellation_included":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
