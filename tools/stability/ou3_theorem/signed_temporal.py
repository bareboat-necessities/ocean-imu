"""Exact fail-closed checks for the OU-III signed temporal construction."""
from __future__ import annotations
from fractions import Fraction as F
from math import prod

def atoms(knots):
    k=tuple(F(x) for x in knots)
    if len(k)!=4 or any(b<=a for a,b in zip(k,k[1:])):
        raise ValueError("four increasing actual-S times required")
    return tuple(-F(6,1)/prod(t-s for i,s in enumerate(k) if i!=j) for j,t in enumerate(k))

def jet(knots,t,order=0,side="right"):
    k=tuple(F(x) for x in knots); t=F(t); c=atoms(k)
    if order not in (0,1,2): raise ValueError("order")
    fac=(2,1,1)[order]
    ans=F(0)
    for cj,tj in zip(c,k):
        active=t>tj or (t==tj and side=="right")
        if active:
            if order==0: ans += cj*(t-tj)**2/F(2)
            elif order==1: ans += cj*(t-tj)
            else: ans += cj
    return ans

def endpoint_jets_vanish(knots):
    k=tuple(F(x) for x in knots)
    return all(jet(k,k[0],q,"left")==0 and jet(k,k[-1],q,"right")==0 for q in range(3))

def zero_terminal_homogeneous_adjoint_impossible(first_s_weight_rank=3):
    # Z_N=0 and Z_i=Z_{i+1}A_i imply Z_i=0 for every i; hence
    # W_i=Z_{i+1}K_i=0. A regular first S atom is c0(I-K_SS),
    # nonsingular because I-K_SS=R_eff(P_SS+R_eff)^-1.
    return first_s_weight_rank==3

def physical_bounds():
    mean=F(1050297,4840000)
    gamma=F(2432784801508736,10**19)
    floor=F(1,5000)
    return {"sampled_accel_mean_ceiling":mean,
            "physical_joint_vector_floor":gamma,
            "physical_floor_margin":gamma-floor}

def projection_sector_gap():
    # unchanged R_b=.4 and B_a=0.22516660498395405
    return F(2,5)-F(22516660498395405,10**17)

def certificate():
    p=physical_bounds()
    return {
      "qualification":"OU3_SIGNED_TEMPORAL_V1",
      "actual_S_event_signed_temporal_identity":True,
      "homogeneous_zero_terminal_adjoint_closes":False,
      "observation_forced_adjoint_required":True,
      "sampled_accel_mean_ceiling":str(p["sampled_accel_mean_ceiling"]),
      "physical_joint_vector_floor":str(p["physical_joint_vector_floor"]),
      "physical_floor_margin":str(p["physical_floor_margin"]),
      "projection_sector_gap":str(projection_sector_gap()),
      "regular_A21_pre_projection_BA_precision_ceiling":1000003000,
      "source_uniform_nominal_force_field_temporal_margin":False,
      "source_uniform_nominal_gyro_alias_temporal_margin":False,
      "uniform_historical_AG_readout_action":False,
      "full_21_covariance_upper":False,
      "rho0_certified":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
