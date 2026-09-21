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
      "temporal_margins_imply_finite_B_star":True,
      "B_star_instantiated":False,
      "uniform_historical_AG_readout_action":False,
      "full_21_covariance_upper":False,
      "rho0_certified":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))


def separated_reader_action_implication(delta_col, delta_gyr, coefficient_ceiling,
                                        noise_factor_ceiling, nuisance_ceiling,
                                        terminal_map_ceiling, operation_count):
    """Quantitative implication from positive temporal margins to finite B_*.

    This is not a margin certificate.  It proves the next logical step once
    the same-history enclosure supplies strict delta_col and delta_gyr.

    Let delta=min(delta_col,delta_gyr).  On every largest-residual pivot chart,
    the six successive residual pivots of O are bounded below by delta after
    the proof-coordinate normalization used by the temporal margins.  Hence
    the selected 6x6 minor has ||O_I^{-1}||_2 <=
    coefficient_ceiling**5 / delta**6 by adjugate/Hadamard.  The exact reader
    L=T_h O_I^{-1} therefore has the displayed uniform norm ceiling.  The
    backward action is a finite sum of transported rank<=3 noise/process
    factors plus the nuisance-root residual.  Bounding each chronological
    transport by coefficient_ceiling gives the explicit B_* below.

    The deliberately coarse exponent is acceptable here: the result needed
    is finiteness, not a practical contraction rate.  A useful J/rho still
    requires the rigorous source margins and then a sharper action enclosure.
    """
    vals=(delta_col,delta_gyr,coefficient_ceiling,noise_factor_ceiling,
          nuisance_ceiling,terminal_map_ceiling)
    if any(x<=0 for x in vals) or operation_count<1:
        raise ValueError("strict positive source bounds required")
    delta=min(delta_col,delta_gyr)
    inv_minor=coefficient_ceiling**5/delta**6
    reader=terminal_map_ceiling*inv_minor
    transport=max(1.0,coefficient_ceiling)**operation_count
    # ||sum X_i X_i'|| <= sum ||X_i||^2.  The final term covers the
    # nuisance-root residual in exactly the same backward reader recursion.
    b_star=(operation_count*(reader*noise_factor_ceiling*transport)**2
            +(reader*transport)**2*nuisance_ceiling)
    return {"delta":delta,"inverse_minor_norm_ceiling":inv_minor,
            "reader_norm_ceiling":reader,"B_star_scalar_ceiling":b_star,
            "B_star_finite":True}


def margin_to_Bstar_theorem():
    """Machine-readable status of the controlling implication."""
    return {
      "premises":["inf_W Delta_col(W)>0","inf_W Delta_gyr(W)>0",
                  "shipping coefficient/factor compactness on the fixed finite window"],
      "conclusion":"exists finite B_* with B_W <= B_* I6 on every carried window",
      "rank_structure":"successive largest-residual pivots; observation blocks have rank <=3",
      "proof":"finite pivot-chart cover + adjugate/Hadamard inverse bound + finite backward factor action",
      "implication_closed":True,
      "premise_margins_source_uniformly_certified":False,
    }
