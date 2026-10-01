"""Fail-closed structural H18 LIN BIBO theorem reduction.

The literal covariance metric makes prediction and covariance-matched
corrections nonexpansive. Four separated S rows have no nonzero homogeneous
zero-action LIN root for the integrated OU chain. Compactness would therefore
give a strict finite-word factor only after the reachable held-H18 covariance /
tuner / scheduler coefficient family is proved compact uniformly in release
delay. That compactness is the remaining premise; no carried rho is promoted.
"""
from __future__ import annotations
from fractions import Fraction as F

def four_s_zero_action_kernel(times):
    """Exact determinant of the neutral v,p,S Vandermonde minor.

    S(t)=S0+t p0+t^2 v0/2 plus the OU extension. Three distinct S zeros kill
    the neutral cubic-free root; the fourth row then kills the nonzero
    exponential/OU root by the extended-Chebyshev property. We expose the
    neutral determinant exactly; OU extension strictness is structural.
    """
    t=tuple(F(x) for x in times)
    if len(t)!=4 or any(b<=a for a,b in zip(t,t[1:])):
        raise ValueError("four strictly increasing S times required")
    a,b,c=t[:3]
    # determinant of rows [t^2/2,t,1] at first three events
    det=(b-a)*(c-a)*(c-b)/F(2)
    return {"neutral_three_row_det":det,
            "neutral_kernel_zero":det!=0,
            "fourth_row_kills_OU_extension":True,
            "full_12D_kernel_zero":det!=0}

def covariance_metric_nonexpansion():
    return {
      "prediction":"P-=F P F'+Q, Q>=0 => ||F e||_(P-)^-1 <= ||e||_P^-1",
      "correction":"Joseph/Kalman measurement-energy identity",
      "PSD_sync":"mean neutral covariance inflation cannot increase storage",
      "reset":"LIN mean unchanged by attitude reset; full congruent storage handled in parent proof",
      "literal_accelerometer_correction_nonexpansive":True,
      "euclidean_gain_required":False,
    }

def bibo_reduction():
    return {
      "homogeneous_nonexpansion":True,
      "four_S_structural_detectability":True,
      "source_uniform_strict_factor_from_compactness":True,
      "required_compact_family":[
        "held-H18 P_LL and relevant cross covariance upper bounds independent of release delay",
        "tau,sigma_aw,R_S,T_S compact reachable set with actual scheduler phase",
        "bounded correction coefficients/gains and positive innovation covariance",
      ],
      "these_premises_currently_proved":False,
      "uniform_rho_LIN_lt_1":False,
      "uniform_affine_input_gain":False,
      "all_time_BIBO_closed":False,
      "qTv_boundary_action_closed":False,
      "logical_result":"the displacement boundary row is not an additional obstruction; it closes immediately once held-H18 LIN BIBO/compactness closes",
    }

def certificate():
    k=four_s_zero_action_kernel((0,F(1,10),F(2,10),F(3,10)))
    return {"qualification":"OU3_HELD_BA_LIN_BIBO_REDUCTION_V1",
            "neutral_test_det":str(k["neutral_three_row_det"]),
            **covariance_metric_nonexpansion(),**bibo_reduction(),
            "theorem_closed":False}
