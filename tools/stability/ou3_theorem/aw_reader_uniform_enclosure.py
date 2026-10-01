"""Dependency-preserving source-uniform AW-reader enclosure contract.

This is proof infrastructure, not a promoted numerical certificate. It encodes
the exact decomposition that a rigorous interval/factor backend must enclose
without independent coefficient boxes.
"""
from fractions import Fraction as F
from .nonrecurring_accel_bridge import certificate as physical_bridge

ALLOW=F(physical_bridge()["required_full_state_reader_budget_mps2"])

def compose_enclosure(*, carried_center, factor_radius, nonlinear_radius,
                      arithmetic_radius):
    """One-ball enclosure after exact chronological composition.

    Inputs are radii for the COMPOSED functional Phi, not independent gains or
    source terms. Hence triangle addition occurs only between rigorously
    separated remainder classes after process/acc/S pairing and endpoint
    telescoping have already happened.
    """
    vals=tuple(F(x) for x in (carried_center,factor_radius,nonlinear_radius,
                              arithmetic_radius))
    if any(x<0 for x in vals): raise ValueError("nonnegative composed radii")
    c,f,n,a=vals; upper=c+f+n+a
    return {"carried_center_mps2":str(c),"factor_enclosure_radius_mps2":str(f),
            "nonlinear_radius_mps2":str(n),"arithmetic_radius_mps2":str(a),
            "source_uniform_upper_mps2":str(upper),
            "allowance_mps2":str(ALLOW),"strict_margin_mps2":str(ALLOW-upper),
            "closed":upper<ALLOW}

def required_factor_radius(carried_center=F("0.371143"),
                           nonlinear_radius=F(0),arithmetic_radius=F(0)):
    c,n,a=map(F,(carried_center,nonlinear_radius,arithmetic_radius))
    return ALLOW-c-n-a

def certificate():
    return {
      "qualification":"OU3_AW_READER_UNIFORM_ENCLOSURE_CONTRACT_V1",
      "functional":"Phi(h)=P_B L_AW(h)b_full(h)",
      "reachable_window_s":"16",
      "allowance_mps2":str(ALLOW),
      "carried_stress_center_mps2":"0.371143",
      "available_total_enclosure_radius_mps2":str(required_factor_radius()),
      "independent_gain_boxes_allowed":False,
      "independent_tuner_extrema_allowed":False,
      "independent_source_norms_allowed":False,
      "process_acc_S_pair_before_wrap":True,
      "physical_endpoint_telescope_before_wrap":True,
      "shared_history_dependency_required":True,
      "source_uniform_factor_radius_certified":False,
      "nonlinear_radius_certified":False,
      "float32_radius_certified":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
