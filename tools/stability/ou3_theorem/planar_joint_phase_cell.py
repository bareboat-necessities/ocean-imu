"""Fail-closed joint phase-cell feasibility thresholds.

This calculation does not certify the joint cell.  It converts the measured
fixed-coefficient covariance contraction into a reported point quantity.
The private Mahony decrement is reported separately and supplies no MEKF budget, and converts the carried
physical 2x2 magnetic floor into the exact Weyl perturbation budget that a
source-uniform joint cell must beat.

No independent marginal interval is promoted to reachability.
"""
from __future__ import annotations
import argparse, json, math
from fractions import Fraction
from pathlib import Path

def calculate(frechet: dict, mahony: dict, service_floor: float) -> dict:
    if frechet.get("result_type") != "FINITE DIAGNOSTIC ONLY":
        raise ValueError("covariance input must remain finite diagnostic")
    if mahony.get("all_time_frontend_tube_verified") is not False:
        raise ValueError("frontend tube must remain unpromoted")
    gains=[float(frechet["parity_blocks"][p]["relative_Frobenius_gain"])
           for p in ("even","odd")]
    a=max(gains)
    delta=Fraction(mahony["squared_norm_decrement_lower"])
    h=Fraction(mahony["h_binary32_exact"])
    steps=round(20.0/float(h))
    # Homogeneous common-quadratic comparison only.  The certified additive
    # charge is NOT discarded; it belongs in the invariant-radius inequality.
    d=math.exp(0.5*steps*math.log1p(-float(delta)))
    # For nonnegative comparison matrix [[a,b],[c,d]], rho<1 iff
    # a<1,d<1 and b*c < (1-a)(1-d).
    # d is a PRIVATE MAHONY gain, not the MEKF mean or quotient gain.
    coupling_product_budget=None
    service_weyl_budget=float(service_floor)-1.0
    if not (0<a<1 and 0<d<1 and service_weyl_budget>0):
        raise ValueError("feasibility margin is not positive")
    return {
      "qualification":"OU3_PLANAR_JOINT_PHASE_CELL_THRESHOLD_V1",
      "result_type":"CONDITIONAL",
      "covariance_partial_gain_max":a,
      "mahony_homogeneous_20s_gain":d,
      "mahony_steps":steps,
      "joint_two_block_small_gain_condition":
        "b_Q*c_Q < (1-rho_P)*(1-rho_Q), with certified linked quotient blocks",
      "coupling_product_budget":coupling_product_budget,
      "symmetric_coupling_norm_budget":None,
      "rho_quotient_uniform":None,
      "mahony_gain_used_as_mekf_mean_gain":False,
      "gauge_forcing_required":"C_Q alpha, with amplitude from the physical SLOW+FAST constraints",
      "carried_physical_2x2_service_floor":float(service_floor),
      "weyl_information_perturbation_budget_to_muM_1":service_weyl_budget,
      "required_information_condition":
        "for each physical +/- 2x2 block, ||I_cell-I_carried||_2 < carried_floor-1",
      "additive_mahony_charge_retained_separately":True,
      "generated_reverse_feedback_ports_zero":True,
      "generated_one_way_coefficient_ranges_all_time_bounded":False,
      "remaining_bidirectional_loop":"MEKF mean <-> covariance through state-dependent H and P-dependent K",
      "generated_mean_coefficient_couplings_bounded":False,
      "joint_phase_cell_forward_invariant":False,
      "every_placed_window_magnetic_service_verified":False,
      "all_time_magnetic_service_verified":False,
      "theorem_closed":False,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frechet",type=Path,required=True)
    p.add_argument("--mahony",type=Path,required=True)
    p.add_argument("--service-floor",type=float,required=True)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    out=calculate(json.loads(a.frechet.read_text()),json.loads(a.mahony.read_text()),a.service_floor)
    text=json.dumps(out,indent=2,sort_keys=True)+"\n"
    if a.output:a.output.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
