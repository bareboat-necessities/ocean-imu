"""Non-promoting comparison algebra for a quotient/covariance tube.

Scalar inputs cannot certify uniform derivatives, branch coverage, inherited
entry, gauge amplitude, every-prefix retention or all future service windows.
This solver therefore NEVER changes an admission/theorem flag to true.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

def solve(rho_p,rho_x,b,c,q_p,q_x,service_budget,information_lipschitz=None,
          *,gauge_injection=None,gauge_amplitude=None,covariance_gauge_injection=None):
    values=(rho_p,rho_x,b,c,q_p,q_x,service_budget,information_lipschitz,
            gauge_injection,gauge_amplitude,covariance_gauge_injection)
    if any(v is not None and (not math.isfinite(v) or v<0) for v in values):
        raise ValueError("finite nonnegative comparison data required")
    out={"qualification":"OU3_PLANAR_COMMON_NORM_RADIUS_V2",
      "result_type":"CONDITIONAL", "scope":"comparison algebra only; input certification is external and OPEN",
      "rho_P":rho_p,"rho_quotient":rho_x,"service_weyl_budget":service_budget,
      "b_quotient_to_P":b,"c_P_to_quotient":c,"P_additive_charge":q_p,
      "transverse_additive_charge_before_gauge":q_x,
      "gauge_injection_norm":gauge_injection,"gauge_amplitude_bound":gauge_amplitude,
      "covariance_gauge_injection_norm":covariance_gauge_injection,
      "information_cell_lipschitz":information_lipschitz,
      "comparison_fixed_point_exists":False,"joint_cell_forward_invariant":False,
      "information_cell_perturbation_certified":False,
      "all_time_magnetic_service_verified":False,"theorem_closed":False}
    if None in (rho_p,rho_x,b,c,q_p,q_x,gauge_injection,gauge_amplitude,covariance_gauge_injection):
        out["missing"]="uniform quotient/covariance gains and charges, C_Q, C_P and physical gauge amplitude"
        return out
    product_budget=(1-rho_p)*(1-rho_x)
    out["cross_product_budget"]=product_budget
    if rho_p>=1 or rho_x>=1 or b*c>=product_budget:
        out["failure_class"]="D_SUFFICIENT_BOUND_FAILURE"
        out["missing"]="valid contracting comparison; failure is not a shipping counterexample"
        return out
    charge=q_x+gauge_injection*gauge_amplitude
    charge_p=q_p+covariance_gauge_injection*gauge_amplitude
    det=product_budget-b*c
    rp=((1-rho_x)*charge_p+b*charge)/det
    rq=(c*charge_p+(1-rho_p)*charge)/det
    if not all(math.isfinite(v) for v in (charge,charge_p,det,rp,rq)):
        raise ValueError("comparison arithmetic overflow")
    out.update({"candidate_P_radius":rp,"candidate_quotient_radius":rq,
                "quotient_additive_charge_including_gauge":charge,
                "P_additive_charge_including_gauge":charge_p,
                "comparison_determinant":det,"comparison_fixed_point_exists":True})
    if information_lipschitz is not None:
        di=information_lipschitz*max(rp,rq)
        if not math.isfinite(di): raise ValueError("information arithmetic overflow")
        out["candidate_information_perturbation"]=di
        out["candidate_information_margin_positive"]=di<service_budget
    out["missing"]="certified same-history cell, inherited entry, all branches/prefixes/phases and every future physical service window"
    return out

def main():
    p=argparse.ArgumentParser(); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    # No private Mahony or point homogeneous factor may supply rho_quotient.
    o=solve(None,None,None,None,None,None,6.024764605642485)
    a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n"); print(json.dumps(o,sort_keys=True))
if __name__=="__main__": main()
