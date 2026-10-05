"""Fail-closed common-norm radius solver for the reduced planar joint cell.

The solver separates theorem inputs from finite diagnostics.  It may report a
candidate radius only after both cross-port gains and additive charges have
rigorous certificates.  Until then it reports the exact missing inequalities.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

def solve(rho_p:float,rho_x:float,b:float|None,c:float|None,q_p:float|None,q_x:float|None,
          service_budget:float,information_lipschitz:float|None=None):
    product_budget=(1-rho_p)*(1-rho_x)
    out={"qualification":"OU3_PLANAR_COMMON_NORM_RADIUS_V1",
      "result_type":"CONDITIONAL radius solver; fail closed on uncertified ports",
      "rho_P":rho_p,"rho_mean":rho_x,"cross_product_budget":product_budget,
      "service_weyl_budget":service_budget,"b_mean_to_P":b,"c_P_to_mean":c,
      "P_additive_charge":q_p,"mean_additive_charge":q_x,
      "information_cell_lipschitz":information_lipschitz,
      "joint_cell_forward_invariant":False,"information_cell_perturbation_certified":False,
      "all_time_magnetic_service_verified":False,"theorem_closed":False}
    if None in (b,c,q_p,q_x):
      out["missing"]="rigorous common-norm cross gains b,c and additive charges q_P,q_x"
      return out
    if b*c>=product_budget:
      out["failure_class"]="D_SUFFICIENT_BOUND_FAILURE"
      out["missing"]="sharper linked cross-port bound; not a shipping counterexample"
      return out
    # Minimal positive fixed point r=(I-A)^-1 q for comparison A.
    det=(1-rho_p)*(1-rho_x)-b*c
    rp=((1-rho_x)*q_p+b*q_x)/det
    rx=(c*q_p+(1-rho_p)*q_x)/det
    out.update({"candidate_P_radius":rp,"candidate_mean_radius":rx,
                "comparison_determinant":det,"joint_cell_forward_invariant":True})
    if information_lipschitz is None:
      out["missing"]="rigorous one-second information Lipschitz bound on the same cell"
      return out
    di=information_lipschitz*max(rp,rx)
    out["information_perturbation_upper"]=di
    out["information_cell_perturbation_certified"]=di<service_budget
    out["all_time_magnetic_service_verified"]=di<service_budget
    return out

def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=solve(.8738212970667967,.5827057639270218,None,None,None,None,6.024764605642485)
 a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
