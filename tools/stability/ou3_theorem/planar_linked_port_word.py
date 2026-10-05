"""Literal 20-s linked-port tangent calculation scaffold.

This consumes the actual operation stream.  It computes exact point tangent
blocks for the covariance-only port and records which mean-dependent row
derivatives still require source-level tangent export.  It never promotes
point derivatives to uniform bounds.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .planar_service_frechet import frechet
from .planar_service_audit import audit

def calculate(stream:Path):
 f=frechet(stream); a=audit(stream)
 rho=max(v["relative_Frobenius_gain"] for v in f["parity_blocks"].values())
 return {"qualification":"OU3_PLANAR_LINKED_PORT_WORD_V1",
   "result_type":"FINITE POINT TANGENT + OPEN UNIFORM PORTS",
   "word_samples":[f["root_sample"],f["end_sample"]],"rho_P_point":rho,
   "b_mean_to_P_uniform":None,"c_P_to_mean_uniform":None,
   "q_P_uniform":None,"q_mean_uniform":None,
   "required_source_tangent_exports":[
     "dH_acc/d(mean) at each literal accepted accelerometer correction",
     "dH_mag/d(mean) at each literal accepted magnetic correction",
     "dK/dP linked through the same S,H,P at each correction",
     "same-boundary affine mean charge under fixed measurement history"],
   "local_operation_residual_max":max(a["maximum_absolute_operation_residuals"].values()),
   "classification":"OPEN; finite point tangent is not a uniform cell bound",
   "joint_cell_forward_invariant":False,"all_time_magnetic_service_verified":False,
   "theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=calculate(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
