"""Literal 20-s linked-port tangent calculation scaffold.

This consumes the actual operation stream.  It computes exact point tangent
blocks for the covariance-only port and records which mean-dependent row
derivatives still require source-level tangent export.  It never promotes
point derivatives to uniform bounds.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .planar_service_frechet import diagnostic as frechet
from .planar_service_audit import audit
from .planar_service_stream import records,expand
import numpy as np

def calculate(stream:Path):
 f=frechet(stream); a=audit(stream)
 rho=max(v["relative_Frobenius_gain"] for v in f["parity_blocks"].values())
 tangents=[]
 for kind,k,z in records(stream):
  if kind==9 and f["root_sample"] < k <= f["end_sample"]:
   P=expand(z[:225]);H=z[225:288].reshape(3,21);S=z[288:297].reshape(3,3)
   K=z[297:360].reshape(21,3);res=z[360:363];x=z[363:384]
   # Exact point norms of the operands entering the linked differential.
   tangents.append({"sample":k,"sensor_record_index":len(tangents),"P_norm":float(np.linalg.norm(P,2)),
     "H_norm":float(np.linalg.norm(H,2)),"Sinv_norm":float(np.linalg.norm(np.linalg.inv(S),2)),
     "K_norm":float(np.linalg.norm(K,2)),"residual_norm":float(np.linalg.norm(res)),
     "state_norm":float(np.linalg.norm(x))})
 if not tangents: raise ValueError("missing read-only mean tangent records")
 maxima={key:max(t[key] for t in tangents) for key in
         ("P_norm","H_norm","Sinv_norm","K_norm","residual_norm","state_norm")}
 return {"qualification":"OU3_PLANAR_LINKED_PORT_WORD_V1",
   "result_type":"FINITE POINT TANGENT + OPEN UNIFORM PORTS",
   "word_samples":[f["root_sample"],f["end_sample"]],"rho_P_point":rho,
   "point_tangent_records":len(tangents),"point_operand_norm_maxima":maxima,"same_operation_operands_exported":True,
   "b_mean_to_P_uniform":None,"c_P_to_mean_uniform":None,
   "q_P_uniform":None,"q_mean_uniform":None,
   "required_source_tangent_exports":[
     "uniformize exact dH_acc/d(mean) over the candidate cell",
     "uniformize exact dH_mag/d(mean) over the candidate cell",
     "propagate dK/dP linked through the same S,H,P at each correction",
     "same-boundary affine mean charge under fixed measurement history"],
   "tangent_block_initialization":{"D_xx":"I","D_Px":"0","D_xP":"0","D_PP":"I"},
   "endpoint_blocks_requested":["D_xx","D_Px","D_xP","D_PP"],
   "local_operation_residual_max":max(a["maximum_absolute_operation_residuals"].values()),
   "classification":"OPEN; finite point tangent is not a uniform cell bound",
   "joint_cell_forward_invariant":False,"all_time_magnetic_service_verified":False,
   "theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=calculate(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
