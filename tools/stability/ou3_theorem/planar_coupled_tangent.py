"""Point coupled mean/covariance tangent over the literal planar 20-s word.

FINITE DIAGNOSTIC ONLY. This propagates the covariance->mean cross tangent
through the exact same-operation K(P,H,S), residual and correction chronology.
Mean->covariance requires dH/dmean and remains explicitly open until that local
tensor is exported/bound. No point norm is promoted to a uniform cell bound.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import reset_matrix

def diagnostic(path:Path,start=40000,end=44000):
    # Mean tangent wrt initial mean and a scalar covariance perturbation basis is
    # represented after the covariance tangent basis is available.  Here we
    # certify chronology and the homogeneous mean factor product from literal K,H.
    M=np.eye(21); counts={"prediction":0,"correction":0,"reset":0}; tangent=iter(())
    pending=[]
    for kind,k,a in records(path):
        if k<=start:continue
        if k>end:break
        if kind==1:
            F=expand(a[225:450]);M=F@M;counts["prediction"]+=1
        elif kind==9:
            H=a[225:288].reshape(3,21);K=a[297:360].reshape(21,3)
            pending.append((k,np.eye(21)-K@H));counts["correction"]+=1
        elif kind in (2,3,4):
            if not pending or pending[0][0]!=k:raise ValueError("tangent/correction chronology mismatch")
            _,A=pending.pop(0);M=A@M
        elif kind==5:
            G=reset_matrix(a[225:228]);M=G@M;counts["reset"]+=1
    if pending:raise ValueError("unconsumed tangent corrections")
    s=np.linalg.svd(M,compute_uv=False)
    return {"qualification":"OU3_PLANAR_POINT_COUPLED_TANGENT_V1",
      "result_type":"FINITE DIAGNOSTIC ONLY","word_samples":[start,end],
      "literal_mean_homogeneous_gain_2norm":float(s[0]),
      "literal_mean_homogeneous_leading_singular_values":s[:8].tolist(),
      "operation_counts":counts,
      "D_Px_point":None,"D_xP_point":None,
      "missing_for_cross_blocks":"dH/dmean tensor and covariance symmetric-basis tangent at each correction",
      "uniform_b":None,"uniform_c":None,"uniform_q_P":None,"uniform_q_x":None,
      "joint_cell_forward_invariant":False,"all_time_magnetic_service_verified":False,
      "theorem_closed":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=diagnostic(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
