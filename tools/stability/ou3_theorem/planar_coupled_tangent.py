"""Homogeneous factor audit over the literal planar 20-s word.

FINITE DIAGNOSTIC ONLY. This does NOT propagate either coupled cross tangent.
The product F/(I-KH)/G omits dK*r, mean-dependent prediction/measurement
coefficients and exact quaternion-injection derivatives. No point norm is a
uniform cell bound or a complete nonlinear mean derivative.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import reset_matrix
from .planar_compatibility_quotient_mean import SCOPE,word

def diagnostic(path:Path,start=40000,end=44000):
    # Independent strict complete-word/paired-operand audit.
    checked,_,_=word(path,start,end)
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
    if not np.array_equal(M,checked):raise ValueError("homogeneous products disagree")
    s=np.linalg.svd(M,compute_uv=False)
    return {"qualification":"OU3_PLANAR_POINT_COUPLED_TANGENT_V1",
      "result_type":"FINITE DIAGNOSTIC ONLY","word_samples":[start,end],
      "operator_scope":SCOPE,"complete_nonlinear_mean_derivative_computed":False,
      "literal_mean_homogeneous_gain_2norm":float(s[0]),
      "literal_mean_homogeneous_leading_singular_values":s[:8].tolist(),
      "operation_counts":counts,
      "D_Px_point":None,"D_xP_point":None,
      "missing_for_cross_blocks":"linked dK*r, dH, dF/dBG, dQ/dBG, injection/reset differentials and branch/forcing enclosure",
      "uniform_b":None,"uniform_c":None,"uniform_q_P":None,"uniform_q_x":None,
      "joint_cell_forward_invariant":False,"all_time_magnetic_service_verified":False,
      "theorem_closed":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=diagnostic(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
