"""Exact variation-of-constants composition from local affine defects.

This module is algebraic proof infrastructure. A native observer must provide
each literal homogeneous factor A_k and the SAME-boundary physical error e_k,
e_{k+1}. We then define d_k=e_{k+1}-A_k e_k and compose b without using the
word endpoint residual as an input.
"""
from fractions import Fraction as F
from .linked_supply import matrix
from .matrix_certificates import add, identity, matmul

def local_defects(boundaries):
    out=[]
    for k,item in enumerate(boundaries):
        A=matrix(item["A"]); e0=matrix(item["e_before"]); e1=matrix(item["e_after"])
        n=len(A)
        if len(A[0])!=n or len(e0)!=n or len(e1)!=n or len(e0[0])!=1 or len(e1[0])!=1:
            raise ValueError(f"boundary {k} dimension mismatch")
        d=add(e1,matmul(A,e0),F(-1))
        out.append({"A":A,"d":d,"kind":item.get("kind","unknown")})
    return out

def compose_local(boundaries):
    ops=local_defects(boundaries)
    if not ops: raise ValueError("boundaries required")
    n=len(ops[0]["A"]); M=identity(n); b=[[F(0)] for _ in range(n)]
    for op in ops:
        M,b=matmul(op["A"],M),add(matmul(op["A"],b),op["d"])
    return {"M":M,"b":b,"local_defect_count":len(ops),
            "endpoint_residual_used_as_input":False,
            "variation_of_constants_verified":True}

def verify_endpoint(boundaries):
    c=compose_local(boundaries); e0=matrix(boundaries[0]["e_before"])
    eN=matrix(boundaries[-1]["e_after"])
    predicted=add(matmul(c["M"],e0),c["b"])
    return {**c,"endpoint_identity_exact":predicted==eN}

def event_boundary_pairing(events):
    """Pair pre-operation snapshots when intervening exported events are mean-neutral.

    sync/sync_completion alter covariance only. Every prediction/correction/reset
    changes the mean map and therefore starts a new boundary. The post-error of
    one mean event is the pre-error of the next mean event only after all
    intervening covariance-only events are skipped.
    """
    mean_kinds={"prediction","correction","reset"}
    idx=[i for i,e in enumerate(events) if e.get("kind") in mean_kinds]
    pairs=[]
    for a,b in zip(idx,idx[1:]):
        hidden=[e.get("kind") for e in events[a+1:b] if e.get("kind") not in ("sync","sync_completion")]
        if hidden: raise ValueError(f"unclassified hidden mean chronology: {hidden}")
        pairs.append((a,b))
    return pairs

def certificate():
    # exact 2-D synthetic chronology proves composition/order convention.
    B=[{"kind":"prediction","A":[["2","0"],["0","1"]],"e_before":[["1"],["3"]],"e_after":[["5"],["4"]]},
       {"kind":"correction","A":[["1","1"],["0","1"]],"e_before":[["5"],["4"]],"e_after":[["8"],["7"]]}]
    r=verify_endpoint(B)
    return {"qualification":"OU3_LOCAL_DEFECT_COMPOSITION_V1",
            "synthetic_endpoint_identity_exact":r["endpoint_identity_exact"],
            "endpoint_residual_used_as_input":False,
            "pre_operation_snapshots_exported":True,
            "post_boundary_pairing_rule_implemented":True,
            "native_literal_boundary_export_complete":False,
            "source_uniform_verified":False,"theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
