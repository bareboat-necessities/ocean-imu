"""Exact variation-of-constants composition from local affine defects.\n

This module is algebraic proof infrastructure. A native observer must provide
each literal homogeneous factor A_k and the SAME-boundary physical error e_k,
e_{k+1}. We then define d_k=e_{k+1}-A_k e_k and compose b without using the
word endpoint residual as an input.
"""
from fractions import Fraction as F
from .linked_supply import matrix
from .matrix_certificates import add, identity, matmul

def quat_mul(a,b):
    # Stored Eigen coeff order [x,y,z,w].
    ax,ay,az,aw=map(float,a); bx,by,bz,bw=map(float,b)
    return [aw*bx+ax*bw+ay*bz-az*by,
            aw*by-ay*0+ay*bw+az*bx-ax*bz,
            aw*bz+az*bw+ax*by-ay*bx,
            aw*bw-ax*bx-ay*by-az*bz]

def quat_conj(q):
    x,y,z,w=map(float,q); return [-x,-y,-z,w]

def left_attitude_error(est_wb,true_bw):
    """Rotation vector d with true W->B = Exp(d) * estimated W->B.

    Inputs use Eigen quaternion coeff order [x,y,z,w]. Carried fixture has B'=B.
    """
    import math
    true_wb=quat_conj(true_bw)
    dq=quat_mul(true_wb,quat_conj(est_wb))
    n=math.sqrt(sum(x*x for x in dq))
    dq=[x/n for x in dq]
    if dq[3]<0: dq=[-x for x in dq]
    v=math.sqrt(sum(x*x for x in dq[:3]))
    if v<1e-14: return [2*x for x in dq[:3]]
    ang=2*math.atan2(v,dq[3])
    return [ang*x/v for x in dq[:3]]

def _col(v):
    return [float(x[0] if isinstance(x,(list,tuple)) else x) for x in v]

def carried_error(event, after=False):
    """21-state physical error in the literal local error-state coordinates."""
    suffix="_after" if after else ""
    def physical(name):
        key=name+suffix
        return _col(event[key] if after and key in event else event[name])
    x=_col(event["estimator_state"+suffix])
    th=left_attitude_error(_col(event["estimator_quaternion"+suffix]),physical("physical_quaternion"))
    th=[a-b for a,b in zip(th,x[:3])]
    bg=[a-b for a,b in zip(physical("physical_bg"),x[3:6])]
    v=[a-b for a,b in zip(physical("physical_v"),x[6:9])]
    p=[a-b for a,b in zip(physical("physical_p"),x[9:12])]
    S=[a-b for a,b in zip(physical("physical_S"),x[12:15])]
    aw=[a-b for a,b in zip(physical("physical_a"),x[15:18])]
    ba=[a-b for a,b in zip(physical("physical_ba"),x[18:21])]
    return [[str(z)] for z in th+bg+v+p+S+aw+ba]

def literal_mean_factor(event):
    kind=event["kind"]; A=[[F(int(i==j)) for j in range(21)] for i in range(21)]
    if kind=="prediction":
        for off,key in ((0,"F_AG"),(6,"F_LIN")):
            B=matrix(event[key])
            for i,row in enumerate(B):
                A[off+i][off:off+len(row)]=row
        phi=F(event["phi_BA"])
        for i in range(18,21): A[i][i]=phi
    elif kind=="correction":
        H=matrix(event["H"]); K=matrix(event["K"])
        A=add(A,matmul(K,H),F(-1))
    elif kind=="reset":
        x,y,z=[F(row[0]) for row in event["d"]]
        cross=[[0,-z,y],[z,0,-x],[-y,x,0]]
        for i in range(3):
            for j in range(3): A[i][j]+=F(cross[i][j],2)
    else:
        raise ValueError("mean-neutral event has no mean factor")
    return A

def carried_boundaries(events):
    """Use literal same-operation pre/post snapshots; never the next event."""
    out=[]
    for event in events:
        if event.get("kind") not in {"prediction","correction","reset"}:
            continue
        if "estimator_state_after" not in event or "estimator_quaternion_after" not in event:
            raise ValueError("mean event missing same-boundary post snapshot")
        out.append({"kind":event["kind"],"sensor":event.get("sensor"),"A":literal_mean_factor(event),
                    "e_before":carried_error(event),
                    "e_after":carried_error(event, after=True)})
    return out

def carried_local_defect_certificate(trace):
    B=carried_boundaries(trace["events"])
    if not B: raise ValueError("no carried mean boundaries")
    c=compose_local(B)
    e0=matrix(B[0]["e_before"]); eN=matrix(B[-1]["e_after"])
    residual=add(eN,matmul(c["M"],e0),F(-1))
    diff=add(c["b"],residual,F(-1))
    maxdiff=max(abs(float(x[0])) for x in diff)
    gaps=[]
    for k,(left,right) in enumerate(zip(B,B[1:])):
        gap=add(matrix(right["e_before"]),matrix(left["e_after"]),F(-1))
        gaps.append((max(abs(float(x[0])) for x in gap),k,left["kind"],right["kind"]))
    worst=max(gaps,default=(0.0,-1,"none","none"))
    return {"local_defect_count":c["local_defect_count"],
            "endpoint_residual_used_as_input":False,
            "local_b_endpoint_residual_max_abs":maxdiff,
            "local_b_endpoint_parity":maxdiff < 5e-10,
            "maximum_interoperation_error_gap":worst[0],
            "worst_gap_after_operation":worst[1],
            "worst_gap_kind_pair":[worst[2],worst[3]],
            "native_literal_boundary_export_complete":True,
            "same_operation_pre_post_snapshots":True,
            "b_local":[[float(x[0])] for x in c["b"]]}

def local_defects(boundaries):
    out=[]
    for k,item in enumerate(boundaries):
        A=matrix(item["A"]); e0=matrix(item["e_before"]); e1=matrix(item["e_after"])
        n=len(A)
        if len(A[0])!=n or len(e0)!=n or len(e1)!=n or len(e0[0])!=1 or len(e1[0])!=1:
            raise ValueError(f"boundary {k} dimension mismatch")
        d=add(e1,matmul(A,e0),F(-1))
        out.append({"A":A,"d":d,"kind":item.get("kind","unknown"),"sensor":item.get("sensor")})
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
            "carried_driver_truth_reconstructible":True,
            "source_uniform_truth_boundary_export_complete":False,
            "source_uniform_verified":False,"theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
