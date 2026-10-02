"""Common-Q carried feasibility for shipping-normalized OU-III factors.

Diagnostic only.  Uses actual exported gains and generated tuner tuples.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .adaptive_shaped_storage import lin_scaling21

def tuple_at(events,t):
    a=[e for e in events if e.get("kind")=="adaptive_state" and float(e["physical_t"])<=t+1e-9]
    if not a: raise ValueError("no causal adaptive state")
    return a[-1]

def normalized_ops(trace):
    out=[]
    for e in trace["events"]:
        k=e.get("kind")
        if k=="adaptive_state" or k in ("sync","sync_completion"): continue
        t=float(e["physical_t"]); u=tuple_at(trace["events"],t)
        D=lin_scaling21(float(u["tau"]),float(u["sigma_aw"]))
        if k=="prediction":
            A=np.eye(21);A[:6,:6]=np.array(e["F_AG"]);A[6:18,6:18]=np.array(e["F_LIN"]);A[18:,18:]=float(e["phi_BA"])*np.eye(3)
        elif k=="correction":
            A=np.eye(21)-np.array(e["K"])@np.array(e["H"])
        elif k=="reset":
            d=np.array(e["d"],float).reshape(3);S=np.array([[0,-d[2],d[1]],[d[2],0,-d[0]],[-d[1],d[0],0.]])
            A=np.eye(21);A[:3,:3]+=S/2
        else: raise ValueError(k)
        out.append((D@A@np.linalg.inv(D),k,e.get("sensor"),u))
    return out

def product(ops):
    M=np.eye(21)
    for A,*_ in ops:M=A@M
    return M

def solve_q(M):
    # Block-diagonal shaped storage: AG,BG, four normalized LIN chain blocks, BA.
    # This restriction is a sufficient-proof ansatz only.
    groups=[slice(0,3),slice(3,6),slice(6,9),slice(9,12),slice(12,15),slice(15,18),slice(18,21)]
    def Q(z):
        q=np.zeros((21,21))
        for w,g in zip(np.exp(z),groups):q[g,g]=w*np.eye(3)
        return q
    def obj(z):
        q=Q(z); L=np.linalg.cholesky(q); C=np.linalg.solve(L,M.T@q@M)@np.linalg.inv(L.T)
        return np.linalg.eigvalsh((C+C.T)/2).max()+1e-8*np.dot(z,z)
    res=minimize(obj,np.zeros(7),method="BFGS",options={"maxiter":2000,"gtol":1e-10})
    q=Q(res.x); L=np.linalg.cholesky(q); C=np.linalg.solve(L,M.T@q@M)@np.linalg.inv(L.T)
    return res,float(np.linalg.eigvalsh((C+C.T)/2).max()),np.exp(res.x)

def analyze(path):
    tr=json.loads(Path(path).read_text());ops=normalized_ops(tr);M=product(ops);res,rho,w=solve_q(M)
    S=[o for o in ops if o[1]=="correction" and o[2]=="S"]
    return {"profile":tr["profile"],"horizon_s":tr.get("horizon_s"),"operations":len(ops),"S_corrections":len(S),
      "common_block_Q_optimizer_success":bool(res.success),"common_block_Q_rho":rho,"common_block_Q_weights":w.tolist(),
      "strict_endpoint_contraction":rho<1,
      "structures_preserved":["literal exported gains","joint generated tuner tuples","literal prediction/correction/reset chronology","persistent history"],
      "relaxations_introduced":["block-diagonal common-Q ansatz","finite carried replay only"],"failure_class_if_not_strict":"D_SUFFICIENT_BOUND_FAILURE","theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--directory",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--horizons",default="60,100");a=p.parse_args()
 out={"kind":"NON_PROMOTING_SHAPED_Q_FEASIBILITY","words":[analyze(a.directory/f"word-{prof}-{h}s.json") for h in map(int,a.horizons.split(",")) for prof in ("0","wave")]}
 a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()
