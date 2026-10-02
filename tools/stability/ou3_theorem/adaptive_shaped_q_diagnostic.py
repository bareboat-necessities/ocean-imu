"""Generated-state, fully coupled LIN shaped-storage diagnostic.

No common-Q requirement and no v/p/S/aw block separation.  At each endpoint
the 12-state normalized LIN storage is allowed to depend on the actual generated
shipping tuner tuple through the literal carried covariance metric. Diagnostic
only; no source-uniform promotion.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .adaptive_shaped_storage import lin_scaling21

def adaptive_events(events):
    return [e for e in events if e.get("kind")=="adaptive_state"]

def tuple_at(events,t):
    a=[e for e in adaptive_events(events) if float(e["physical_t"])<=t+1e-9]
    if not a: raise ValueError("no causal generated tuner state")
    return a[-1]

def literal_factor(e):
    k=e["kind"]
    if k=="prediction":
        A=np.eye(21);A[:6,:6]=np.asarray(e["F_AG"],float);A[6:18,6:18]=np.asarray(e["F_LIN"],float);A[18:,18:]=float(e["phi_BA"])*np.eye(3);return A
    if k=="correction": return np.eye(21)-np.asarray(e["K"],float)@np.asarray(e["H"],float)
    if k=="reset":
        d=np.asarray(e["d"],float).reshape(3);S=np.array([[0,-d[2],d[1]],[d[2],0,-d[0]],[-d[1],d[0],0.]])
        A=np.eye(21);A[:3,:3]+=S/2;return A
    raise ValueError(k)

def normalized_product(trace):
    M=np.eye(21); count=0
    for e in trace["events"]:
        if e.get("kind") not in ("prediction","correction","reset"): continue
        u=tuple_at(trace["events"],float(e["physical_t"]))
        D=lin_scaling21(float(u["tau"]),float(u["sigma_aw"]))
        M=(D@literal_factor(e)@np.linalg.inv(D))@M;count+=1
    return M,count

def generated_Q(P,tau,sigma):
    # Exact covariance metric transported into shipping-normalized coordinates.
    # Its LIN block is fully dense 12x12 and changes with the generated state.
    D=lin_scaling21(tau,sigma)
    Pn=D@np.asarray(P,float)@D.T
    Jn=np.linalg.inv((Pn+Pn.T)/2)
    return Jn

def analyze(path):
    tr=json.loads(Path(path).read_text()); ae=adaptive_events(tr["events"])
    if not ae: raise ValueError("adaptive state trace absent")
    u0=ae[0];uN=ae[-1]
    Q0=generated_Q(tr["root_covariance"],float(u0["tau"]),float(u0["sigma_aw"]))
    QN=generated_Q(tr["terminal_covariance"],float(uN["tau"]),float(uN["sigma_aw"]))
    M,n=normalized_product(tr)
    L=np.linalg.cholesky(Q0)
    C=np.linalg.solve(L,M.T@QN@M)@np.linalg.inv(L.T)
    rho=float(np.linalg.eigvalsh((C+C.T)/2).max())
    lin=QN[6:18,6:18]
    off=lin-np.diag(np.diag(lin))
    return {"profile":tr["profile"],"horizon_s":tr.get("horizon_s"),"operations":n,
      "generated_Q_endpoint_rho":rho,"strict_endpoint_contraction":rho<1,
      "Q0_LIN_condition":float(np.linalg.cond(Q0[6:18,6:18])),
      "QN_LIN_condition":float(np.linalg.cond(lin)),
      "QN_LIN_offdiagonal_frobenius_fraction":float(np.linalg.norm(off)/np.linalg.norm(lin)),
      "root_generated_tuple":{k:u0[k] for k in ("tau","sigma_aw","R_S","T_S")},
      "terminal_generated_tuple":{k:uN[k] for k in ("tau","sigma_aw","R_S","T_S")},
      "structures_preserved":["fully coupled 12x12 LIN storage","generated-state-dependent endpoint storage","literal covariance metric","literal gains and chronology","joint generated tuner history","persistent physical/filter history"],
      "relaxations_introduced":["finite carried replay only; source-uniform generated-state storage family not yet enclosed"],
      "failure_class_if_not_strict":"D_SUFFICIENT_BOUND_FAILURE","theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--directory",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--horizons",default="60,100");a=p.parse_args()
 out={"kind":"NON_PROMOTING_GENERATED_FULL_LIN_STORAGE","words":[analyze(a.directory/f"word-{prof}-{h}s.json") for h in map(int,a.horizons.split(",")) for prof in ("0","wave")]}
 a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()
