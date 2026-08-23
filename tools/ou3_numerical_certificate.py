#!/usr/bin/env python3
"""Numerical OU-III source/path certificate for the unchanged adaptive filter.

The certificate uses exact closed-loop maps emitted by the estimator.  It never
fits a transition matrix from the noisy truth trajectory and never defines a
Lyapunov matrix from empirical state covariance.

Claim levels are deliberately separate:
  * filter_regression: existing eight-sea RMS/quality gates;
  * exact_linear_source_certificate: group-compatible path LMIs on every
    executed ordinary-Live exact word;
  * numerical_certificate: linear + nonlinear/group + handoff/hybrid + noise;
  * deployment_theorem_certificate: validated continuous-source enclosure.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import struct
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
try:
    import cvxpy as cp
except ImportError:  # validation imports this module without solver dependency
    cp = None

from ou_sweep_common import PATTERNS, RECORDS

REPO = Path(__file__).resolve().parents[1]
TEST_DIR = REPO / "tests" / "kalman_ou_iii"
DEFAULT_DATA_DIR = REPO / "plots" / "kalman_ou_ii"
DEFAULT_OUT = REPO / "reports" / "results" / "ou3_numerical_certificate"
EPS = 1e-10
NX = 21
MAP_MAGIC = b"OU3MAP3\0"
MAP_VERSION = 3
MAP_PREFIX = struct.Struct("<ddIiii7f")
HORIZONS_S = (0.25, 0.5, 1.0, 2.0, 4.0)

# Numerical coordinate conditioning only; it never enters the estimator.
SCALE_ACTIVE = np.array(
    [1.0]*3 + [0.02]*3 + [5.0]*3 + [10.0]*3 +
    [100.0]*3 + [2.0]*3 + [0.2]*3, dtype=float
)
SCALE_HELD = SCALE_ACTIVE[:-3]
TAU_BINS = np.array([0.0,1.5,2.5,4.0,6.0,math.inf])
SIGMA_BINS = np.array([0.0,0.5,1.0,1.5,2.5,math.inf])
RS_BINS = np.array([0.0,0.25,2.0,10.0,25.0,80.0,math.inf])


@dataclass(frozen=True, order=True)
class SourceNode:
    mode: str
    mag_lock: int
    mag_refined: int
    tau_cell: int
    sigma_cell: int
    rs_cell: int

    def label(self) -> str:
        return (f"{self.mode}:m{self.mag_lock}:r{self.mag_refined}:"
                f"t{self.tau_cell}:s{self.sigma_cell}:q{self.rs_cell}")


@dataclass
class MapBlock:
    record: str
    index: int
    t0: float
    t1: float
    valid: bool
    hybrid_jump: bool
    start_live: bool
    end_live: bool
    start_active: bool
    end_active: bool
    start_lock: bool
    end_lock: bool
    start_refined: bool
    end_refined: bool
    acc_count: int
    mag_count: int
    pseudo_count: int
    tau0: float
    sigma0: float
    rs0: float
    tau1: float
    sigma1: float
    rs1: float
    linearization_residual: float
    phi: np.ndarray

    @property
    def mode(self) -> str:
        return "A" if self.start_active else "H"

    def start_node(self) -> SourceNode:
        return make_node(self.start_active,self.start_lock,self.start_refined,
                         self.tau0,self.sigma0,self.rs0)

    def end_node(self) -> SourceNode:
        return make_node(self.end_active,self.end_lock,self.end_refined,
                         self.tau1,self.sigma1,self.rs1)


@dataclass
class Word:
    record: str
    start_index: int
    end_index: int
    start_node: SourceNode
    end_node: SourceNode
    acc_count: int
    mag_count: int
    pseudo_count: int
    phi: np.ndarray
    prefix_norm_max: float

    def summary(self) -> dict:
        return {
            "record":self.record,"start_block":self.start_index,
            "end_block":self.end_index,"start":self.start_node.label(),
            "end":self.end_node.label(),"acc_count":self.acc_count,
            "mag_count":self.mag_count,"pseudo_count":self.pseudo_count,
        }


def _cell(v: float, bins: np.ndarray) -> int:
    return int(np.searchsorted(bins,v,side="right")-1)


def make_node(active: bool, lock: bool, refined: bool,
              tau: float, sigma: float, rs: float) -> SourceNode:
    return SourceNode("A" if active else "H",int(lock),int(refined),
                      _cell(tau,TAU_BINS),_cell(sigma,SIGMA_BINS),_cell(rs,RS_BINS))


def zyx_matrix(roll: float,pitch: float,yaw: float) -> np.ndarray:
    r,p,y=np.deg2rad([roll,pitch,yaw])
    cr,sr,cp,sp,cy,sy=math.cos(r),math.sin(r),math.cos(p),math.sin(p),math.cos(y),math.sin(y)
    return (np.array([[cy,-sy,0],[sy,cy,0],[0,0,1.]]) @
            np.array([[cp,0,sp],[0,1,0],[-sp,0,cp]]) @
            np.array([[1.,0,0],[0,cr,-sr],[0,sr,cr]]))


def so3_log(R: np.ndarray) -> np.ndarray:
    c=float(np.clip((np.trace(R)-1)*0.5,-1,1)); th=math.acos(c)
    vee=np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]],float)
    if th<1e-8: return 0.5*vee
    if math.pi-th<1e-5:
        A=0.5*(R+np.eye(3)); axis=np.sqrt(np.maximum(np.diag(A),0.0))
        s=np.sign(vee); s[s==0]=1; axis*=s; n=np.linalg.norm(axis)
        return th*(np.array([1.,0,0]) if n<EPS else axis/n)
    return th*vee/(2*math.sin(th))


def group_energy(rotvec: np.ndarray) -> float:
    return 1.0-math.cos(min(float(np.linalg.norm(rotvec)),math.pi))


def zu_to_ned(v: np.ndarray) -> np.ndarray:
    return np.array([v[1],v[0],-v[2]],float)


def output_csv_for(path: Path) -> Path:
    name=path.name
    if name.startswith("wave_data_"): name="w3d_"+name[len("wave_data_"):]
    stem=name[:-4] if name.endswith(".csv") else name
    return path.with_name(stem+"_fusion_ou3_cert.csv")


def parse_metrics(stdout: str) -> dict:
    ans={}
    for k,p in PATTERNS.items():
        m=p.search(stdout); ans[k]=float(m.group(1)) if m else None
    return ans


def load_columns(path: Path,names: list[str]) -> dict[str,np.ndarray]:
    values={n:[] for n in names}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            for n in names:
                try: values[n].append(float(row[n]))
                except (TypeError,ValueError): values[n].append(float("nan"))
    return {k:np.asarray(v,float) for k,v in values.items()}


def nearest_indices(t: np.ndarray,target: np.ndarray) -> np.ndarray:
    i=np.clip(np.searchsorted(t,target),0,len(t)-1); p=np.clip(i-1,0,len(t)-1)
    return np.where(np.abs(t[p]-target)<np.abs(t[i]-target),p,i)


def build_error_states(trace: np.ndarray,timeseries: Path) -> tuple[np.ndarray,np.ndarray]:
    cols=["time","roll_ref","pitch_ref","yaw_ref","roll_est","pitch_est","yaw_est",
          "disp_ref_x","disp_ref_y","disp_ref_z","vel_ref_x","vel_ref_y","vel_ref_z",
          "acc_ref_x","acc_ref_y","acc_ref_z","acc_bias_x","acc_bias_y","acc_bias_z",
          "gyro_bias_x","gyro_bias_y","gyro_bias_z"]
    d=load_columns(timeseries,cols); t=d["time"]
    ii=nearest_indices(t,np.asarray(trace["time_s"],float))
    pt=np.column_stack([d["disp_ref_y"],d["disp_ref_x"],-d["disp_ref_z"]])
    St=np.zeros_like(pt)
    if len(t)>1: St[1:]=np.cumsum(0.5*(pt[1:]+pt[:-1])*np.diff(t)[:,None],axis=0)
    E=np.zeros((len(trace),21),float); theta=np.zeros(len(trace))
    for k,j in enumerate(ii):
        rv=so3_log(zyx_matrix(d["roll_est"][j],d["pitch_est"][j],d["yaw_est"][j]) @
                   zyx_matrix(d["roll_ref"][j],d["pitch_ref"][j],d["yaw_ref"][j]).T)
        E[k,:3]=rv; theta[k]=np.linalg.norm(rv)
        bg=zu_to_ned(np.array([d["gyro_bias_x"][j],d["gyro_bias_y"][j],d["gyro_bias_z"][j]]))
        ba=zu_to_ned(np.array([d["acc_bias_x"][j],d["acc_bias_y"][j],d["acc_bias_z"][j]]))
        vt=zu_to_ned(np.array([d["vel_ref_x"][j],d["vel_ref_y"][j],d["vel_ref_z"][j]]))
        at=zu_to_ned(np.array([d["acc_ref_x"][j],d["acc_ref_y"][j],d["acc_ref_z"][j]]))
        E[k,3:6]=np.array([trace["bg_x"][k],trace["bg_y"][k],trace["bg_z"][k]])-bg
        E[k,6:9]=np.array([trace["v_x"][k],trace["v_y"][k],trace["v_z"][k]])-vt
        E[k,9:12]=np.array([trace["p_x"][k],trace["p_y"][k],trace["p_z"][k]])-pt[j]
        E[k,12:15]=np.array([trace["S_x"][k],trace["S_y"][k],trace["S_z"][k]])-St[j]
        E[k,15:18]=np.array([trace["aw_x"][k],trace["aw_y"][k],trace["aw_z"][k]])-at
        E[k,18:21]=np.array([trace["ba_x"][k],trace["ba_y"][k],trace["ba_z"][k]])-ba
    return E,theta


def load_exact_maps(path: Path,record: str) -> tuple[list[MapBlock],dict]:
    blocks=[]
    with path.open("rb") as f:
        if f.read(8)!=MAP_MAGIC: raise RuntimeError(f"bad exact-map magic: {path}")
        version,nx,stride=struct.unpack("<III",f.read(12))
        if version!=MAP_VERSION or nx!=NX: raise RuntimeError(f"unsupported exact-map format {version}/{nx}")
        i=0
        while True:
            p=f.read(MAP_PREFIX.size)
            if not p: break
            if len(p)!=MAP_PREFIX.size: raise RuntimeError(f"truncated map prefix: {path}")
            vals=MAP_PREFIX.unpack(p)
            raw=f.read(NX*NX*4)
            if len(raw)!=NX*NX*4: raise RuntimeError(f"truncated map matrix: {path}")
            t0,t1,flags,ac,mc,pc,tau0,sigma0,rs0,tau1,sigma1,rs1,resid=vals
            blocks.append(MapBlock(record,i,t0,t1,bool(flags&1),bool(flags&(1<<7)),
                bool(flags&(1<<1)),bool(flags&(1<<2)),bool(flags&(1<<3)),bool(flags&(1<<4)),
                bool(flags&(1<<5)),bool(flags&(1<<6)),bool(flags&(1<<8)),bool(flags&(1<<9)),
                ac,mc,pc,tau0,sigma0,rs0,tau1,sigma1,rs1,resid,
                np.frombuffer(raw,dtype="<f4").reshape(NX,NX).astype(float)))
            i+=1
    return blocks,{"version":version,"nx":nx,"stride_samples":stride,
                   "base_period_s":stride/200.0,"block_count":len(blocks)}


def normalized_map(phi: np.ndarray,mode: str) -> np.ndarray:
    dim=21 if mode=="A" else 18; scale=SCALE_ACTIVE if mode=="A" else SCALE_HELD
    A=phi[:dim,:dim]
    return A*scale[np.newaxis,:]/scale[:,np.newaxis]


def compose_words(blocks: list[MapBlock],mode: str,horizon_s: float) -> list[Word]:
    if not blocks: return []
    base=float(np.median([b.t1-b.t0 for b in blocks if b.t1>b.t0])); n=max(1,int(round(horizon_s/base)))
    words=[]
    for a in range(len(blocks)-n+1):
        seq=blocks[a:a+n]
        if any(not b.valid or b.hybrid_jump or b.mode!=mode for b in seq): continue
        if any(abs(seq[k+1].t0-seq[k].t1)>5e-3 for k in range(n-1)): continue
        Phi=np.eye(21); pre=1.0
        for b in seq:
            Phi=b.phi@Phi; pre=max(pre,float(np.linalg.norm(normalized_map(Phi,mode),2)))
        words.append(Word(seq[0].record,seq[0].index,seq[-1].index,seq[0].start_node(),seq[-1].end_node(),
                          sum(b.acc_count for b in seq),sum(b.mag_count for b in seq),
                          sum(b.pseudo_count for b in seq),normalized_map(Phi,mode),pre))
    return words


def generalized_lambda(phi: np.ndarray,Pi: np.ndarray,Pj: np.ndarray) -> float:
    Pi=0.5*(Pi+Pi.T); Pj=0.5*(Pj+Pj.T)
    lam,U=np.linalg.eigh(Pi)
    if np.min(lam)<=0: return math.inf
    invsqrt=U@np.diag(1/np.sqrt(lam))@U.T
    M=invsqrt@phi.T@Pj@phi@invsqrt
    return float(np.max(np.linalg.eigvalsh(0.5*(M+M.T))))


def group_compatible_numeric_metric(a_R: float,P_xi: np.ndarray) -> np.ndarray:
    """Small-angle quadratic matching a_R*(1-cos(theta))+xi'P_xi xi."""
    n=3+P_xi.shape[0]; M=np.zeros((n,n),float)
    M[:3,:3]=0.5*a_R*np.eye(3); M[3:,3:]=P_xi
    return M


def initial_design(words: list[Word],cap: int=160) -> list[int]:
    by_edge=defaultdict(list); score=[]
    for i,w in enumerate(words):
        by_edge[(w.start_node,w.end_node)].append(i); score.append(float(np.linalg.norm(w.phi,"fro")))
    chosen=set()
    for ids in by_edge.values(): chosen.update(sorted(ids,key=lambda i:score[i],reverse=True)[:2])
    chosen.update(sorted(range(len(words)),key=lambda i:score[i],reverse=True)[:cap])
    return sorted(chosen)


def solve_feasible(words: list[Word],design: list[int],dim: int,rho: float):
    """Solve the path LMI with the exact group-compatible local metric.

    The lifted theorem uses W=a_R V_R+xi'P_xi xi.  Since
    V_R=1-cos(theta)=theta^2/2+O(theta^4), its local quadratic is
    diag((a_R/2)I3,P_xi).  Enforcing this structure here prevents an easier
    full-state quadratic with attitude/xi cross terms from masquerading as the
    metric used by the nonlinear theorem.
    """
    if cp is None: return None,"CVXPY_UNAVAILABLE"
    nodes=sorted({words[i].start_node for i in design}|{words[i].end_node for i in design})
    if not nodes: return None,"NO_NODES"
    dxi=dim-3; I3=np.eye(3); Ix=np.eye(dxi); Z3x=np.zeros((3,dxi)); Zx3=np.zeros((dxi,3))
    a={n:cp.Variable(nonneg=True) for n in nodes}
    Px={n:cp.Variable((dxi,dxi),symmetric=True) for n in nodes}
    def Pexpr(n):
        return cp.bmat([[0.5*a[n]*I3,Z3x],[Zx3,Px[n]]])
    constraints=[]
    for n in nodes:
        constraints += [a[n]>=4e-5,Px[n]>>2e-5*Ix]
    constraints.append(cp.sum([cp.trace(Pexpr(n)) for n in nodes])==dim*len(nodes))
    I=np.eye(dim)
    for i in design:
        w=words[i]; Pi=Pexpr(w.start_node); Pj=Pexpr(w.end_node)
        constraints.append(w.phi.T@Pj@w.phi-rho*Pi << -1e-7*I)
    prob=cp.Problem(cp.Minimize(0),constraints)
    try:
        prob.solve(solver=cp.CLARABEL,verbose=False)
    except Exception:
        try: prob.solve(solver=cp.SCS,eps=2e-5,max_iters=12000,verbose=False)
        except Exception as exc: return None,f"SOLVER_ERROR:{type(exc).__name__}"
    if prob.status not in (cp.OPTIMAL,cp.OPTIMAL_INACCURATE): return None,str(prob.status)
    out={}
    for n in nodes:
        if a[n].value is None or Px[n].value is None: return None,"NO_VALUE"
        ax=float(a[n].value); X=0.5*(np.asarray(Px[n].value,float)+np.asarray(Px[n].value,float).T)
        M=group_compatible_numeric_metric(ax,X)
        if not np.all(np.isfinite(M)) or np.min(np.linalg.eigvalsh(M))<=1e-8: return None,"NON_SPD_VALUE"
        out[n]=M
    return out,str(prob.status)


def evaluate_metrics(words: list[Word],metrics: dict[SourceNode,np.ndarray]):
    vals=np.full(len(words),math.inf)
    for i,w in enumerate(words):
        if w.start_node in metrics and w.end_node in metrics:
            vals[i]=generalized_lambda(w.phi,metrics[w.start_node],metrics[w.end_node])
    return vals,int(np.argmax(vals)) if len(vals) else -1


def solve_path_metrics(words: list[Word],dim: int):
    if len(words)<4: return {"status":"INSUFFICIENT_WORDS","word_count":len(words),"linear_exact_replay_pass":False},None
    design=initial_design(words); target=0.9999; metrics=None; vals=None; status="not_run"
    for _ in range(6):
        metrics,status=solve_feasible(words,design,dim,target)
        if metrics is None: break
        vals,worst_i=evaluate_metrics(words,metrics)
        if np.all(np.isfinite(vals)) and float(np.max(vals))<1: break
        order=np.argsort(vals)[::-1]
        new=[int(i) for i in order if np.isfinite(vals[i]) and i not in design][:24]
        if not new: break
        design=sorted(set(design)|set(new))
    passed=metrics is not None and vals is not None and np.all(np.isfinite(vals)) and float(np.max(vals))<1
    rho_feasible=target if metrics is not None else None
    if not passed:
        for rho in (1.001,1.01,1.05,1.1,1.25,1.5,2.,4.,8.,16.):
            cand,s=solve_feasible(words,design,dim,rho)
            if cand is not None:
                metrics,status,rho_feasible=cand,s,rho; vals,worst_i=evaluate_metrics(words,metrics); break
    if vals is None: vals=np.full(len(words),math.inf); worst_i=0
    finite=vals[np.isfinite(vals)]; worst=float(np.max(vals)) if len(vals) else math.inf
    return {
        "status":status,"metric_structure":"diag((a_R/2) I3, P_xi)","word_count":len(words),
        "node_count":len({w.start_node for w in words}|{w.end_node for w in words}),
        "design_word_count":len(design),"rho_target":target,"rho_feasible_design":rho_feasible,
        "lambda_worst_exact_replay":worst,
        "lambda_p99_exact_replay":float(np.quantile(finite,0.99)) if len(finite) else None,
        "linear_exact_replay_pass":bool(passed),
        "worst_word":words[worst_i].summary() if words else None,
        "worst_prefix_euclidean_gain":float(max(w.prefix_norm_max for w in words)),
    },metrics


def choose_mode_certificate(all_blocks: dict[str,list[MapBlock]],mode: str):
    dim=21 if mode=="A" else 18; attempts=[]; best=None; best_metrics=None
    for h in HORIZONS_S:
        words=[]
        for blocks in all_blocks.values(): words.extend(compose_words(blocks,mode,h))
        r,m=solve_path_metrics(words,dim); r["horizon_s"]=h; attempts.append(r)
        if r.get("linear_exact_replay_pass"): best,best_metrics=r,m; break
        if m is not None and (best is None or r.get("lambda_worst_exact_replay",math.inf)<best.get("lambda_worst_exact_replay",math.inf)):
            best,best_metrics=r,m
    if best is None: best={"status":"NO_WORDS","linear_exact_replay_pass":False,"horizon_s":None}
    return {"mode":mode,"selected":best,"attempts":attempts},best_metrics


def record_linear_result(blocks,mode_result,metrics):
    h=mode_result["selected"].get("horizon_s")
    if h is None or metrics is None: return {"status":"NO_METRIC","word_count":0,"lambda_worst":None}
    words=compose_words(blocks,mode_result["mode"],float(h))
    if not words: return {"status":"NOT_EXERCISED","word_count":0,"lambda_worst":None}
    vals,worst=evaluate_metrics(words,metrics)
    return {"status":"PASS" if np.all(np.isfinite(vals)) and np.max(vals)<1 else "FAIL",
            "word_count":len(words),"lambda_worst":float(np.max(vals)),"worst_word":words[worst].summary()}


def handoff_hybrid(trace: np.ndarray,E: np.ndarray) -> dict:
    live=np.asarray(trace["live"],int); bias=np.asarray(trace["bias_active"],int)
    lock=np.asarray(trace["mag_lock"],int); refine=np.asarray(trace["mag_refined"],int)
    rising=lambda x:np.flatnonzero(x[1:]>x[:-1])+1
    li,bi,mi,ri=map(rising,(live,bias,lock,refine)); h=int(li[0]) if len(li) else (0 if len(live) and live[0] else -1)
    return {"live_handoff_time_s":float(trace["time_s"][h]) if h>=0 else None,
            "handoff_theta_deg":math.degrees(float(np.linalg.norm(E[h,:3]))) if h>=0 else None,
            "bias_release_time_s":float(trace["time_s"][bi[0]]) if len(bi) else None,
            "mag_lock_time_s":float(trace["time_s"][mi[0]]) if len(mi) else None,
            "mag_refine_time_s":float(trace["time_s"][ri[0]]) if len(ri) else None,
            "bias_release_events":int(len(bi)),"mag_lock_events":int(len(mi)),"mag_refine_events":int(len(ri))}


def stochastic_diagnostic(n_words: int) -> dict:
    d=9; wstar=6*math.sqrt(d); lo,hi=0.,wstar*wstar
    for _ in range(80):
        t=0.5*(lo+hi); q=d+2*math.sqrt(d*t)+2*t
        if q<=wstar*wstar: lo=t
        else: hi=t
    return {"gaussian_dimension":d,"w_star_normalized":wstar,"t_star":lo,
            "localization_union_bound":min(1.,max(1,n_words)*math.exp(-lo)),
            "qualification":"WAITING_FOR_METRIC_DEPENDENT_BW_VW_ENCLOSURE"}


def run_record(exe: Path,data: Path,out: Path):
    trace=(out/(data.stem+"_certificate_trace.csv")).resolve(); maps=(out/(data.stem+"_exact_maps.bin")).resolve()
    env=os.environ.copy(); env["OU3_CERT_TRACE"]=str(trace); env["OU3_CERT_MAP_TRACE"]=str(maps)
    env.setdefault("OU3_CERT_TRACE_STRIDE","10"); env.setdefault("OU3_CERT_MAP_STRIDE","50")
    env["W3D_WRITE_TIMESERIES"]="1"; env["W3D_VALIDATION_WINDOW_SEC"]="900"
    p=subprocess.run([str(exe),"--input",str(data)],cwd=data.parent,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    ts=output_csv_for(data)
    if not trace.exists() or not maps.exists() or not ts.exists(): raise RuntimeError(f"missing certificate outputs for {data.name}\n{p.stdout[-5000:]}")
    return trace,maps,ts,parse_metrics(p.stdout),p.returncode==0,p.stdout


def markdown(report: dict) -> str:
    lin=report["exact_linear_source_certificate"]
    x=["# OU-III exact-map numerical stability certificate","",
       f"Filter regression: **{report['filter_regression']}**",
       f"Exact executed-word linear source/path certificate: **{lin['status']}**",
       f"Numerical full source-funnel certificate: **{report['numerical_certificate']}**",
       f"Deployment theorem certificate: **{report['deployment_theorem_certificate']}**","",
       "Path metrics are constrained to `diag((a_R/2) I3, P_xi)`, the local quadratic of the exact group metric.","",
       f"Held: horizon {lin['held']['selected'].get('horizon_s')} s, worst lambda {lin['held']['selected'].get('lambda_worst_exact_replay')}",
       f"Active: horizon {lin['active']['selected'].get('horizon_s')} s, worst lambda {lin['active']['selected'].get('lambda_worst_exact_replay')}",
       f"Maximum map reconstruction residual: {report['map_integrity']['max_linearization_residual']:.3e}","",
       "| Sea | RMS | exact linear | max theta | handoff theta | invalid blocks |","|---|---:|---:|---:|---:|---:|"]
    for r in report["records"]:
        x.append(f"| {r['family']} {r['Hs_m']:.2f} | {'PASS' if r['rms_regression_pass'] else 'FAIL'} | "
                 f"{'PASS' if r['exact_linear_pass'] else 'FAIL'} | {r['theta_max_deg']:.2f} deg | "
                 f"{r['handoff_hybrid']['handoff_theta_deg']} | {r['invalid_map_blocks']} |")
    x += ["","The exact-linear status covers executed ordinary-Live words only. Full numerical certification still requires the SO(3) nonlinear margin, source-shaped handoff/invariant funnel, hybrid jumps and metric-dependent stochastic constants. Continuous-source theorem status remains separate."]
    return "\n".join(x)


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--data-dir",type=Path,default=DEFAULT_DATA_DIR); ap.add_argument("--output-dir",type=Path,default=DEFAULT_OUT)
    ap.add_argument("--sim",type=Path,default=TEST_DIR/"ou3-certificate-sim"); ap.add_argument("--no-build",action="store_true"); args=ap.parse_args()
    if cp is None: raise RuntimeError("cvxpy is required for path-LMI certificate")
    data_dir=args.data_dir.resolve(); out=args.output_dir.resolve(); exe=args.sim.resolve(); out.mkdir(parents=True,exist_ok=True); (out/"logs").mkdir(exist_ok=True)
    if not args.no_build: subprocess.run(["make","-C",str(TEST_DIR),"ou3-certificate-sim"],check=True)
    raw=[]; all_blocks={}; meta={}
    for family,hs,name in RECORDS:
        data=(data_dir/name).resolve()
        if not data.exists(): raise FileNotFoundError(data)
        tr,mp,ts,met,ok,log=run_record(exe,data,out); slug=f"{family.lower().replace('-','_')}_{hs:.2f}".replace(".","_")
        (out/"logs"/f"{slug}.log").write_text(log)
        trace=np.genfromtxt(tr,delimiter=",",names=True,dtype=None,encoding=None); E,theta=build_error_states(trace,ts); blocks,m=load_exact_maps(mp,slug)
        all_blocks[slug]=blocks; meta[slug]=m; raw.append((family,hs,slug,trace,E,theta,met,ok,blocks))
    held,Ph=choose_mode_certificate(all_blocks,"H"); active,Pa=choose_mode_certificate(all_blocks,"A")
    linear=bool(held["selected"].get("linear_exact_replay_pass") and active["selected"].get("linear_exact_replay_pass"))
    arrays={}
    for prefix,P in (("H",Ph),("A",Pa)):
        if P:
            labels=sorted(P); arrays[f"{prefix}_labels"]=np.asarray([n.label() for n in labels]); arrays[f"{prefix}_P_local"]=np.stack([P[n] for n in labels]); arrays[f"{prefix}_a_R"]=np.asarray([2*P[n][0,0] for n in labels])
    if arrays: np.savez_compressed(out/"path_metrics.npz",**arrays)
    records=[]
    for family,hs,slug,trace,E,theta,met,ok,blocks in raw:
        hr=record_linear_result(blocks,held,Ph); ar=record_linear_result(blocks,active,Pa); ex=[x for x in (hr,ar) if x["status"]!="NOT_EXERCISED"]
        residual=max((b.linearization_residual for b in blocks),default=math.nan)
        records.append({"family":family,"Hs_m":hs,"slug":slug,"rms_regression_pass":bool(ok),"rms_metrics":met,
            "theta_max_deg":math.degrees(float(np.max(theta))),"group_energy_max":float(max(group_energy(x[:3]) for x in E)),
            "all_attitude_inside_pi":bool(np.max(theta)<math.pi),"held_linear":hr,"active_linear":ar,
            "exact_linear_pass":bool(ex) and all(x["status"]=="PASS" for x in ex),"map_blocks":len(blocks),
            "invalid_map_blocks":sum(not b.valid for b in blocks),"hybrid_map_blocks":sum(b.hybrid_jump for b in blocks),
            "max_linearization_residual":residual,"handoff_hybrid":handoff_hybrid(trace,E),"stochastic":stochastic_diagnostic(len(blocks))})
    rms=all(r["rms_regression_pass"] for r in records); maxres=max(r["max_linearization_residual"] for r in records); integrity=math.isfinite(maxres) and maxres<5e-3
    lstatus="PASS" if linear and integrity else "FAIL"
    report={"schema":4,"scope":"eight_noisy_reference_replays_exact_filter_maps_group_compatible_metrics","record_count":len(records),
        "filter_regression":"PASS" if rms else "FAIL","map_integrity":{"max_linearization_residual":maxres,"pass":integrity,"per_record":meta},
        "exact_linear_source_certificate":{"status":lstatus,"metric_structure":"diag((a_R/2) I3, P_xi)","held":held,"active":active},
        "numerical_certificate":"BLOCKED_AFTER_LINEAR" if lstatus=="PASS" else "FAIL_AT_LINEAR_GATE",
        "numerical_missing":["nodewise exact SO(3) theta_star","exact nonlinear word infimum mu_W","handoff/invariant funnel c0,b,N_H,T_H","held-active/magnetic/tilt-cooldown jump inequalities","metric-dependent non-empirical b_W,v_W"],
        "deployment_theorem_certificate":"NOT_ESTABLISHED","deployment_missing":"validated continuous-source enclosure","records":records}
    (out/"certificate.json").write_text(json.dumps(report,indent=2,sort_keys=True)); text=markdown(report); (out/"certificate.md").write_text(text); print(text)
    return 1 if (not rms or not integrity) else 0


if __name__=="__main__": raise SystemExit(main())
