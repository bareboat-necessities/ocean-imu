#!/usr/bin/env python3
"""Numerically instantiate the OU-III source/path stability certificate.

This tool deliberately separates three questions:

* filter_regression_pass: do the unchanged eight noisy reference simulations pass
  the existing performance gates?
* exact_linear_replay_certificate: do the exact closed-loop error maps executed
  by those runs admit source/path Lyapunov metrics with rho < 1?
* deployment_theorem_certificate: have the continuous source cells, nonlinear
  word maps, funnel jumps and stochastic constants been rigorously enclosed?

Unlike the first version of this tool, no state-transition matrix is identified
from noisy error trajectories and no Lyapunov matrix is defined as an inverse
empirical covariance.  The C++ certificate executable emits exact composed
closed-loop maps from the estimator's own prediction and Kalman update matrices;
this script solves path LMIs for those maps and then evaluates every executed
word against the solved metrics.
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
from typing import Iterable

import numpy as np

try:
    import cvxpy as cp
except ImportError:  # pragma: no cover - CI installs python3-cvxpy
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

# Coordinate normalization only.  The estimator never sees these scales.
SCALE_ACTIVE = np.array(
    [1.0]*3 + [0.02]*3 + [5.0]*3 + [10.0]*3 +
    [100.0]*3 + [2.0]*3 + [0.2]*3, dtype=float
)
SCALE_HELD = SCALE_ACTIVE[:-3]
TAU_BINS = np.array([0.0, 1.5, 2.5, 4.0, 6.0, math.inf])
SIGMA_BINS = np.array([0.0, 0.5, 1.0, 1.5, 2.5, math.inf])
RS_BINS = np.array([0.0, 0.25, 2.0, 10.0, 25.0, 80.0, math.inf])


@dataclass(frozen=True, order=True)
class SourceNode:
    mode: str
    mag_lock: int
    mag_refined: int
    tau_cell: int
    sigma_cell: int
    rs_cell: int

    def label(self) -> str:
        return (
            f"{self.mode}:m{self.mag_lock}:r{self.mag_refined}:"
            f"t{self.tau_cell}:s{self.sigma_cell}:q{self.rs_cell}"
        )


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
        return make_node(self.start_active, self.start_lock, self.start_refined,
                         self.tau0, self.sigma0, self.rs0)

    def end_node(self) -> SourceNode:
        return make_node(self.end_active, self.end_lock, self.end_refined,
                         self.tau1, self.sigma1, self.rs1)


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
            "record": self.record,
            "start_block": self.start_index,
            "end_block": self.end_index,
            "start": self.start_node.label(),
            "end": self.end_node.label(),
            "acc_count": self.acc_count,
            "mag_count": self.mag_count,
            "pseudo_count": self.pseudo_count,
        }


def _cell(v: float, bins: np.ndarray) -> int:
    return int(np.searchsorted(bins, v, side="right") - 1)


def make_node(active: bool, lock: bool, refined: bool,
              tau: float, sigma: float, rs: float) -> SourceNode:
    return SourceNode(
        "A" if active else "H", int(lock), int(refined),
        _cell(float(tau), TAU_BINS),
        _cell(float(sigma), SIGMA_BINS),
        _cell(float(rs), RS_BINS),
    )


def zyx_matrix(roll: float, pitch: float, yaw: float) -> np.ndarray:
    r, p, y = np.deg2rad([roll, pitch, yaw])
    cr, sr, cp, sp, cy, sy = (
        math.cos(r), math.sin(r), math.cos(p), math.sin(p), math.cos(y), math.sin(y)
    )
    Rx = np.array([[1,0,0],[0,cr,-sr],[0,sr,cr]], float)
    Ry = np.array([[cp,0,sp],[0,1,0],[-sp,0,cp]], float)
    Rz = np.array([[cy,-sy,0],[sy,cy,0],[0,0,1]], float)
    return Rz @ Ry @ Rx


def so3_log(R: np.ndarray) -> np.ndarray:
    c = float(np.clip((np.trace(R)-1.0)*0.5, -1.0, 1.0))
    th = math.acos(c)
    vee = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]], float)
    if th < 1e-8:
        return 0.5 * vee
    if math.pi - th < 1e-5:
        A = 0.5 * (R + np.eye(3))
        axis = np.sqrt(np.maximum(np.diag(A), 0.0))
        signs = np.sign(vee); signs[signs == 0] = 1
        axis *= signs
        n = np.linalg.norm(axis)
        return th * (np.array([1.0,0,0]) if n < EPS else axis/n)
    return th * vee / (2.0 * math.sin(th))


def group_energy(rotvec: np.ndarray) -> float:
    return 1.0 - math.cos(min(float(np.linalg.norm(rotvec)), math.pi))


def zu_to_ned(v: np.ndarray) -> np.ndarray:
    return np.array([v[1], v[0], -v[2]], float)


def output_csv_for(input_path: Path) -> Path:
    name = input_path.name
    if name.startswith("wave_data_"):
        name = "w3d_" + name[len("wave_data_"):]
    stem = name[:-4] if name.endswith(".csv") else name
    return input_path.with_name(stem + "_fusion_ou3_cert.csv")


def parse_metrics(stdout: str) -> dict[str, float | None]:
    ans = {}
    for k, p in PATTERNS.items():
        m = p.search(stdout)
        ans[k] = float(m.group(1)) if m else None
    return ans


def load_columns(path: Path, names: list[str]) -> dict[str, np.ndarray]:
    values = {n: [] for n in names}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            for n in names:
                try:
                    values[n].append(float(row[n]))
                except (TypeError, ValueError):
                    values[n].append(float("nan"))
    return {k: np.asarray(v, float) for k,v in values.items()}


def nearest_indices(t: np.ndarray, target: np.ndarray) -> np.ndarray:
    i = np.clip(np.searchsorted(t, target), 0, len(t)-1)
    p = np.clip(i-1, 0, len(t)-1)
    return np.where(np.abs(t[p]-target) < np.abs(t[i]-target), p, i)


def build_error_states(trace: np.ndarray, timeseries: Path) -> tuple[np.ndarray, np.ndarray]:
    cols = [
        "time","roll_ref","pitch_ref","yaw_ref","roll_est","pitch_est","yaw_est",
        "disp_ref_x","disp_ref_y","disp_ref_z","vel_ref_x","vel_ref_y","vel_ref_z",
        "acc_ref_x","acc_ref_y","acc_ref_z",
        "acc_bias_x","acc_bias_y","acc_bias_z","gyro_bias_x","gyro_bias_y","gyro_bias_z",
    ]
    d = load_columns(timeseries, cols)
    t = d["time"]
    ii = nearest_indices(t, np.asarray(trace["time_s"], float))

    p_truth = np.column_stack([d["disp_ref_y"], d["disp_ref_x"], -d["disp_ref_z"]])
    S_truth = np.zeros_like(p_truth)
    if len(t) > 1:
        S_truth[1:] = np.cumsum(
            0.5*(p_truth[1:]+p_truth[:-1])*np.diff(t)[:,None], axis=0
        )

    E = np.zeros((len(trace), 21), float)
    theta = np.zeros(len(trace), float)
    for k,j in enumerate(ii):
        Re = zyx_matrix(d["roll_est"][j], d["pitch_est"][j], d["yaw_est"][j])
        Rr = zyx_matrix(d["roll_ref"][j], d["pitch_ref"][j], d["yaw_ref"][j])
        rv = so3_log(Re @ Rr.T)
        E[k,:3] = rv
        theta[k] = np.linalg.norm(rv)
        bg = zu_to_ned(np.array([d["gyro_bias_x"][j],d["gyro_bias_y"][j],d["gyro_bias_z"][j]]))
        ba = zu_to_ned(np.array([d["acc_bias_x"][j],d["acc_bias_y"][j],d["acc_bias_z"][j]]))
        vt = zu_to_ned(np.array([d["vel_ref_x"][j],d["vel_ref_y"][j],d["vel_ref_z"][j]]))
        at = zu_to_ned(np.array([d["acc_ref_x"][j],d["acc_ref_y"][j],d["acc_ref_z"][j]]))
        E[k,3:6] = np.array([trace["bg_x"][k],trace["bg_y"][k],trace["bg_z"][k]]) - bg
        E[k,6:9] = np.array([trace["v_x"][k],trace["v_y"][k],trace["v_z"][k]]) - vt
        E[k,9:12] = np.array([trace["p_x"][k],trace["p_y"][k],trace["p_z"][k]]) - p_truth[j]
        E[k,12:15] = np.array([trace["S_x"][k],trace["S_y"][k],trace["S_z"][k]]) - S_truth[j]
        E[k,15:18] = np.array([trace["aw_x"][k],trace["aw_y"][k],trace["aw_z"][k]]) - at
        E[k,18:21] = np.array([trace["ba_x"][k],trace["ba_y"][k],trace["ba_z"][k]]) - ba
    return E, theta


def load_exact_maps(path: Path, record: str) -> tuple[list[MapBlock], dict]:
    blocks: list[MapBlock] = []
    with path.open("rb") as f:
        magic = f.read(8)
        if magic != MAP_MAGIC:
            raise RuntimeError(f"bad exact-map magic in {path}: {magic!r}")
        version, nx, stride = struct.unpack("<III", f.read(12))
        if version != MAP_VERSION or nx != NX:
            raise RuntimeError(f"unsupported exact-map format {version}/{nx} in {path}")
        i = 0
        while True:
            prefix = f.read(MAP_PREFIX.size)
            if not prefix:
                break
            if len(prefix) != MAP_PREFIX.size:
                raise RuntimeError(f"truncated exact-map prefix in {path}")
            (t0,t1,flags,ac,mc,pc,
             tau0,sigma0,rs0,tau1,sigma1,rs1,resid) = MAP_PREFIX.unpack(prefix)
            raw = f.read(NX*NX*4)
            if len(raw) != NX*NX*4:
                raise RuntimeError(f"truncated exact-map matrix in {path}")
            phi = np.frombuffer(raw, dtype="<f4").reshape(NX,NX).astype(float)
            blocks.append(MapBlock(
                record=record,index=i,t0=t0,t1=t1,
                valid=bool(flags & (1<<0)),hybrid_jump=bool(flags & (1<<7)),
                start_live=bool(flags & (1<<1)),end_live=bool(flags & (1<<2)),
                start_active=bool(flags & (1<<3)),end_active=bool(flags & (1<<4)),
                start_lock=bool(flags & (1<<5)),end_lock=bool(flags & (1<<6)),
                start_refined=bool(flags & (1<<8)),end_refined=bool(flags & (1<<9)),
                acc_count=ac,mag_count=mc,pseudo_count=pc,
                tau0=tau0,sigma0=sigma0,rs0=rs0,
                tau1=tau1,sigma1=sigma1,rs1=rs1,
                linearization_residual=resid,phi=phi,
            ))
            i += 1
    return blocks, {"version":version,"nx":nx,"stride_samples":stride,
                    "base_period_s":stride/200.0,"block_count":len(blocks)}


def normalized_map(phi: np.ndarray, mode: str) -> np.ndarray:
    dim = 21 if mode == "A" else 18
    scale = SCALE_ACTIVE if mode == "A" else SCALE_HELD
    A = phi[:dim,:dim]
    return A * scale[np.newaxis,:] / scale[:,np.newaxis]


def compose_words(blocks: list[MapBlock], mode: str, horizon_s: float) -> list[Word]:
    if not blocks:
        return []
    base = float(np.median([b.t1-b.t0 for b in blocks if b.t1>b.t0]))
    n = max(1, int(round(horizon_s/base)))
    words: list[Word] = []
    for a in range(0, len(blocks)-n+1):
        seq = blocks[a:a+n]
        if any(not b.valid or b.hybrid_jump or b.mode != mode for b in seq):
            continue
        if any(abs(seq[k+1].t0-seq[k].t1) > 5e-3 for k in range(len(seq)-1)):
            continue
        Phi = np.eye(21)
        prefix = 1.0
        for b in seq:
            Phi = b.phi @ Phi
            prefix = max(prefix, float(np.linalg.norm(normalized_map(Phi, mode), 2)))
        words.append(Word(
            record=seq[0].record,start_index=seq[0].index,end_index=seq[-1].index,
            start_node=seq[0].start_node(),end_node=seq[-1].end_node(),
            acc_count=sum(b.acc_count for b in seq),
            mag_count=sum(b.mag_count for b in seq),
            pseudo_count=sum(b.pseudo_count for b in seq),
            phi=normalized_map(Phi,mode),prefix_norm_max=prefix,
        ))
    return words


def generalized_lambda(phi: np.ndarray, Pi: np.ndarray, Pj: np.ndarray) -> float:
    Pi = 0.5*(Pi+Pi.T); Pj = 0.5*(Pj+Pj.T)
    lam, U = np.linalg.eigh(Pi)
    if np.min(lam) <= 0:
        return math.inf
    invsqrt = U @ np.diag(1.0/np.sqrt(lam)) @ U.T
    M = invsqrt @ phi.T @ Pj @ phi @ invsqrt
    return float(np.max(np.linalg.eigvalsh(0.5*(M+M.T))))


def initial_design(words: list[Word], cap: int = 180) -> list[int]:
    by_edge: dict[tuple[SourceNode,SourceNode], list[int]] = defaultdict(list)
    scores = []
    for i,w in enumerate(words):
        by_edge[(w.start_node,w.end_node)].append(i)
        scores.append(float(np.linalg.norm(w.phi, "fro")))
    chosen: set[int] = set()
    for ids in by_edge.values():
        ids = sorted(ids, key=lambda i:scores[i], reverse=True)
        chosen.update(ids[:2])
    chosen.update(sorted(range(len(words)), key=lambda i:scores[i], reverse=True)[:cap])
    return sorted(chosen)


def solve_feasible(words: list[Word], design: list[int], dim: int,
                   rho: float) -> tuple[dict[SourceNode,np.ndarray] | None, str]:
    if cp is None:
        return None, "CVXPY_UNAVAILABLE"
    nodes = sorted({words[i].start_node for i in design} | {words[i].end_node for i in design})
    if not nodes:
        return None, "NO_NODES"
    P = {n: cp.Variable((dim,dim), symmetric=True) for n in nodes}
    I = np.eye(dim)
    constraints = [P[n] >> 2e-5*I for n in nodes]
    constraints.append(cp.sum([cp.trace(P[n]) for n in nodes]) == dim*len(nodes))
    for i in design:
        w=words[i]
        constraints.append(w.phi.T @ P[w.end_node] @ w.phi - rho*P[w.start_node] << -1e-7*I)
    prob=cp.Problem(cp.Minimize(0),constraints)
    try:
        prob.solve(solver=cp.SCS,eps=2e-5,max_iters=15000,verbose=False,
                   acceleration_lookback=10)
    except Exception as exc:
        return None, f"SOLVER_ERROR:{type(exc).__name__}"
    if prob.status not in (cp.OPTIMAL,cp.OPTIMAL_INACCURATE):
        return None, str(prob.status)
    out={}
    for n in nodes:
        if P[n].value is None:
            return None, "NO_VALUE"
        M=np.asarray(P[n].value,float)
        M=0.5*(M+M.T)
        if not np.all(np.isfinite(M)) or np.min(np.linalg.eigvalsh(M)) <= 1e-8:
            return None, "NON_SPD_VALUE"
        out[n]=M
    return out,str(prob.status)


def evaluate_metrics(words: list[Word], metrics: dict[SourceNode,np.ndarray]) -> tuple[np.ndarray,int]:
    vals=np.full(len(words),math.inf,float)
    for i,w in enumerate(words):
        if w.start_node in metrics and w.end_node in metrics:
            vals[i]=generalized_lambda(w.phi,metrics[w.start_node],metrics[w.end_node])
    return vals,int(np.argmax(vals)) if len(vals) else -1


def solve_path_metrics(words: list[Word], dim: int) -> tuple[dict,dict[SourceNode,np.ndarray] | None]:
    if len(words) < 4:
        return {"status":"INSUFFICIENT_WORDS","word_count":len(words),
                "linear_exact_replay_pass":False},None
    design=initial_design(words)
    rho_target=0.9999
    metrics=None; solver_status="not_run"; vals=None; worst_i=-1
    for iteration in range(6):
        metrics,solver_status=solve_feasible(words,design,dim,rho_target)
        if metrics is None:
            break
        vals,worst_i=evaluate_metrics(words,metrics)
        finite=vals[np.isfinite(vals)]
        worst=float(np.max(finite)) if len(finite) else math.inf
        if worst < 1.0:
            break
        order=np.argsort(vals)[::-1]
        new=[int(i) for i in order if np.isfinite(vals[i]) and i not in design][:24]
        if not new:
            break
        design=sorted(set(design)|set(new))

    pass_metrics = metrics is not None and vals is not None and np.all(np.isfinite(vals)) \
                   and float(np.max(vals)) < 1.0

    rho_feasible=rho_target if metrics is not None else None
    if not pass_metrics:
        # Find a useful noncontracting metric too.  This is diagnostic only; a
        # rho >= 1 result is never promoted to certificate success.
        for rho in (1.001,1.01,1.05,1.1,1.25,1.5,2.0,4.0,8.0,16.0):
            cand,status=solve_feasible(words,design,dim,rho)
            if cand is not None:
                metrics=cand; solver_status=status; rho_feasible=rho
                vals,worst_i=evaluate_metrics(words,metrics)
                break

    if vals is None:
        vals=np.full(len(words),math.inf)
        worst_i=0
    worst=float(np.max(vals)) if len(vals) else math.inf
    finite=vals[np.isfinite(vals)]
    result={
        "status":solver_status,
        "word_count":len(words),
        "node_count":len({w.start_node for w in words}|{w.end_node for w in words}),
        "design_word_count":len(design),
        "rho_target":rho_target,
        "rho_feasible_design":rho_feasible,
        "lambda_worst_exact_replay":worst,
        "lambda_p99_exact_replay":float(np.quantile(finite,0.99)) if len(finite) else None,
        "linear_exact_replay_pass":bool(pass_metrics),
        "worst_word":words[worst_i].summary() if words and worst_i>=0 else None,
        "worst_prefix_euclidean_gain":float(max(w.prefix_norm_max for w in words)),
    }
    return result,metrics


def choose_mode_certificate(all_blocks: dict[str,list[MapBlock]], mode: str) -> tuple[dict,dict | None]:
    dim=21 if mode=="A" else 18
    attempts=[]; best=None; best_metrics=None
    for horizon in HORIZONS_S:
        words=[]
        for blocks in all_blocks.values():
            words.extend(compose_words(blocks,mode,horizon))
        r,m=solve_path_metrics(words,dim)
        r["horizon_s"]=horizon
        attempts.append(r)
        if r.get("linear_exact_replay_pass"):
            best=r; best_metrics=m; break
        if m is not None and (best is None or r.get("lambda_worst_exact_replay",math.inf)
                              < best.get("lambda_worst_exact_replay",math.inf)):
            best=r; best_metrics=m
    if best is None:
        best={"status":"NO_WORDS","linear_exact_replay_pass":False,"horizon_s":None}
    return {"mode":mode,"selected":best,"attempts":attempts},best_metrics


def record_linear_result(blocks: list[MapBlock], mode_result: dict,
                         metrics: dict[SourceNode,np.ndarray] | None) -> dict:
    selected=mode_result["selected"]
    horizon=selected.get("horizon_s")
    if horizon is None or metrics is None:
        return {"status":"NO_METRIC","word_count":0,"lambda_worst":None}
    words=compose_words(blocks,mode_result["mode"],float(horizon))
    if not words:
        return {"status":"NOT_EXERCISED","word_count":0,"lambda_worst":None}
    vals,worst_i=evaluate_metrics(words,metrics)
    return {
        "status":"PASS" if np.all(np.isfinite(vals)) and float(np.max(vals))<1.0 else "FAIL",
        "word_count":len(words),"lambda_worst":float(np.max(vals)),
        "worst_word":words[worst_i].summary(),
    }


def handoff_hybrid(trace: np.ndarray, E: np.ndarray) -> dict:
    live=np.asarray(trace["live"],int); bias=np.asarray(trace["bias_active"],int)
    lock=np.asarray(trace["mag_lock"],int); refine=np.asarray(trace["mag_refined"],int)
    rising=lambda x: np.flatnonzero(x[1:]>x[:-1])+1
    li,bi,mi,ri=map(rising,(live,bias,lock,refine))
    h=int(li[0]) if len(li) else (0 if len(live) and live[0] else -1)
    theta=np.linalg.norm(E[:,:3],axis=1)
    return {
        "live_handoff_time_s":float(trace["time_s"][h]) if h>=0 else None,
        "handoff_theta_deg":math.degrees(float(theta[h])) if h>=0 else None,
        "bias_release_time_s":float(trace["time_s"][bi[0]]) if len(bi) else None,
        "mag_lock_time_s":float(trace["time_s"][mi[0]]) if len(mi) else None,
        "mag_refine_time_s":float(trace["time_s"][ri[0]]) if len(ri) else None,
        "bias_release_events":int(len(bi)),"mag_lock_events":int(len(mi)),
        "mag_refine_events":int(len(ri)),
    }


def stochastic_diagnostic(theta: np.ndarray, n_words: int) -> dict:
    # The raw reference test model has nine Gaussian sensor coordinates.  Keep
    # this stage explicitly diagnostic until the theorem's W-Lipschitz and
    # martingale constants are enclosed in the solved path metric.
    d=9; wstar=6.0*math.sqrt(d); lo,hi=0.0,wstar*wstar
    for _ in range(80):
        t=0.5*(lo+hi)
        q=d+2*math.sqrt(d*t)+2*t
        if q<=wstar*wstar: lo=t
        else: hi=t
    return {
        "noise_coordinate":"per-sensor-sigma-normalized pre-gate Gaussian increment",
        "gaussian_dimension":d,"w_star_normalized":wstar,"t_star":lo,
        "localization_union_bound":min(1.0,max(1,n_words)*math.exp(-lo)),
        "qualification":"WAITING_FOR_METRIC_DEPENDENT_BW_VW_ENCLOSURE",
    }


def run_record(exe: Path, data: Path, out: Path) -> tuple[Path,Path,Path,dict,bool,str]:
    trace=(out/(data.stem+"_certificate_trace.csv")).resolve()
    maps=(out/(data.stem+"_exact_maps.bin")).resolve()
    env=os.environ.copy()
    env["OU3_CERT_TRACE"]=str(trace)
    env["OU3_CERT_MAP_TRACE"]=str(maps)
    env.setdefault("OU3_CERT_TRACE_STRIDE","10")
    env.setdefault("OU3_CERT_MAP_STRIDE","50")
    env["W3D_WRITE_TIMESERIES"]="1"
    env["W3D_VALIDATION_WINDOW_SEC"]="900"
    p=subprocess.run([str(exe),"--input",str(data)],cwd=data.parent,env=env,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False)
    ts=output_csv_for(data)
    if not trace.exists() or not maps.exists() or not ts.exists():
        raise RuntimeError(f"certificate replay outputs missing for {data.name}\n{p.stdout[-6000:]}")
    return trace,maps,ts,parse_metrics(p.stdout),p.returncode==0,p.stdout


def markdown(report: dict) -> str:
    lin=report["exact_linear_source_certificate"]
    x=[
        "# OU-III exact-map numerical stability certificate","",
        f"Filter regression on eight noisy seas: **{report['filter_regression']}**",
        f"Exact executed-word linear source/path certificate: **{lin['status']}**",
        f"Numerical full source-funnel certificate: **{report['numerical_certificate']}**",
        f"Deployment theorem certificate: **{report['deployment_theorem_certificate']}**","",
        "The closed-loop matrices in this report are emitted from the estimator's actual "
        "prediction/Kalman/reset operations. They are not fitted from noisy trajectories.","",
        "## Linear path result","",
        f"- held-bias mode: horizon={lin['held']['selected'].get('horizon_s')} s, "
        f"worst lambda={lin['held']['selected'].get('lambda_worst_exact_replay')}",
        f"- active-bias mode: horizon={lin['active']['selected'].get('horizon_s')} s, "
        f"worst lambda={lin['active']['selected'].get('lambda_worst_exact_replay')}",
        f"- maximum P H^T reconstruction residual: {report['map_integrity']['max_linearization_residual']:.3e}","",
        "| Sea | RMS | exact linear | max theta | handoff theta | invalid map blocks |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in report["records"]:
        exact="PASS" if r["exact_linear_pass"] else "FAIL"
        x.append(
            f"| {r['family']} {r['Hs_m']:.2f} | {'PASS' if r['rms_regression_pass'] else 'FAIL'} | "
            f"{exact} | {r['theta_max_deg']:.2f} deg | "
            f"{r['handoff_hybrid']['handoff_theta_deg'] if r['handoff_hybrid']['handoff_theta_deg'] is not None else 'n/a'} | "
            f"{r['invalid_map_blocks']} |"
        )
    x += ["","## Claim boundary","",
          "`exact_linear_source_certificate` covers every valid ordinary-Live word executed by the eight noisy runs using solved path LMIs.",
          "Hybrid-transition blocks are reported separately rather than being forced into a fixed-dimensional ordinary-Live word.",
          "`numerical_certificate` remains BLOCKED until the exact nonlinear SO(3) word margin, handoff/invariant funnel levels, hybrid jump inequalities, and metric-dependent stochastic constants are numerically closed.",
          "`deployment_theorem_certificate` additionally requires validated continuous-source enclosure; replay coverage alone cannot establish it.",""]
    return "\n".join(x)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,default=DEFAULT_DATA_DIR)
    ap.add_argument("--output-dir",type=Path,default=DEFAULT_OUT)
    ap.add_argument("--sim",type=Path,default=TEST_DIR/"ou3-certificate-sim")
    ap.add_argument("--no-build",action="store_true")
    args=ap.parse_args()
    data_dir=args.data_dir.resolve(); out=args.output_dir.resolve(); exe=args.sim.resolve()
    out.mkdir(parents=True,exist_ok=True); (out/"logs").mkdir(exist_ok=True)
    if cp is None:
        raise RuntimeError("cvxpy is required for the exact path-LMI certificate")
    if not args.no_build:
        subprocess.run(["make","-C",str(TEST_DIR),"ou3-certificate-sim"],check=True)

    raw_records=[]; all_blocks={}; map_meta={}
    for family,hs,name in RECORDS:
        data=(data_dir/name).resolve()
        if not data.exists(): raise FileNotFoundError(f"missing reference record: {data}")
        trace_path,map_path,ts,metrics,rms_ok,log=run_record(exe,data,out)
        slug=f"{family.lower().replace('-','_')}_{hs:.2f}".replace(".","_")
        (out/"logs"/f"{slug}.log").write_text(log)
        trace=np.genfromtxt(trace_path,delimiter=",",names=True,dtype=None,encoding=None)
        E,theta=build_error_states(trace,ts)
        blocks,meta=load_exact_maps(map_path,slug)
        all_blocks[slug]=blocks; map_meta[slug]=meta
        raw_records.append((family,hs,slug,trace,E,theta,metrics,rms_ok,blocks))

    held,Ph=choose_mode_certificate(all_blocks,"H")
    active,Pa=choose_mode_certificate(all_blocks,"A")
    linear_pass=bool(held["selected"].get("linear_exact_replay_pass")
                     and active["selected"].get("linear_exact_replay_pass"))

    # Save the solved numerical metrics separately from the compact JSON report.
    arrays={}
    if Ph:
        labels=sorted(Ph); arrays["H_labels"]=np.asarray([n.label() for n in labels])
        arrays["H_P"]=np.stack([Ph[n] for n in labels])
    if Pa:
        labels=sorted(Pa); arrays["A_labels"]=np.asarray([n.label() for n in labels])
        arrays["A_P"]=np.stack([Pa[n] for n in labels])
    if arrays: np.savez_compressed(out/"path_metrics.npz",**arrays)

    records=[]
    for family,hs,slug,trace,E,theta,metrics,rms_ok,blocks in raw_records:
        Hrec=record_linear_result(blocks,held,Ph)
        Arec=record_linear_result(blocks,active,Pa)
        exercised=[x for x in (Hrec,Arec) if x["status"]!="NOT_EXERCISED"]
        rec_linear=bool(exercised) and all(x["status"]=="PASS" for x in exercised)
        invalid=sum(1 for b in blocks if not b.valid)
        records.append({
            "family":family,"Hs_m":hs,"slug":slug,
            "rms_regression_pass":bool(rms_ok),"rms_metrics":metrics,
            "theta_max_deg":math.degrees(float(np.max(theta))),
            "group_energy_max":float(max(group_energy(e[:3]) for e in E)),
            "all_attitude_inside_pi":bool(float(np.max(theta))<math.pi),
            "held_linear":Hrec,"active_linear":Arec,"exact_linear_pass":rec_linear,
            "map_blocks":len(blocks),"invalid_map_blocks":invalid,
            "hybrid_map_blocks":sum(1 for b in blocks if b.hybrid_jump),
            "max_linearization_residual":max((b.linearization_residual for b in blocks),default=math.nan),
            "handoff_hybrid":handoff_hybrid(trace,E),
            "stochastic":stochastic_diagnostic(theta,max(1,len(blocks))),
        })

    rms_all=all(r["rms_regression_pass"] for r in records)
    map_resid=max(r["max_linearization_residual"] for r in records)
    integrity=math.isfinite(map_resid) and map_resid < 5e-3
    linear_status="PASS" if linear_pass and integrity else "FAIL"
    report={
        "schema":3,
        "scope":"eight_reference_noisy_replays_default_test_seeds_exact_filter_maps",
        "record_count":len(records),
        "filter_regression":"PASS" if rms_all else "FAIL",
        "map_integrity":{"max_linearization_residual":map_resid,
                         "pass":integrity,"per_record":map_meta},
        "exact_linear_source_certificate":{
            "status":linear_status,"held":held,"active":active,
            "meaning":"solved path LMIs evaluated on every executed valid ordinary-Live word",
        },
        "numerical_certificate":"BLOCKED_AFTER_LINEAR" if linear_status=="PASS" else "FAIL_AT_LINEAR_GATE",
        "numerical_missing":[
            "nodewise exact SO(3) large-angle sector theta_star",
            "exact nonlinear source-word infimum mu_W",
            "source-shaped handoff/invariant funnel c0,b,N_H,T_H",
            "held-to-active, magnetic and tilt/cooldown jump inequalities",
            "metric-dependent non-empirical b_W,v_W stochastic constants",
        ],
        "deployment_theorem_certificate":"NOT_ESTABLISHED",
        "deployment_missing":"validated continuous-source enclosure of the numerical objects",
        "records":records,
    }
    (out/"certificate.json").write_text(json.dumps(report,indent=2,sort_keys=True))
    text=markdown(report); (out/"certificate.md").write_text(text); print(text)

    # Scientific certificate failure is reported, not hidden as a CI failure.
    # Existing filter regression and map-extraction integrity remain hard gates.
    return 1 if (not rms_all or not integrity) else 0


if __name__=="__main__":
    raise SystemExit(main())
