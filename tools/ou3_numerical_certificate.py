#!/usr/bin/env python3
"""Numerical source-funnel diagnostics for the deployed OU-III filter.

There are intentionally two claim levels:
  * finite_replay_certificate: check of the eight executed noisy reference runs;
  * deployment_theorem_certificate: only established by validated continuous
    source-family enclosure, never by sampled trajectories.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ou_sweep_common import PATTERNS, RECORDS

REPO = Path(__file__).resolve().parents[1]
TEST_DIR = REPO / "tests" / "kalman_ou_iii"
DEFAULT_DATA_DIR = REPO / "plots" / "kalman_ou_ii"
DEFAULT_OUT = REPO / "reports" / "results" / "ou3_numerical_certificate"
WORD_SEC = 1.0
EPS = 1e-10

# Coordinate normalization only; these do not change filter tuning.
SCALE_ACTIVE = np.array(
    [1.0]*3 + [0.02]*3 + [5.0]*3 + [10.0]*3 +
    [100.0]*3 + [2.0]*3 + [0.2]*3, dtype=float
)
SCALE_HELD = SCALE_ACTIVE[:-3]
TAU_BINS = np.array([0.0, 1.5, 2.5, 4.0, 6.0, math.inf])
SIGMA_BINS = np.array([0.0, 0.5, 1.0, 1.5, 2.5, math.inf])
RS_BINS = np.array([0.0, 0.25, 2.0, 10.0, 25.0, math.inf])


@dataclass(frozen=True)
class SourceNode:
    mode: str
    mag_lock: int
    tau_cell: int
    sigma_cell: int
    rs_cell: int

    def label(self) -> str:
        return f"{self.mode}:m{self.mag_lock}:t{self.tau_cell}:s{self.sigma_cell}:r{self.rs_cell}"


def _cell(v: float, bins: np.ndarray) -> int:
    return int(np.searchsorted(bins, v, side="right") - 1)


def source_node(row) -> SourceNode:
    return SourceNode(
        "A" if int(row["bias_active"]) else "H",
        int(row["mag_lock"]),
        _cell(float(row["tau_applied"]), TAU_BINS),
        _cell(float(row["sigma_applied"]), SIGMA_BINS),
        _cell(float(row["rs_applied"]), RS_BINS),
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
        signs = np.sign(vee)
        signs[signs == 0] = 1
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


def metric_from_samples(X: np.ndarray) -> np.ndarray:
    if len(X) < 2:
        return np.eye(X.shape[1])
    C = X.T @ X / len(X)
    ridge = max(1e-6, 1e-4*float(np.trace(C))/C.shape[0])
    return np.linalg.inv(C + ridge*np.eye(C.shape[0]))


def fit_map(X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = X.shape[1]
    phi = np.linalg.solve(X.T@X + 1e-6*np.eye(n), X.T@Y).T
    return phi, Y - X@phi.T


def generalized_lambda(phi: np.ndarray, Pi: np.ndarray, Pj: np.ndarray) -> float:
    L = np.linalg.cholesky(Pi)
    invL = np.linalg.solve(L, np.eye(L.shape[0]))
    M = invL @ phi.T @ Pj @ phi @ invL.T
    return float(np.max(np.linalg.eigvalsh(0.5*(M+M.T))))


def analyze_mode(trace: np.ndarray, E: np.ndarray, mode: str, word_rows: int) -> dict:
    dim = 21 if mode == "A" else 18
    scale = SCALE_ACTIVE if mode == "A" else SCALE_HELD
    Z = E[:,:dim]/scale
    nodes = [source_node(r) for r in trace]
    eligible = np.array([bool(int(trace["live"][k])) and nodes[k].mode == mode
                         for k in range(len(trace))])

    samples = defaultdict(list)
    for k,ok in enumerate(eligible):
        if ok:
            samples[nodes[k]].append(Z[k])
    P = {n: metric_from_samples(np.asarray(x)) for n,x in samples.items()}

    groups = defaultdict(list)
    valid_starts = set()
    for a in range(len(trace)-word_rows):
        b = a + word_rows
        if not (eligible[a] and eligible[b]):
            continue
        if any(nodes[k].mode != mode for k in range(a,b+1)):
            continue
        valid_starts.add(a)
        signature = (
            nodes[a], nodes[b],
            int(np.sum(trace["acc_accepted"][a:b])),
            int(np.sum(trace["mag_accepted"][a:b])),
            int(np.sum(trace["pseudo_due_mirror"][a:b])),
        )
        groups[signature].append((a,b))

    words, lambdas, mus, gammas = [], [], [], []
    modeled_starts = set()
    Bpre = 1.0
    for sig,pairs in groups.items():
        if len(pairs) < 2:
            continue
        ni,nj,ac,mc,pc = sig
        X = np.asarray([Z[a] for a,_ in pairs])
        Y = np.asarray([Z[b] for _,b in pairs])
        phi,resid = fit_map(X,Y)
        lam = generalized_lambda(phi,P[ni],P[nj])
        Wi = np.einsum("bi,ij,bj->b",X,P[ni],X)
        Wj = np.einsum("bi,ij,bj->b",Y,P[nj],Y)
        den = np.maximum(np.sum(X*X,axis=1),1e-12)
        mu = (Wi-Wj)/den
        gamma = float(max(0.0,np.max(Wj-min(lam,0.999999)*Wi)))
        lambdas.append(lam); mus.extend(mu.tolist()); gammas.append(gamma)
        modeled_starts.update(a for a,_ in pairs)

        for a,b in pairs:
            base = max(float(Z[a]@P[ni]@Z[a]),1e-12)
            for k in range(a+1,b):
                nk=nodes[k]
                if nk in P:
                    Bpre=max(Bpre,float(Z[k]@P[nk]@Z[k])/base)

        words.append({
            "start":ni.label(),"end":nj.label(),"samples":len(pairs),
            "acc_count":ac,"mag_count":mc,"pseudo_count":pc,
            "lambda_gen":lam,"residual_inf":float(np.max(np.abs(resid))),
            "gamma_replay":gamma,"mu_min_replay":float(np.min(mu)),
            "mu_p05_replay":float(np.quantile(mu,0.05)),
        })

    coverage = len(modeled_starts)/len(valid_starts) if valid_starts else 0.0
    worst = max(lambdas) if lambdas else None
    return {
        "mode":mode,"dimension":dim,"node_count":len(P),"word_family_count":len(words),
        "word_start_coverage":coverage,"lambda_worst":worst,"B_pre_replay":Bpre if words else None,
        "mu_min_replay":min(mus) if mus else None,
        "mu_p05_replay":float(np.quantile(mus,0.05)) if mus else None,
        "gamma_worst_replay":max(gammas) if gammas else None,
        "linear_candidate_pass":worst is not None and worst < 1.0,
        "replay_word_coverage_pass":coverage > 0.98,
        "words":words,
    }


def handoff_hybrid(trace: np.ndarray, E: np.ndarray, word_rows: int) -> dict:
    live=np.asarray(trace["live"],int); bias=np.asarray(trace["bias_active"],int)
    lock=np.asarray(trace["mag_lock"],int); refine=np.asarray(trace["mag_refined"],int)
    rising=lambda x: np.flatnonzero(x[1:]>x[:-1])+1
    li,bi,mi,ri=map(rising,(live,bias,lock,refine))
    h=int(li[0]) if len(li) else (0 if len(live) and live[0] else -1)
    z=E/SCALE_ACTIVE
    W=np.sum(z*z,axis=1)
    tail=max(0,int(0.8*len(W)))
    inner=float(np.quantile(W[tail:],0.999)) if len(W) else None
    capture=None
    if h>=0 and inner is not None:
        running_max=np.maximum.accumulate(W[::-1])[::-1]
        q=np.flatnonzero(running_max[h:] <= inner*(1+1e-9))
        capture=h+int(q[0]) if len(q) else None
    dt=float(np.median(np.diff(trace["time_s"]))) if len(trace)>1 else 0.05
    return {
        "live_handoff_time_s":float(trace["time_s"][h]) if h>=0 else None,
        "handoff_theta_deg":math.degrees(float(np.linalg.norm(E[h,:3]))) if h>=0 else None,
        "bias_release_time_s":float(trace["time_s"][bi[0]]) if len(bi) else None,
        "mag_lock_time_s":float(trace["time_s"][mi[0]]) if len(mi) else None,
        "mag_refine_time_s":float(trace["time_s"][ri[0]]) if len(ri) else None,
        "inner_replay_level":inner,
        "capture_word_count_replay":(
            int(math.ceil((capture-h)/word_rows)) if capture is not None and h>=0 else None
        ),
        "capture_time_s_replay":float((capture-h)*dt) if capture is not None and h>=0 else None,
        "hybrid_events_observed":{
            "bias_release":int(len(bi)),"mag_lock":int(len(mi)),"mag_refine":int(len(ri))
        },
    }


def stochastic_diagnostic(dW: np.ndarray, N: int) -> dict:
    # Nine raw sensor channels, normalized by their own test standard deviations.
    d=9
    wstar=6.0*math.sqrt(d)
    lo,hi=0.0,wstar*wstar
    for _ in range(80):
        t=0.5*(lo+hi)
        q=d+2*math.sqrt(d*t)+2*t
        if q<=wstar*wstar: lo=t
        else: hi=t
    localization=min(1.0,N*math.exp(-lo))
    if len(dW):
        centered=dW-np.mean(dW)
        b=float(np.max(np.maximum(centered,0.0)))
        v=float(np.var(centered))
        x=max(float(np.quantile(dW,0.999)-np.mean(dW)),0.0)
        freed=min(1.0,N*math.exp(-x*x/(2*(v+b*x/3)))) if x>0 and (v>0 or b>0) else 0.0
    else:
        b=v=freed=None
    return {
        "noise_coordinate":"per-sensor-sigma-normalized pre-gate Gaussian increment",
        "gaussian_dimension":d,"w_star_normalized":wstar,"t_star":lo,
        "localization_union_bound":localization,"b_W_replay":b,"v_W_replay":v,
        "freedman_replay_bound":freed,
        "combined_replay_bound":None if freed is None else min(1.0,localization+freed),
        "qualification":"DIAGNOSTIC_EMPIRICAL_BW_VW",
    }


DEPLOYMENT_MISSING = [
    "validated continuous-source word-family enclosure",
    "robust verified path LMIs over every source cell",
    "rigorous large-angle sector theta_star per node",
    "rigorous nonlinear source-word infimum mu_W",
    "rigorous handoff/hybrid funnel inequalities",
    "non-empirical martingale b_W and v_W bounds",
]


def analyze_record(family: str, hs: float, trace_path: Path, ts_path: Path,
                   rms_metrics: dict, rms_pass: bool) -> dict:
    trace=np.genfromtxt(trace_path,delimiter=",",names=True,dtype=None,encoding=None)
    E,theta=build_error_states(trace,ts_path)
    dt=float(np.median(np.diff(trace["time_s"]))) if len(trace)>1 else 0.05
    wr=max(1,int(round(WORD_SEC/dt)))
    H=analyze_mode(trace,E,"H",wr)
    A=analyze_mode(trace,E,"A",wr)
    modes=[m for m in (H,A) if m["word_family_count"]]
    linear=bool(modes) and all(m["linear_candidate_pass"] for m in modes)
    coverage=bool(modes) and all(m["replay_word_coverage_pass"] for m in modes)
    V=np.asarray([group_energy(x) for x in E[:,:3]])
    dW=V[wr:]-V[:-wr] if len(V)>wr else np.array([])
    theta_max=float(np.max(theta)) if len(theta) else math.nan
    replay_pass=bool(
        rms_pass and np.all(np.isfinite(E)) and theta_max<math.pi and coverage and linear
    )
    return {
        "family":family,"Hs_m":hs,"trace_rows":int(len(trace)),
        "theta_max_deg":math.degrees(theta_max),"group_energy_max":float(np.max(V)),
        "all_attitude_inside_pi":bool(theta_max<math.pi),
        "rms_regression_pass":rms_pass,"rms_metrics":rms_metrics,
        "held_mode":H,"active_mode":A,"handoff_hybrid":handoff_hybrid(trace,E,wr),
        "stochastic":stochastic_diagnostic(dW,max(1,len(V)//wr)),
        "linear_candidate_pass":linear,"finite_replay_source_coverage_pass":coverage,
        "finite_replay_certificate":"PASS" if replay_pass else "FAIL",
        "deployment_theorem_certificate":"NOT_ESTABLISHED",
        "deployment_missing":DEPLOYMENT_MISSING,
    }


def run_record(exe: Path, data: Path, out: Path) -> tuple[Path,Path,dict,bool,str]:
    trace=(out/(data.stem+"_certificate_trace.csv")).resolve()
    env=os.environ.copy()
    env["OU3_CERT_TRACE"]=str(trace)
    env.setdefault("OU3_CERT_TRACE_STRIDE","10")
    env["W3D_WRITE_TIMESERIES"]="1"
    env["W3D_VALIDATION_WINDOW_SEC"]="900"
    p=subprocess.run([str(exe),"--input",str(data)],cwd=data.parent,env=env,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False)
    ts=output_csv_for(data)
    if not trace.exists() or not ts.exists():
        raise RuntimeError(f"certificate replay outputs missing for {data.name}\n{p.stdout[-4000:]}")
    return trace,ts,parse_metrics(p.stdout),p.returncode==0,p.stdout


def markdown(report: dict) -> str:
    x=[
        "# OU-III numerical source-funnel certificate","",
        f"Overall eight-replay status: **{report['finite_replay_certificate']}**","",
        "The replay status requires RMS regression success, SO(3) chart safety, >98% "
        "modeled source-word coverage, and contraction of every fitted held/active "
        "word family. It is not a deployment theorem certificate.","",
        "| Sea | RMS | replay | max theta | linear | coverage H/A | theorem |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in report["records"]:
        x.append(
            f"| {r['family']} {r['Hs_m']:.2f} | "
            f"{'PASS' if r['rms_regression_pass'] else 'FAIL'} | {r['finite_replay_certificate']} | "
            f"{r['theta_max_deg']:.3f} deg | {'PASS' if r['linear_candidate_pass'] else 'FAIL'} | "
            f"{r['held_mode']['word_start_coverage']:.3f}/{r['active_mode']['word_start_coverage']:.3f} | "
            f"{r['deployment_theorem_certificate']} |"
        )
    x += [
        "","## Claim boundary","",
        "- `finite_replay_certificate` is a check of the eight executed noisy traces.",
        "- `linear_candidate_pass` uses fitted source-word maps and path metrics; it is not interval proof.",
        "- `deployment_theorem_certificate` remains `NOT_ESTABLISHED` until the validated continuous-source word-family enclosure, nonlinear word maps, funnel jumps and martingale constants are independently enclosed.",
        "",
    ]
    return "\n".join(x)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,default=DEFAULT_DATA_DIR)
    ap.add_argument("--output-dir",type=Path,default=DEFAULT_OUT)
    ap.add_argument("--sim",type=Path,default=TEST_DIR/"ou3-certificate-sim")
    ap.add_argument("--no-build",action="store_true")
    args=ap.parse_args()
    data_dir=args.data_dir.resolve()
    out=args.output_dir.resolve()
    exe=args.sim.resolve()
    out.mkdir(parents=True,exist_ok=True)
    (out/"logs").mkdir(exist_ok=True)
    if not args.no_build:
        subprocess.run(["make","-C",str(TEST_DIR),"ou3-certificate-sim"],check=True)

    recs=[]
    for family,hs,name in RECORDS:
        data=(data_dir/name).resolve()
        if not data.exists():
            raise FileNotFoundError(f"missing reference record: {data}")
        trace,ts,met,ok,log=run_record(exe,data,out)
        slug=f"{family.lower().replace('-','_')}_{hs:.2f}".replace(".","_")
        (out/"logs"/f"{slug}.log").write_text(log)
        recs.append(analyze_record(family,hs,trace,ts,met,ok))

    status="PASS" if all(r["finite_replay_certificate"]=="PASS" for r in recs) else "FAIL"
    report={
        "schema":2,"scope":"eight_reference_noisy_replays_default_test_seeds",
        "record_count":len(recs),"finite_replay_certificate":status,
        "deployment_theorem_certificate":"NOT_ESTABLISHED","records":recs,
    }
    (out/"certificate.json").write_text(json.dumps(report,indent=2,sort_keys=True))
    text=markdown(report)
    (out/"certificate.md").write_text(text)
    print(text)
    # Candidate certificate failure is a scientific result, not CI breakage.
    # Existing RMS regressions remain hard failures.
    return 1 if any(not r["rms_regression_pass"] for r in recs) else 0


if __name__=="__main__":
    raise SystemExit(main())
