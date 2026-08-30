#!/usr/bin/env python3
"""Source-dynamic reachability backend for the OU-III P4 path certificate.

The microscopic P4 frontier is driven by source cells that combine very small
applied sigma_aw with very large applied R_S.  Those combinations can occur
after a large-to-small sea transition, so simply deleting them from the source
box would be unsound.  They are, however, *transient states of the deployed
adaptation dynamics*: tau/sigma and R_S are first-order lag states driven by a
common target, and R_S_target is a deterministic function of the same
(tau_target,sigma_target) pair.

This producer constructs a conservative finite transition graph for the three
shipping adaptation states

    (tau_applied, sigma_applied, R_S_applied)

at the 0.1 s parameter-commit cadence.  Edges are admitted only when the exact
first-order EMA solution over one commit interval intersects the destination
cell for at least one source target cell.  The R_S target interval is obtained
from the deployed default SpectralMSE law with outward-rounded endpoint
arithmetic.  No replay extrema are used.

The graph is not itself a nonlinear P4 certificate.  Its purpose is to replace
an arbitrary Cartesian product of source cells by a source-reachable path
language.  In particular it detects whether the low-sigma/high-R_S corner that
dominates the old scalar P3/P4 margin can form a recurrent strongly connected
component.  If it cannot, a path Lyapunov certificate is allowed to charge that
weak cell only for its finite residence time rather than forever.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from collections import deque

import ou3_source_reachable_matrix_p3 as P3
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _member(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed member {name}")
    return float(m.group(1))


def _constexpr(text: str, name: str) -> float:
    return float(SOURCE.parse_const(text, name))


def _constants() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")
    if "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" not in text:
        raise RuntimeError("path backend requires deployed SpectralMSE law")
    return {
        "dt": _constexpr(text, "FREQ_SMOOTHER_DT"),
        "commit": _constexpr(text, "ADAPT_EVERY_SECS"),
        "tau_coeff": _member(text, "tau_coeff_"),
        "sigma_coeff": _member(text, "sigma_coeff_"),
        "adapt_tau_sea_periods": _member(text, "adapt_tau_sea_periods_"),
        "adapt_RS_mult": _member(text, "adapt_RS_mult_"),
        "min_tau": _member(text, "min_tau_s_"),
        "max_tau": _member(text, "max_tau_s_"),
        "max_sigma": _member(text, "max_sigma_a_"),
        "min_RS": _member(text, "min_R_S_"),
        "max_RS": _member(text, "max_R_S_"),
        "min_freq": _member(text, "min_tune_freq_hz_"),
        "max_freq": _member(text, "max_tune_freq_hz_"),
        "pseudo_ratio": _member(text, "pseudo_update_tau_ratio_"),
        "pseudo_min": _member(text, "pseudo_update_period_min_s_"),
        "pseudo_max": _member(text, "pseudo_update_period_max_s_"),
        "mse_coeff": _member(text, "rs_mse_coeff_"),
        "noise_density": _member(text, "rs_accel_noise_density_"),
    }


def _target_rs(tau: float, sigma: float, c: dict) -> float:
    """Deployed SpectralMSE target evaluated at one target point."""
    tau = min(max(float(tau), c["min_tau"]), c["max_tau"])
    sigma = min(max(float(sigma), 1.0e-6), c["max_sigma"])
    TS = min(max(c["pseudo_ratio"] * tau, c["pseudo_min"]), c["pseudo_max"])
    qpow = (2.0 * c["noise_density"]) ** (1.0 / 14.0)
    sigma_aB = sigma / c["sigma_coeff"]
    raw = c["mse_coeff"] * qpow * (sigma_aB ** (6.0 / 7.0)) * (tau ** (24.0 / 7.0)) / math.sqrt(TS)
    return min(max(raw, c["min_RS"]), c["max_RS"])


def _rs_interval(tau: tuple[float,float], sigma: tuple[float,float], c: dict) -> tuple[float,float]:
    # The deployed law is monotone increasing in tau and sigma on the whole
    # positive domain, including each cadence-clamp branch.
    lo = _target_rs(tau[0], sigma[0], c)
    hi = _target_rs(tau[1], sigma[1], c)
    return down(lo), up(hi)


def _ema_image(x: tuple[float,float], target: tuple[float,float], horizon: tuple[float,float], dt: float) -> tuple[float,float]:
    if not (horizon[0] > 0.0 and horizon[1] >= horizon[0]):
        raise RuntimeError("invalid EMA horizon")
    # x+ = a*x + (1-a)*u, a=exp(-dt/T).  For positive T, a increases with T.
    alo = down(math.exp(-dt / horizon[0]))
    ahi = up(math.exp(-dt / horizon[1]))
    vals = []
    for a in (alo, ahi):
        for xx in x:
            for uu in target:
                vals.append(a * xx + (1.0 - a) * uu)
    return down(min(vals)), up(max(vals))


def _overlap(a: tuple[float,float], b: tuple[float,float]) -> bool:
    return not (a[1] < b[0] or b[1] < a[0])


def _edges(lo: float, hi: float, n: int) -> list[float]:
    return P3.geom_edges(float(lo), float(hi), int(n))


def _cells(edges: list[float]) -> list[tuple[float,float]]:
    return [(down(edges[i]), up(edges[i+1])) for i in range(len(edges)-1)]


def _tar_tau_from_freq(freq: tuple[float,float], c: dict) -> tuple[float,float]:
    lo = c["tau_coeff"] * 0.5 / freq[1]
    hi = c["tau_coeff"] * 0.5 / freq[0]
    return down(max(c["min_tau"], min(c["max_tau"], lo))), up(max(c["min_tau"], min(c["max_tau"], hi)))


def _scc(graph: list[list[int]]) -> list[list[int]]:
    n=len(graph); order=[]; seen=[False]*n
    def dfs(v):
        seen[v]=True
        for w in graph[v]:
            if not seen[w]: dfs(w)
        order.append(v)
    for v in range(n):
        if not seen[v]: dfs(v)
    rg=[[] for _ in range(n)]
    for v,ws in enumerate(graph):
        for w in ws: rg[w].append(v)
    comps=[]; seen=[False]*n
    def rdfs(v,acc):
        seen[v]=True; acc.append(v)
        for w in rg[v]:
            if not seen[w]: rdfs(w,acc)
    for v in reversed(order):
        if not seen[v]:
            acc=[]; rdfs(v,acc); comps.append(acc)
    return comps


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain=json.loads(Path(domain_path).read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("path reachability domain must not be trajectory fitted")
    c=_constants()

    # Keep the graph intentionally moderate for CI.  These are proof cells, not
    # replay bins.  Geometric partitions retain the low-end resolution that
    # matters for the pathological corner.
    tau_cells=_cells(_edges(max(c["min_tau"],c["tau_coeff"]*0.5/c["max_freq"]), c["max_tau"], 10))
    sig_cells=_cells(_edges(0.05, c["max_sigma"], 8))
    rs_cells=_cells(_edges(c["min_RS"], c["max_RS"], 10))
    freq_cells=_cells(_edges(c["min_freq"], c["max_freq"], 8))

    states=[]; index={}
    for it,t in enumerate(tau_cells):
        for isg,s in enumerate(sig_cells):
            for ir,r in enumerate(rs_cells):
                index[(it,isg,ir)]=len(states)
                states.append({"tau":t,"sigma":s,"RS":r,"ijk":[it,isg,ir]})

    graph=[set() for _ in states]
    dt=c["commit"]
    # Targets are exogenous measurement-derived quantities, but coupled through
    # the deployed law.  Enumerate target frequency and sigma cells, compute the
    # induced tau and R_S targets, then propagate all three lag states.
    target_boxes=[]
    for f in freq_cells:
        tt=_tar_tau_from_freq(f,c)
        for ss in sig_cells:
            rr=_rs_interval(tt,ss,c)
            # tau/sigma common horizon = 0.4*Tsea = 0.4*tau_target/tau_coeff.
            ht=(down(c["adapt_tau_sea_periods"]*tt[0]/c["tau_coeff"]), up(c["adapt_tau_sea_periods"]*tt[1]/c["tau_coeff"]))
            # R_S horizon = 1.5*tau_target (slew-log is deployed zero).
            hr=(down(c["adapt_RS_mult"]*tt[0]), up(c["adapt_RS_mult"]*tt[1]))
            target_boxes.append((tt,ss,rr,ht,hr))

    for i,st in enumerate(states):
        candidate=set()
        for tt,ss,rr,ht,hr in target_boxes:
            ti=_ema_image(st["tau"],tt,ht,dt)
            si=_ema_image(st["sigma"],ss,ht,dt)
            ri=_ema_image(st["RS"],rr,hr,dt)
            for jt,tcell in enumerate(tau_cells):
                if not _overlap(ti,tcell): continue
                for js,scell in enumerate(sig_cells):
                    if not _overlap(si,scell): continue
                    for jr,rcell in enumerate(rs_cells):
                        if _overlap(ri,rcell): candidate.add(index[(jt,js,jr)])
        graph[i]=candidate

    gl=[sorted(x) for x in graph]
    comps=_scc(gl)
    comp_of={v:k for k,cc in enumerate(comps) for v in cc}
    recurrent=set()
    for k,cc in enumerate(comps):
        if len(cc)>1 or (len(cc)==1 and cc[0] in graph[cc[0]]): recurrent.update(cc)

    # Exact old worst-cell bands from CI diagnostic, expressed as intersection
    # predicates.  This is a theorem-source region, not a replay-derived point.
    bad=[]
    for i,st in enumerate(states):
        xlo=c["dt"]/st["tau"][1]; xhi=c["dt"]/st["tau"][0]
        if _overlap(st["sigma"],(0.05,0.13025855423486765)) and _overlap(st["RS"],(149.21548743644342,400.0)) and _overlap((xlo,xhi),(0.00041666665735344083,0.0004837652693428343)):
            bad.append(i)
    bad_recurrent=[i for i in bad if i in recurrent]
    bad_sccs=sorted(set(comp_of[i] for i in bad))

    # Longest residence entirely inside the bad induced subgraph, unless it has
    # a cycle.  A cycle means no finite residence claim may be made yet.
    badset=set(bad)
    bad_graph={v:[w for w in gl[v] if w in badset] for v in bad}
    cycle=any(v in recurrent for v in bad)
    max_steps=None
    if not cycle:
        memo={}
        def longest(v):
            if v in memo: return memo[v]
            memo[v]=1+max((longest(w) for w in bad_graph[v]),default=0)
            return memo[v]
        max_steps=max((longest(v) for v in bad),default=0)

    return {
        "qualification":"OU3_P4_SOURCE_DYNAMIC_PATH_REACHABILITY",
        "source_only":True,
        "trajectory_replay_used":False,
        "deployed_default_law":"SpectralMSE",
        "commit_period_s":dt,
        "partition":{"tau":len(tau_cells),"sigma":len(sig_cells),"R_S":len(rs_cells),"states":len(states),"target_boxes":len(target_boxes)},
        "transition_edges":sum(len(x) for x in gl),
        "strongly_connected_components":len(comps),
        "recurrent_states":len(recurrent),
        "old_worst_corner_state_count":len(bad),
        "old_worst_corner_recurrent_state_count":len(bad_recurrent),
        "old_worst_corner_has_recurrent_cycle":bool(cycle),
        "old_worst_corner_max_consecutive_commit_steps_upper":max_steps,
        "old_worst_corner_max_residence_s_upper":None if max_steps is None else up(max_steps*dt),
        "old_worst_corner_scc_count":len(bad_sccs),
        "path_graph_ready":True,
        "usable_P4_promoted":False,
        "next_obligation":"propagate complete-word Phi/Omega and exact nonlinear return map on this source-reachable graph; weak transient cells must be charged by finite residence rather than a global min(delta)",
        "failures":[],
    }


def validate(d: dict) -> list[str]:
    f=list(d.get("failures",[]))
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False: f.append("path graph is not source-only")
    if d.get("path_graph_ready") is not True: f.append("path graph not ready")
    if not int(d.get("partition",{}).get("states",0))>0: f.append("empty path state partition")
    if not int(d.get("transition_edges",0))>0: f.append("empty path transition relation")
    if d.get("usable_P4_promoted") is not False: f.append("reachability stage prematurely promoted P4")
    return f


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(a.domain.resolve()); f=validate(d); d["validation_failures"]=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({k:d[k] for k in ("partition","transition_edges","strongly_connected_components","recurrent_states","old_worst_corner_state_count","old_worst_corner_recurrent_state_count","old_worst_corner_has_recurrent_cycle","old_worst_corner_max_residence_s_upper","failures")},indent=2,sort_keys=True))
    return 0 if not f else 2

if __name__=="__main__": raise SystemExit(main())
