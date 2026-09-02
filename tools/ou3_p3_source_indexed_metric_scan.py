#!/usr/bin/env python3
"""Source-indexed translation metric diagnostic from the frozen P2 interface.

This design scan attaches a diagonal translation covariance seed to each P2
physical source node and measures seed changes across the physical projection of
legal pair-state transitions.  It is useful for sizing the eventual switched
metric, but it is not itself a P3 covariance/information propagation: the frozen
P2 contract requires the actual proof to retain pair state (c,s), or a separately
certified sufficient quotient.

No replay, domain shrink, filter change, powf/sqrtf target tightening, or theorem
gate relaxation is used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as P2
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def _node_seed(node, domain, sched):
    h=float(sched["dt_s"]); tau=Interval(*map(float,node["tau_s"])); sigma=Interval(*map(float,node["sigma_filter_committed_mps2"])); rs=Interval(*map(float,node["R_S_filter_std"]))
    Tpe=BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"],"PE recurrence")
    upper,timing=BASE.translation_upper(tau,sigma,rs,Tpe,sched); upper=list(map(float,upper))
    scales=[sigma.lo*h,sigma.lo*h*h,sigma.lo*h*h*h,sigma.lo]; scales2=[x*x for x in scales]
    norm=[BASE.up(u/BASE.down(s)) for u,s in zip(upper,scales2)]
    if any(not(math.isfinite(x) and x>0) for x in upper+norm): raise RuntimeError("metric seed lost positivity")
    return {"Sigma_translation_diagonal_upper":upper,"Sigma_over_Dh2_diagonal_upper":norm,"word_horizon_s_lower":float(timing["word_horizon_s_lower"]),"word_horizon_s_upper":float(timing["word_horizon_s_upper"])}


def _ratio(a,b): return max(BASE.up(float(x)/BASE.down(float(y))) for x,y in zip(a,b))


def _summary(v):
    x=sorted(map(float,v))
    if not x: raise RuntimeError("empty switching population")
    def q(p): return x[min(len(x)-1,max(0,int(math.ceil(p*len(x)))-1))]
    return {"count":len(x),"min":x[0],"p50":q(.5),"p90":q(.9),"p99":q(.99),"max":x[-1],"count_le_1p01":sum(z<=1.01 for z in x),"count_le_1p05":sum(z<=1.05 for z in x),"count_le_1p25":sum(z<=1.25 for z in x),"count_le_2":sum(z<=2 for z in x),"count_le_4":sum(z<=4 for z in x),"count_le_10":sum(z<=10 for z in x),"count_le_100":sum(z<=100 for z in x)}


def build(domain_path: Path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve(); domain=json.loads(path.read_text())
    if domain.get("trajectory_fit") is not False: raise RuntimeError("metric scan must not be trajectory fitted")
    rt=P2.runtime(path); nodes=rt["nodes"]; graph=rt["union_successors"]
    if len(nodes)!=NODES.EXPECTED_STATES: raise RuntimeError("P2 state count changed")
    sched=BASE.source_schedule(); seeds=[_node_seed(n,domain,sched) for n in nodes]
    covs=[]; norms=[]; jumps={"tau_index":0,"sigma_raw_index":0,"R_S_index":0}; worst_cov=worst_norm=None; edge_count=0
    for i,outs in enumerate(graph):
        for j0 in outs:
            j=int(j0); edge_count+=1
            for k in jumps: jumps[k]=max(jumps[k],abs(int(nodes[i][k])-int(nodes[j][k])))
            c=_ratio(seeds[i]["Sigma_translation_diagonal_upper"],seeds[j]["Sigma_translation_diagonal_upper"]); n=_ratio(seeds[i]["Sigma_over_Dh2_diagonal_upper"],seeds[j]["Sigma_over_Dh2_diagonal_upper"])
            covs.append(c); norms.append(n)
            if worst_cov is None or c>worst_cov["ratio"]: worst_cov={"ratio":c,"start":i,"end":j}
            if worst_norm is None or n>worst_norm["ratio"]: worst_norm={"ratio":n,"start":i,"end":j}
    spread=[]; nspread=[]
    for k in range(4):
        v=[s["Sigma_translation_diagonal_upper"][k] for s in seeds]; nv=[s["Sigma_over_Dh2_diagonal_upper"][k] for s in seeds]
        spread.append(BASE.up(max(v)/BASE.down(min(v)))); nspread.append(BASE.up(max(nv)/BASE.down(min(nv))))
    return {
        "schema":SCHEMA,"qualification":"OU3_P3_SOURCE_INDEXED_METRIC_DIAGNOSTIC_FROM_P2_CORRELATION_INTERFACE","source_only":True,"diagnostic_only":True,"trajectory_replay_used":False,"filter_changed":False,"declared_domain_changed":False,
        "P2_correlation_interface_version":P2.INTERFACE_VERSION,"P2_pair_state_interface_consumed":True,"physical_transition_projection_only_adds_edges":True,"physical_transition_projection_used_as_P3_covariance_information_proof":False,
        "translation_diagonal_seed_is_not_full_information_metric":True,"powf_sqrtf_target_tightening_used":False,"P3_PROMOTED":False,"P4_PROMOTED":False,
        "source_states":len(nodes),"finite_speed_transition_edges":edge_count,"base_untimed_transition_edges":len(nodes)**2,"finite_speed_graph_all_to_all":edge_count==len(nodes)**2,"clock_gap_samples":[min(rt["gaps"]),max(rt["gaps"])],"state_order":["v","p","S","a_w"],
        "max_partition_index_jump_per_finite_stage":jumps,"global_node_covariance_spread_by_coordinate":spread,"global_node_normalized_spread_by_coordinate":nspread,
        "directed_metric_switch_summary":{"physical_covariance_seed":_summary(covs),"dimensionless_Dh_seed":_summary(norms)},"worst_directed_metric_switch_edges":{"physical_covariance_seed":worst_cov,"dimensionless_Dh_seed":worst_norm},
        "next_obligation":"build the actual pair-state invariant covariance/information propagation on P2 segment kernels; this projected edge scan is diagnostic only","failures":[]}


def validate(d):
    f=list(d.get("failures",[]))
    if d.get("schema")!=SCHEMA:f.append("schema mismatch")
    if d.get("P2_correlation_interface_version")!=P2.INTERFACE_VERSION:f.append("P2 interface binding changed")
    for k in ("source_only","diagnostic_only","P2_pair_state_interface_consumed","physical_transition_projection_only_adds_edges","translation_diagonal_seed_is_not_full_information_metric"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in ("trajectory_replay_used","filter_changed","declared_domain_changed","physical_transition_projection_used_as_P3_covariance_information_proof","powf_sqrtf_target_tightening_used","P3_PROMOTED","P4_PROMOTED","finite_speed_graph_all_to_all"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("source_states",0))!=NODES.EXPECTED_STATES:f.append("source state count changed")
    if d.get("clock_gap_samples")!=[13,26]:f.append("clock gap changed")
    for k,row in d.get("directed_metric_switch_summary",{}).items():
        x=row.get("max")
        if not isinstance(x,(int,float)) or not(math.isfinite(float(x)) and float(x)>0):f.append(f"{k} invalid max")
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({"P2_interface":d["P2_correlation_interface_version"],"source_states":d["source_states"],"finite_speed_transition_edges":d["finite_speed_transition_edges"],"directed_metric_switch_summary":d["directed_metric_switch_summary"],"validation_failures":f},indent=2,sort_keys=True)); return 0 if not f else 2

if __name__=="__main__": raise SystemExit(main())
