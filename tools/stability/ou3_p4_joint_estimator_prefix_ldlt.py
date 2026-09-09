#!/usr/bin/env python3
"""Same-history joint-estimator lineage -> endpoint/every-prefix LDLT bridge.

This module is intentionally non-generative: it does not invent coefficient
boxes or master matrices. A caller must supply one ancestry-ordered sequence
whose estimator image was emitted by ``ou3_p4_joint_estimator_transition`` and
whose augmented master/sectors/source/finite-precision maps were constructed
from that same source cell. It then applies the existing full interval ISS LDLT
to every literal prefix and the endpoint.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import argparse,json
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval
import ou3_p4_joint_estimator_transition as EST
import ou3_p4_joint_iss_augmented_master as ISS
SCHEMA=3
QUALIFICATION="OU3_P4_JOINT_ESTIMATOR_EVERY_PREFIX_LDLT_BRIDGE_V3"

@dataclass(frozen=True)
class JointPrefixInput:
    source_token:str;predecessor_token:str|None;estimator_image:EST.JointImage
    master:Sequence[Sequence[Interval]];sectors:Sequence[Sequence[Sequence[Interval]]];multipliers:Sequence[float]
    source_map:Sequence[Sequence[Interval]]|None=None;gamma_s:float=0.0
    fp_map:Sequence[Sequence[Interval]]|None=None;gamma_n:float=0.0

def advance_estimator_with_lineage(state,dt,input_accel,*,child_prefix,acc_noise_floor_sigma=None,still_attenuation=None):
    if not child_prefix:raise ValueError("child prefix required")
    images=EST.step_joint(state,dt,input_accel,acc_noise_floor_sigma=acc_noise_floor_sigma,still_attenuation=still_attenuation)
    out=[]
    for i,image in enumerate(images):
        token=f"{child_prefix}:b{i}";child=replace(image.state,source_token=token,predecessor_token=state.source_token)
        out.append(replace(image,state=child))
    return out

def _shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):raise ValueError("ragged matrix")
    return r,c

def validate_prefix_sequence(cells):
    f=[]
    if not cells:return ["nonempty prefix sequence required"]
    seen=set()
    for i,c in enumerate(cells):
        if not c.source_token or c.source_token in seen:f.append(f"prefix {i}: missing/duplicate source token")
        seen.add(c.source_token)
        if i>0 and c.predecessor_token!=cells[i-1].source_token:f.append(f"prefix {i}: predecessor is not literal previous prefix")
        if i==0 and not c.predecessor_token:f.append("prefix 0: explicit root predecessor required")
        if c.estimator_image.state.source_token!=c.source_token:f.append(f"prefix {i}: estimator source token detached")
        if c.estimator_image.state.predecessor_token!=c.predecessor_token:f.append(f"prefix {i}: estimator predecessor detached")
        if c.estimator_image.status!="joint_target":f.append(f"prefix {i}: estimator did not emit joint target")
        if any(x is None for x in (c.estimator_image.frequency_hz,c.estimator_image.sigma_target_raw_mps2,c.estimator_image.tau_target_s,c.estimator_image.pseudo_period_s,c.estimator_image.rs_target)):f.append(f"prefix {i}: incomplete joint target tuple")
        n,m=_shape(c.master)
        if n==0 or n!=m:f.append(f"prefix {i}: master must be nonempty square")
        if len(c.sectors)!=len(c.multipliers) or not c.sectors:f.append(f"prefix {i}: sector/multiplier family missing")
    return list(dict.fromkeys(f))
def certify_literal_prefixes(cells):
    failures=validate_prefix_sequence(cells)
    if failures:return {"closed":False,"validation_failures":failures,"prefixes":[]}
    records=[]
    for i,c in enumerate(cells):
        ok,pivots=ISS.certify_iss(c.master,c.sectors,c.multipliers,c.source_map,c.gamma_s,c.fp_map,c.gamma_n)
        records.append({"prefix_length":i+1,"source_token":c.source_token,"predecessor_token":c.predecessor_token,"ldlt_closed":bool(ok),"pivot_lowers":pivots})
    every=all(r["ldlt_closed"] for r in records);endpoint=bool(records[-1]["ldlt_closed"])
    return {"closed":bool(every and endpoint),"validation_failures":[],"endpoint_closed":endpoint,"every_prefix_closed":every,"prefixes":records}
def build():
    s=EST.zero_state("root");children=advance_estimator_with_lineage(s,Interval.outward_bounds(.004,.006),Interval.outward_bounds(-.1,.1),child_prefix="k0")
    lineage=bool(children) and all(x.state.predecessor_token=="root" and x.state.source_token.startswith("k0:b") for x in children)
    noise=EST.deployed_acc_noise_floor_sigma()
    return {"schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD","joint_estimator_child_ancestry_materialized":lineage,"explicit_root_predecessor_retained":True,"deployed_noise_floor_inherited_from_estimator":noise.contains(.12),"literal_previous_prefix_relation_required":True,"joint_target_tuple_required_at_each_certified_prefix":True,"same_cell_augmented_master_and_ISS_maps_required":True,"endpoint_cannot_substitute_for_every_prefix":True,"full_interval_ISS_LDLT_reused":True,"independent_f_sigma_target_forbidden":True,"production_source_cover_attached_here":False,"production_endpoint_LDLT_closed_here":False,"production_every_prefix_LDLT_closed_here":False,"P4_MOTION_PASS":False,"P4_PASS":False,"P5_MAY_START":False,"next_obligation":"emit the actual admitted BRMM/radial source-cover lineage with one joint estimator image and one same-cell augmented master per literal prefix, including BIAS1 and binary32 ISS maps, then call certify_literal_prefixes on every retained lineage"}
def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("joint_estimator_child_ancestry_materialized","explicit_root_predecessor_retained","deployed_noise_floor_inherited_from_estimator","literal_previous_prefix_relation_required","joint_target_tuple_required_at_each_certified_prefix","same_cell_augmented_master_and_ISS_maps_required","endpoint_cannot_substitute_for_every_prefix","full_interval_ISS_LDLT_reused","independent_f_sigma_target_forbidden"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("production_source_cover_attached_here","production_endpoint_LDLT_closed_here","production_every_prefix_LDLT_closed_here","P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(json.dumps({"lineage":d["joint_estimator_child_ancestry_materialized"],"noise":d["deployed_noise_floor_inherited_from_estimator"],"endpoint":d["production_endpoint_LDLT_closed_here"],"every_prefix":d["production_every_prefix_LDLT_closed_here"],"P4":d["P4_PASS"],"failures":f},sort_keys=True));return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
