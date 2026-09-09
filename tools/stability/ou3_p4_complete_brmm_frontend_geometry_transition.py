#!/usr/bin/env python3
"""Frontend/geometry transition contract for the complete-BRMM P4 source cover.

This module connects the existing WavePeriodEstimator source-causality lemma and
conditional vector-PE certificate to the theorem-facing P4 source cover.  It is
strict about what is and is not currently materialized.

Closed here:

* the two shipping frequency modes (fixed prior, then one-way usable-period
  takeover) and their one-sample tuner causal edge;
* the exact fixed-prior frequency cell;
* the monotone/clamped shipping frequency -> tau_target map on any supplied
  same-source frequency interval;
* the required same-history vector-geometry payload and the configured PE norm,
  separation, packet-gap and body-rate hypotheses that every vector cell must
  satisfy.

Not closed here:

* a source-uniform transition for finite WavePeriodEstimator EW moments and the
  log-period EMA over every multimodal BRMM continuation;
* the corresponding sigma_target and SpectralMSE R_S_target histories; and
* a source-uniform physical rotation/vector trajectory cell producing the
  actual f_hat/R_hat/m_body used by every Joseph event.

Accordingly this module is a transition *interface plus partial constructor*,
not a complete source cover.  A point RAO trace, a surface Tz interval, or an
independent geometry box may not fill the open coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
import math
from pathlib import Path
import re
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_wave_period_frontend as FRONTEND
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_vector_uco_certificate as VECTOR
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_FRONTEND_GEOMETRY_TRANSITION_V1"


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _finite(x: Interval, positive: bool=False) -> bool:
    return isinstance(x,Interval) and math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo <= x.hi and (not positive or x.lo > 0.0)


def _shape(A) -> tuple[int,int]:
    r=len(A); c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A): raise ValueError("ragged matrix")
    return r,c


def _member_float(text:str,name:str)->float:
    m=re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f\b",text)
    if not m: raise RuntimeError("cannot extract deployed member "+name)
    return float(m.group(1))


def _clamp_interval(x:Interval,lo:float,hi:float)->Interval:
    if not _finite(x): raise ValueError("nonfinite interval")
    if not (math.isfinite(lo) and math.isfinite(hi) and lo<=hi): raise ValueError("invalid clamp")
    return Interval.outward_bounds(max(lo,min(hi,x.lo)),max(lo,min(hi,x.hi)))


def frequency_screen(frontend:dict)->Interval:
    x=frontend["declared_inputs"]["screened_tuning_frequency_hz"]
    return Interval.outward_bounds(float(x[0]),float(x[1]))


def prior_frequency(frontend:dict)->Interval:
    x=frontend["declared_inputs"]["fixed_tuning_frequency_prior_hz"]
    return Interval.outward_bounds(float(x[0]),float(x[1]))


def tau_target_from_frequency(freq:Interval)->Interval:
    """Literal monotone/clamped target tau = tau_coeff*0.5/f."""
    if not _finite(freq,positive=True): raise ValueError("frequency must be finite positive")
    text=WRAPPER.read_text(encoding="utf-8")
    coeff=_member_float(text,"tau_coeff_")
    tmin=float(SOURCE.parse_const(text,"MIN_TAU_S")); tmax=float(SOURCE.parse_const(text,"MAX_TAU_S"))
    raw=I(0.5*coeff)/freq
    return _clamp_interval(raw,tmin,tmax)


@dataclass(frozen=True)
class FrontendModeCell:
    source_token:str
    predecessor_token:str|None
    mode:str  # fixed_prior | usable_period
    frequency_for_current_tuner:Interval
    newly_usable_after_current_tuner:bool=False


def validate_frontend_cell(cell:FrontendModeCell,frontend:dict)->list[str]:
    f=[]
    if not isinstance(cell.source_token,str) or not cell.source_token:f.append("missing source token")
    if cell.predecessor_token is not None and not isinstance(cell.predecessor_token,str):f.append("invalid predecessor token")
    if cell.mode not in ("fixed_prior","usable_period"):f.append("invalid frontend mode")
    if not _finite(cell.frequency_for_current_tuner,positive=True):f.append("invalid tuner frequency interval")
    elif cell.mode=="fixed_prior":
        p=prior_frequency(frontend)
        if cell.frequency_for_current_tuner.lo < p.lo or cell.frequency_for_current_tuner.hi > p.hi:
            f.append("fixed-prior cell detached from shipping prior")
    else:
        s=frequency_screen(frontend)
        if cell.frequency_for_current_tuner.lo < s.lo or cell.frequency_for_current_tuner.hi > s.hi:
            f.append("usable-period frequency outside shipping screen")
    return f


def takeover_successors(cell:FrontendModeCell,frontend:dict)->list[FrontendModeCell]:
    """Preserve the one-sample causal takeover edge.

    If the current sample makes the estimator newly usable, the *current* tuner
    still consumes the predecessor frequency.  The usable-period mode can begin
    only in the successor cell.  The actual estimator frequency interval for
    that successor must come from the still-open EW/log-period source transition;
    this helper therefore emits no fabricated frequency value.
    """
    if validate_frontend_cell(cell,frontend): raise ValueError("invalid frontend predecessor")
    if cell.mode=="usable_period": return [cell]
    if not cell.newly_usable_after_current_tuner:
        return [FrontendModeCell(cell.source_token,cell.source_token,"fixed_prior",prior_frequency(frontend),False)]
    # This marker represents the required edge only; a caller must attach the
    # actual same-source estimator frequency before making a theorem cell.
    return []


@dataclass(frozen=True)
class VectorGeometryCell:
    source_token:str
    predecessor_token:str|None
    f_hat:Sequence[Interval]
    R_hat:Sequence[Sequence[Interval]]
    m_body:Sequence[Interval]
    packet_gap_s:Interval
    body_rate_norm_deg_s:Interval
    sine_separation:Interval


def validate_geometry_cell(cell:VectorGeometryCell,vector:dict)->list[str]:
    f=[]
    if not isinstance(cell.source_token,str) or not cell.source_token:f.append("missing geometry source token")
    if len(cell.f_hat)!=3 or any(not _finite(x) for x in cell.f_hat):f.append("invalid f_hat interval vector")
    if len(cell.m_body)!=3 or any(not _finite(x) for x in cell.m_body):f.append("invalid m_body interval vector")
    if _shape(cell.R_hat)!=(3,3) or any(not _finite(x) for row in cell.R_hat for x in row):f.append("invalid R_hat interval matrix")
    env=vector["operating_envelope"]
    if not _finite(cell.packet_gap_s,positive=True):f.append("invalid packet gap")
    else:
        lo,hi=map(float,env["packet_gap_s"])
        if cell.packet_gap_s.lo < lo or cell.packet_gap_s.hi > hi:f.append("packet gap outside vector-PE envelope")
    if not _finite(cell.body_rate_norm_deg_s) or cell.body_rate_norm_deg_s.lo<0 or cell.body_rate_norm_deg_s.hi>float(env["body_rate_norm_upper_deg_s"]):f.append("body rate outside vector-PE envelope")
    if not _finite(cell.sine_separation) or cell.sine_separation.lo<float(env["vector_sine_separation_lower"]) or cell.sine_separation.hi>1.0:f.append("vector separation outside PE envelope")
    # Norm floors are conditions on the actual same-source vectors.  Interval
    # component boxes alone cannot prove them without a correlated norm object;
    # do not manufacture a componentwise substitute here.
    return list(dict.fromkeys(f))


def build(domain_path:Path=DEFAULT_DOMAIN)->dict:
    frontend=FRONTEND.build(REPO); ff=FRONTEND.validate(frontend)
    dynamic=DYNAMIC.build(domain_path); df=DYNAMIC.validate(dynamic)
    vector=VECTOR.build(); vf=VECTOR.validate(vector)
    bad={k:v for k,v in (("frontend",ff),("dynamic",df),("vector",vf)) if v}
    if bad: raise RuntimeError("frontend/geometry prerequisites failed: "+repr(bad))

    prior=prior_frequency(frontend)
    tau_prior=tau_target_from_frequency(prior)
    tau_screen=tau_target_from_frequency(frequency_screen(frontend))
    inv=dynamic["dynamic_invariant"]
    tlo,thi=map(float,inv["tau_target_s"])
    tau_screen_inside=tau_screen.lo>=tlo and tau_screen.hi<=thi
    prior_cell=FrontendModeCell("PRIOR_BRANCH_SMOKE",None,"fixed_prior",prior,False)
    prior_ok=not validate_frontend_cell(prior_cell,frontend)
    startup=frontend["startup_source_language"]
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "frontend_subcertificate_consumed":True,
        "dynamic_adaptive_contract_consumed":True,
        "conditional_vector_PE_contract_consumed":True,
        "fixed_prior_branch_materialized":prior_ok,
        "fixed_prior_frequency_hz":prior.as_list(),
        "fixed_prior_tau_target_s":tau_prior.as_list(),
        "screened_frequency_to_tau_target_map_materialized":tau_screen_inside,
        "screened_frequency_tau_target_image_s":tau_screen.as_list(),
        "prior_to_usable_takeover_is_one_way":bool(startup["wave_period_startup_takeover_is_one_way_latched"]),
        "tuner_consumes_previous_sample_period_state":bool(startup["tuner_consumes_previous_sample_wave_period_state"]),
        "newly_usable_period_affects_tuner_no_earlier_than_next_sample":bool(startup["first_newly_usable_wave_period_can_affect_tuner_no_earlier_than_next_valid_sample"]),
        "same_history_vector_geometry_cell_API":"VectorGeometryCell + validate_geometry_cell",
        "vector_PE_norm_separation_gap_rate_requirements_attached":True,
        "component_boxes_may_not_replace_correlated_vector_norms":True,
        "surface_Tz_may_replace_estimator_period":False,
        "independent_geometry_box_may_promote_cover":False,
        "point_RAO_trace_may_promote_cover":False,
        "finite_EW_moment_transition_materialized_here":False,
        "log_period_EMA_transition_materialized_here":False,
        "sigma_target_same_source_transition_materialized_here":False,
        "SpectralMSE_RS_target_same_source_transition_materialized_here":False,
        "physical_rotation_vector_geometry_transition_materialized_here":False,
        "complete_frontend_geometry_transition_closed_here":False,
        "P4_promoted_here":False,
        "next_obligation":(
            "materialize the finite EW moment/log-period estimator transition over the admitted multimodal BRMM continuation and a correlated physical rotation/vector-geometry transition; feed their same-source frequency/sigma/R_S targets and f_hat/R_hat/m_body cells into the adaptive/source-cover recurrences"
        ),
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("frontend_subcertificate_consumed","dynamic_adaptive_contract_consumed","conditional_vector_PE_contract_consumed","fixed_prior_branch_materialized","screened_frequency_to_tau_target_map_materialized","prior_to_usable_takeover_is_one_way","tuner_consumes_previous_sample_period_state","newly_usable_period_affects_tuner_no_earlier_than_next_sample","vector_PE_norm_separation_gap_rate_requirements_attached","component_boxes_may_not_replace_correlated_vector_norms"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("surface_Tz_may_replace_estimator_period","independent_geometry_box_may_promote_cover","point_RAO_trace_may_promote_cover","finite_EW_moment_transition_materialized_here","log_period_EMA_transition_materialized_here","sigma_target_same_source_transition_materialized_here","SpectralMSE_RS_target_same_source_transition_materialized_here","physical_rotation_vector_geometry_transition_materialized_here","complete_frontend_geometry_transition_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    for key in ("fixed_prior_tau_target_s","screened_frequency_tau_target_image_s"):
        x=d.get(key,[])
        if len(x)!=2 or not all(math.isfinite(float(v)) for v in x) or not 0<float(x[0])<=float(x[1]):f.append(key+" invalid")
    return f


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args()
    d=build(a.domain);f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"prior":d["fixed_prior_branch_materialized"],"frequency_tau":d["screened_frequency_to_tau_target_map_materialized"],"EW_transition":d["finite_EW_moment_transition_materialized_here"],"geometry_transition":d["physical_rotation_vector_geometry_transition_materialized_here"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
