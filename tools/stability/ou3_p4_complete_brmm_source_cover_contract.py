#!/usr/bin/env python3
"""Theorem-facing same-history source-cover contract for nonlinear OU-III P4.

The nonlinear event machinery consumes one outward SAME-CELL tuple

    (state, P, H-geometry, R, true bias, committed tuner/scheduler state).

For P4, tuner coefficients are not free coordinates. The same physical signal
history must generate BOTH the WavePeriodEstimator statistic and the sigma-band
variance statistic. Only then may the dependent chain

    T -> f=1/T -> tau(f) -> T_S(tau)
    signal,f -> sigma -> R_S(tau,sigma,T_S)

feed the active EMA/commit/scheduler recurrence. A clamp-defined (f,sigma)
rectangle is retained only as a coarse sanity enclosure and cannot close the
coefficient-source obligation.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1
import ou3_p4_complete_brmm_universal_target_relation as TARGET_REL

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_SAME_HISTORY_SOURCE_COVER_CONTRACT_V3"
EVENT_KINDS = ("prediction", "aw_floor", "S_zero", "accelerometer", "magnetometer", "H_to_A")


def _shape(A) -> tuple[int, int]:
    r=len(A); c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A): raise ValueError("ragged source-cover matrix")
    return r,c


def _interval_vector(x: Sequence[Interval], n: int, name: str) -> None:
    if len(x)!=n or any(not isinstance(v,Interval) for v in x):
        raise ValueError(f"{name} must be an outward Interval vector of length {n}")


def _interval_matrix(A,r:int,c:int,name:str)->None:
    if _shape(A)!=(r,c) or any(not isinstance(v,Interval) for row in A for v in row):
        raise ValueError(f"{name} must be an outward {r}x{c} Interval matrix")


@dataclass(frozen=True)
class SourceCoverCell:
    source_token:str
    predecessor_token:str|None
    mode:str
    sample_index:int
    event_ordinal:int
    kind:str
    state:Sequence[Interval]
    P:Sequence[Sequence[Interval]]
    dt_s:Interval
    tau_applied_s:Interval
    sigma_aw_mps2:Interval
    pseudo_elapsed_s:Interval
    R:Sequence[Sequence[Interval]]|None=None
    R_provenance:str|None=None
    f_hat:Sequence[Interval]|None=None
    R_hat:Sequence[Sequence[Interval]]|None=None
    m_body:Sequence[Interval]|None=None
    true_bias:Sequence[Interval]|None=None
    bias_projection_limit:float|None=None
    radial_scale:Interval|None=None


def validate_cell(cell:SourceCoverCell)->list[str]:
    f=[]
    if not isinstance(cell.source_token,str) or not cell.source_token:f.append("missing common source-history token")
    if cell.predecessor_token is not None and not isinstance(cell.predecessor_token,str):f.append("invalid predecessor token")
    if cell.mode not in ("H","A"):return f+["mode must be H/A"]
    n=18 if cell.mode=="H" else 21
    if cell.kind not in EVENT_KINDS:f.append("unsupported event kind")
    if not isinstance(cell.sample_index,int) or cell.sample_index<0:f.append("invalid sample index")
    if not isinstance(cell.event_ordinal,int) or cell.event_ordinal<0:f.append("invalid event ordinal")
    try:_interval_vector(cell.state,n,"state");_interval_matrix(cell.P,n,n,"P")
    except ValueError as exc:f.append(str(exc))
    for name,x in (("dt_s",cell.dt_s),("tau_applied_s",cell.tau_applied_s),("sigma_aw_mps2",cell.sigma_aw_mps2),("pseudo_elapsed_s",cell.pseudo_elapsed_s)):
        if not isinstance(x,Interval) or not(math.isfinite(x.lo) and math.isfinite(x.hi)):f.append(name+" is not a finite outward interval")
    if isinstance(cell.dt_s,Interval) and cell.dt_s.lo<=0:f.append("dt_s lost positivity")
    if isinstance(cell.tau_applied_s,Interval) and cell.tau_applied_s.lo<=0:f.append("tau_applied_s lost positivity")
    if isinstance(cell.sigma_aw_mps2,Interval) and cell.sigma_aw_mps2.lo<=0:f.append("sigma_aw_mps2 lost positivity")
    if cell.radial_scale is None or not isinstance(cell.radial_scale,Interval):f.append("finite-error radial scale interval missing")
    elif cell.radial_scale.lo<0 or cell.radial_scale.hi>1:f.append("radial scale must lie inside [0,1]")
    if cell.kind in ("S_zero","accelerometer","magnetometer"):
        if cell.R is None:f.append("Joseph event missing same-cell R")
        else:
            try:_interval_matrix(cell.R,3,3,"R")
            except ValueError as exc:f.append(str(exc))
    if cell.kind=="S_zero" and cell.R_provenance!=EVENTS.ACTUAL_RS_PROVENANCE:f.append("S=0 event lacks actual applied SpectralMSE R_S provenance")
    if cell.kind=="accelerometer":
        try:_interval_vector(cell.f_hat or (),3,"f_hat");_interval_matrix(cell.R_hat or (),3,3,"R_hat")
        except ValueError as exc:f.append(str(exc))
    if cell.kind=="magnetometer":
        try:_interval_vector(cell.m_body or (),3,"m_body")
        except ValueError as exc:f.append(str(exc))
    if cell.mode=="A":
        try:_interval_vector(cell.true_bias or (),3,"true_bias")
        except ValueError as exc:f.append(str(exc))
        if not(cell.bias_projection_limit is not None and math.isfinite(float(cell.bias_projection_limit)) and float(cell.bias_projection_limit)>0):f.append("A21 cell missing deployed bias projection radius")
    return list(dict.fromkeys(f))


def joseph_event_kwargs(cell:SourceCoverCell)->dict:
    failures=validate_cell(cell)
    if failures:raise ValueError("invalid source-cover cell: "+repr(failures))
    if cell.kind not in ("S_zero","accelerometer","magnetometer"):raise ValueError("cell is not a Joseph event")
    return {"mode":cell.mode,"state":cell.state,"P":cell.P,"R":cell.R,"kind":cell.kind,"f_hat":cell.f_hat,"R_hat":cell.R_hat,"m_body":cell.m_body,"R_provenance":cell.R_provenance,"bias_true":cell.true_bias,"bias_projection_limit":cell.bias_projection_limit}


def build(domain_path:Path=DEFAULT_DOMAIN)->dict:
    path=Path(domain_path).resolve()
    source=COMPLETE.build(path);literal=WORD.build(path);entry=ENTRY.build();bias1=BIAS1.build();target_rel=TARGET_REL.build(path)
    bad={"complete_source":COMPLETE.validate(source),"literal_word":WORD.validate(literal),"hard_entry":ENTRY.validate(entry),"BIAS1":BIAS1.validate(bias1),"coarse_target_relation":TARGET_REL.validate(target_rel)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError("source-cover contract prerequisites failed: "+repr(bad))
    realization=source["BRMM_dynamic_realization"]
    return {
      "schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "cover_cell_API":"SourceCoverCell + validate_cell + joseph_event_kwargs",
      "same_history_source_token_required":True,"literal_predecessor_relation_required":True,
      "reachable_shipping_Riccati_P_required":True,"same_cell_state_P_geometry_R_tuner_scheduler_required":True,
      "actual_applied_RS_provenance_required":True,"A21_same_source_true_bias_and_projection_required":True,
      "all_radial_segments_zero_to_full_hard_entry_required":True,
      "every_valid_IMU_prediction_and_accelerometer_required":bool(literal["every_valid_imu_sample_requires_prediction"] and literal["every_valid_imu_sample_requires_accelerometer_Joseph"]),
      "all_due_S_and_asynchronous_vector_events_required":bool(literal["S_scheduler_is_executed_not_replaced_by_selected_four"] and literal["magnetometer_is_asynchronous_external_event_family"]),
      "complete_BRMM_single_history_required":bool(realization["single_history_required"]),
      "same_history_drives_translation_rotation_frontend_tuner_geometry":bool(realization["same_realization_drives_translation_rotation_frontend_tuner_geometry"]),
      "hard_entry_full_declared_scale_required":bool(entry["full_declared_scale_enforced"]),
      "BIAS1_one_common_root_driver_history_required":bool(bias1["one_root_one_parameter_history_required"]),
      "independent_P_H_R_K_boxes_forbidden":True,"independent_tuner_RS_schedule_forbidden":True,
      "independent_frequency_sigma_coordinates_forbidden":bool(target_rel["independent_frequency_sigma_coordinates_forbidden"]),
      "same_signal_history_must_generate_period_and_sigma":bool(target_rel["same_signal_history_must_generate_period_and_sigma"]),
      "wave_period_moments_must_be_retained":bool(target_rel["wave_period_moments_must_be_retained"]),
      "period_scaled_sigma_band_must_be_retained":bool(target_rel["period_scaled_sigma_band_must_be_retained"]),
      "sigma_variance_statistic_must_be_retained":bool(target_rel["sigma_variance_statistic_must_be_retained"]),
      "period_frequency_reciprocal_relation_must_be_retained":bool(target_rel["period_frequency_reciprocal_relation_must_be_retained"]),
      "coarse_frequency_sigma_rectangle_may_promote_P4":False,
      "coefficient_target_inclusion_closed":bool(target_rel["coefficient_target_inclusion_closed"]),
      "joint_estimator_relation_materialized":bool(target_rel["joint_estimator_relation_materialized_here"]),
      "tau_TS_SpectralMSE_RS_same_cell_relation_available":bool(target_rel["tau_TS_RS_are_correlated_same_cell_images"]),
      "raw_effective_sigma_split_preserved":bool(target_rel["raw_and_effective_sigma_coordinates_preserved"]),
      "trajectory_replay_or_pinned_generator_may_establish_uniform_cover":False,"finite_harmonic_source_may_replace_BRMM":False,
      "marginal_Pbar_only_cover_forbidden":True,"point_trace_schema_can_regress_interface_only":True,"point_trace_can_promote_source_uniform_cover":False,
      "source_cover_transition_operator_materialized":False,"source_cover_all_BRMM_continuations_covered":False,
      "source_cover_every_radial_segment_covered":False,"source_cover_all_Joseph_cells_correlated":False,
      "SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED":False,"P4_promoted_here":False,
      "next_obligation":"materialize one joint same-signal estimator transition: physical/private vertical signal -> WavePeriodEstimator moments -> T/f and the period-scaled AdaptiveWaveBandPass + variance/noise subtraction -> sigma; then derive tau,T_S,R_S and propagate the active schedule before attempting further P4 contraction closure",
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("same_history_source_token_required","literal_predecessor_relation_required","reachable_shipping_Riccati_P_required","same_cell_state_P_geometry_R_tuner_scheduler_required","actual_applied_RS_provenance_required","A21_same_source_true_bias_and_projection_required","all_radial_segments_zero_to_full_hard_entry_required","every_valid_IMU_prediction_and_accelerometer_required","all_due_S_and_asynchronous_vector_events_required","complete_BRMM_single_history_required","same_history_drives_translation_rotation_frontend_tuner_geometry","hard_entry_full_declared_scale_required","BIAS1_one_common_root_driver_history_required","independent_P_H_R_K_boxes_forbidden","independent_tuner_RS_schedule_forbidden","independent_frequency_sigma_coordinates_forbidden","same_signal_history_must_generate_period_and_sigma","wave_period_moments_must_be_retained","period_scaled_sigma_band_must_be_retained","sigma_variance_statistic_must_be_retained","period_frequency_reciprocal_relation_must_be_retained","tau_TS_SpectralMSE_RS_same_cell_relation_available","raw_effective_sigma_split_preserved","marginal_Pbar_only_cover_forbidden","point_trace_schema_can_regress_interface_only"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("coarse_frequency_sigma_rectangle_may_promote_P4","coefficient_target_inclusion_closed","joint_estimator_relation_materialized","trajectory_replay_or_pinned_generator_may_establish_uniform_cover","finite_harmonic_source_may_replace_BRMM","point_trace_can_promote_source_uniform_cover","source_cover_transition_operator_materialized","source_cover_all_BRMM_continuations_covered","source_cover_every_radial_segment_covered","source_cover_all_Joseph_cells_correlated","SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    return list(dict.fromkeys(f))


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"contract_ready":not f,"joint_estimator_relation":d["joint_estimator_relation_materialized"],"coefficient_target_inclusion":d["coefficient_target_inclusion_closed"],"uniform_cover_closed":d["SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED"],"next":d["next_obligation"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
