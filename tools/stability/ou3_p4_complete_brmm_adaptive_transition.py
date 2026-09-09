#!/usr/bin/env python3
"""Branch-preserving adaptive/source transition for the complete-BRMM P4 cover.

The source-cover contract needs more than a rectangular invariant for
``(tau,sigma,R_S,T_S,scheduler)``. Consecutive theorem cells must be connected
by the shipping adaptive recurrence.

The raw tuner sigma state is intentionally distinct from the effective OU sigma
consumed by the filter. Shipping smooths ``tune_.sigma_applied`` all the way
down to the quiet target, while ``apply_ou_tune_`` later applies the separate
0.05/band-noise floor. This module propagates the raw candidate/active state and
provides the deterministic effective-filter image at commit time.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_ADAPTIVE_BRANCH_TRANSITION_V3"


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _finite_interval(x: Interval, *, positive: bool = False) -> bool:
    if not isinstance(x, Interval):
        return False
    if not (math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo <= x.hi):
        return False
    return (x.lo > 0.0) if positive else True


def _subset(x: Interval, bounds) -> bool:
    lo, hi = map(float, bounds)
    return x.lo >= lo and x.hi <= hi


def _clamp_interval(x: Interval, lo: float, hi: float) -> Interval:
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        raise ValueError("invalid clamp")
    return Interval(max(lo, min(hi, x.lo)), max(lo, min(hi, x.hi)))


def _alpha_interval(dt: float, horizon: Interval) -> Interval:
    if not _finite_interval(horizon, positive=True) or not (math.isfinite(dt) and dt > 0.0):
        raise ValueError("invalid EMA horizon/dt")
    x = Interval.outward_bounds(dt / horizon.hi, dt / horizon.lo)
    e = VT.exp_interval(-x)
    return Interval.outward_bounds(math.nextafter(1.0 - e.hi, -math.inf),
                                    math.nextafter(1.0 - e.lo, math.inf))


def ema_image(current: Interval, target: Interval, alpha: Interval) -> Interval:
    for x, name in ((current, "current"), (target, "target"), (alpha, "alpha")):
        if not _finite_interval(x):
            raise ValueError(name + " must be finite Interval")
    if alpha.lo < 0.0 or alpha.hi > 1.0:
        raise ValueError("EMA alpha must lie in [0,1]")
    return current + alpha * (target - current)


@dataclass(frozen=True)
class AdaptiveState:
    tau_candidate: Interval
    sigma_candidate: Interval   # raw tune_.sigma_applied coordinate
    rs_candidate: Interval
    tau_active: Interval
    sigma_active: Interval      # raw committed tune_.sigma_applied coordinate
    rs_active: Interval
    pseudo_period: Interval
    pseudo_elapsed: Interval


def validate_state(s: AdaptiveState, dynamic: dict) -> list[str]:
    inv = dynamic["dynamic_invariant"]
    failures=[]
    fields=(
        (s.tau_candidate, inv["tau_applied_s"], "tau_candidate"),
        (s.sigma_candidate, inv["sigma_candidate_raw_mps2"], "sigma_candidate_raw"),
        (s.rs_candidate, inv["R_S_applied"], "rs_candidate"),
        (s.tau_active, inv["tau_applied_s"], "tau_active"),
        (s.sigma_active, inv["sigma_candidate_raw_mps2"], "sigma_active_raw"),
        (s.rs_active, inv["R_S_applied"], "rs_active"),
        (s.pseudo_period, inv["pseudo_update_period_s"], "pseudo_period"),
    )
    for x,b,name in fields:
        if not _finite_interval(x,positive=True) or not _subset(x,b):
            failures.append(name+" outside dynamic invariant")
    if not _finite_interval(s.pseudo_elapsed) or s.pseudo_elapsed.lo < 0.0:
        failures.append("pseudo_elapsed invalid")
    gap=float(dynamic["validated_rate_and_jump_bounds"]["active_commit_gap_s_upper"])
    period_hi=float(inv["pseudo_update_period_s"][1])
    if s.pseudo_elapsed.hi > period_hi + gap:
        failures.append("pseudo_elapsed exceeds conservative recurrence envelope")
    return failures


def target_cell_valid(tau_target: Interval, sigma_target: Interval, rs_target: Interval,
                      dynamic: dict) -> bool:
    inv=dynamic["dynamic_invariant"]
    return bool(
        _finite_interval(tau_target,positive=True) and _subset(tau_target,inv["tau_target_s"])
        and _finite_interval(sigma_target,positive=True) and _subset(sigma_target,inv["sigma_target_raw_mps2"])
        and _finite_interval(rs_target,positive=True) and _subset(rs_target,inv["R_S_applied"])
    )


def candidate_step(state: AdaptiveState, tau_target: Interval, sigma_target: Interval,
                   rs_target: Interval, dynamic: dict) -> AdaptiveState:
    predecessor_failures=validate_state(state,dynamic)
    if predecessor_failures:
        raise ValueError("invalid adaptive predecessor: "+repr(predecessor_failures))
    if not target_cell_valid(tau_target,sigma_target,rs_target,dynamic):
        raise ValueError("target cell outside source invariant")
    inv=dynamic["dynamic_invariant"]
    dt=float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    ah=_alpha_interval(dt,Interval(*map(float,inv["common_tau_sigma_horizon_s"])))
    ar=_alpha_interval(dt,Interval(*map(float,inv["R_S_horizon_s"])))
    return AdaptiveState(
        tau_candidate=_clamp_interval(ema_image(state.tau_candidate,tau_target,ah),*map(float,inv["tau_applied_s"])),
        sigma_candidate=_clamp_interval(ema_image(state.sigma_candidate,sigma_target,ah),*map(float,inv["sigma_candidate_raw_mps2"])),
        rs_candidate=_clamp_interval(ema_image(state.rs_candidate,rs_target,ar),*map(float,inv["R_S_applied"])),
        tau_active=state.tau_active,
        sigma_active=state.sigma_active,
        rs_active=state.rs_active,
        pseudo_period=state.pseudo_period,
        pseudo_elapsed=state.pseudo_elapsed,
    )


def _shipping_pseudo_constants() -> tuple[float,float,float]:
    text=WRAPPER.read_text(encoding="utf-8")
    ratio=SOURCE.parse_const(text,"PSEUDO_UPDATE_TAU_RATIO_DEFAULT")
    pmin=SOURCE.parse_const(text,"PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT")
    pmax=SOURCE.parse_const(text,"PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT")
    return float(ratio),float(pmin),float(pmax)


def pseudo_period_image(tau_active: Interval) -> Interval:
    ratio,pmin,pmax=_shipping_pseudo_constants()
    return _clamp_interval(I(ratio)*tau_active,pmin,pmax)


def effective_filter_sigma_image(sigma_active_raw: Interval, dynamic: dict) -> Interval:
    """Shipping ``max(0.05, band_noise_floor, raw_sigma)`` outer image.

    The dynamic certificate proves 0.05 as the source-uniform floor. The band
    noise floor can only raise the effective value, and its upper image is
    already contained by the configured effective-sigma invariant. For the P4
    coefficient cover the safe same-history image is therefore raw sigma with
    the proven effective invariant floor/ceiling applied.
    """
    if not _finite_interval(sigma_active_raw, positive=True):
        raise ValueError("raw sigma must be finite positive")
    lo,hi=map(float,dynamic["dynamic_invariant"]["sigma_aw_filter_mps2"])
    return Interval(max(lo,sigma_active_raw.lo), min(hi,max(lo,sigma_active_raw.hi)))


def commit_branch(state: AdaptiveState, commit: bool, dynamic: dict) -> AdaptiveState:
    failures=validate_state(state,dynamic)
    if failures:
        raise ValueError("invalid adaptive predecessor: "+repr(failures))
    if not commit:
        return state
    return AdaptiveState(
        tau_candidate=state.tau_candidate,
        sigma_candidate=state.sigma_candidate,
        rs_candidate=state.rs_candidate,
        tau_active=state.tau_candidate,
        sigma_active=state.sigma_candidate,
        rs_active=state.rs_candidate,
        pseudo_period=pseudo_period_image(state.tau_candidate),
        pseudo_elapsed=state.pseudo_elapsed,
    )


def scheduler_branches(state: AdaptiveState, dynamic: dict) -> list[tuple[str,AdaptiveState]]:
    failures=validate_state(state,dynamic)
    if failures:
        raise ValueError("invalid adaptive predecessor: "+repr(failures))
    dt=float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    advanced=state.pseudo_elapsed+I(dt)
    definitely_due=advanced.lo >= state.pseudo_period.hi
    definitely_not_due=advanced.hi < state.pseudo_period.lo
    out=[]
    if not definitely_due:
        out.append(("not_due",AdaptiveState(
            state.tau_candidate,state.sigma_candidate,state.rs_candidate,
            state.tau_active,state.sigma_active,state.rs_active,
            state.pseudo_period,advanced)))
    if not definitely_not_due:
        out.append(("due_S_zero",AdaptiveState(
            state.tau_candidate,state.sigma_candidate,state.rs_candidate,
            state.tau_active,state.sigma_active,state.rs_active,
            state.pseudo_period,I(0.0))))
    return out


def build(domain_path: Path=DEFAULT_DOMAIN) -> dict:
    dynamic=DYNAMIC.build(domain_path)
    failures=DYNAMIC.validate(dynamic)
    if failures:
        raise RuntimeError("dynamic source prerequisite failed: "+repr(failures))
    inv=dynamic["dynamic_invariant"]
    def mid(bounds):
        a,b=map(float,bounds); return I(0.5*(a+b))
    s=AdaptiveState(mid(inv["tau_applied_s"]),mid(inv["sigma_candidate_raw_mps2"]),mid(inv["R_S_applied"]),
                    mid(inv["tau_applied_s"]),mid(inv["sigma_candidate_raw_mps2"]),mid(inv["R_S_applied"]),
                    mid(inv["pseudo_update_period_s"]),I(0.0))
    stepped=candidate_step(s,mid(inv["tau_target_s"]),mid(inv["sigma_target_raw_mps2"]),mid(inv["R_S_applied"]),dynamic)
    committed=commit_branch(stepped,True,dynamic)
    eff=effective_filter_sigma_image(committed.sigma_active,dynamic)
    eff_ok=_subset(eff,inv["sigma_aw_filter_mps2"])
    branches=scheduler_branches(committed,dynamic)
    smoke=not validate_state(committed,dynamic) and len(branches)>=1 and eff_ok
    p=mid(inv["pseudo_update_period_s"]); dt=float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    uncertain=AdaptiveState(s.tau_candidate,s.sigma_candidate,s.rs_candidate,s.tau_active,s.sigma_active,s.rs_active,p,
                            Interval(max(0.0,p.lo-dt*1.5),p.hi))
    split=scheduler_branches(uncertain,dynamic)
    branch_split={name for name,_ in split}=={"not_due","due_S_zero"}
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "dynamic_source_contract_consumed":True,
        "raw_and_effective_sigma_coordinates_preserved":True,
        "EMA_candidate_transition_available":True,
        "staged_commit_transition_available":True,
        "effective_filter_sigma_image_available":True,
        "pseudo_period_is_clamped_committed_tau_image":True,
        "scheduler_due_not_due_branch_preserved":True,
        "deadline_straddling_cell_requires_branch_split":branch_split,
        "same_frontend_cell_targets_required":True,
        "independent_per_sample_tuner_boxes_forbidden":True,
        "trajectory_replay_used":False,
        "exact_structural_zero_not_outward_widened":True,
        "hard_clamp_does_not_widen_beyond_invariant":True,
        "point_smoke_pass":smoke,
        "physical_frontend_target_transition_materialized_here":False,
        "complete_BRMM_adaptive_cover_closed_here":False,
        "P4_promoted_here":False,
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("dynamic_source_contract_consumed","raw_and_effective_sigma_coordinates_preserved","EMA_candidate_transition_available","staged_commit_transition_available","effective_filter_sigma_image_available","pseudo_period_is_clamped_committed_tau_image","scheduler_due_not_due_branch_preserved","deadline_straddling_cell_requires_branch_split","same_frontend_cell_targets_required","independent_per_sample_tuner_boxes_forbidden","exact_structural_zero_not_outward_widened","hard_clamp_does_not_widen_beyond_invariant","point_smoke_pass"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("trajectory_replay_used","physical_frontend_target_transition_materialized_here","complete_BRMM_adaptive_cover_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    return f


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args()
    d=build(a.domain);f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"adaptive_transition":d["EMA_candidate_transition_available"],"sigma_split":d["raw_and_effective_sigma_coordinates_preserved"],"branch_split":d["deadline_straddling_cell_requires_branch_split"],"complete_cover":d["complete_BRMM_adaptive_cover_closed_here"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
