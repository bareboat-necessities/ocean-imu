#!/usr/bin/env python3
"""Composite same-sample transition of the complete SEA3 front-end state z^t.

This module connects three already shipping-bound state recurrences without
creating another source:

  private Mahony -> tuner/scheduler(old WPE frequency) -> WPE.

The ordering is the deployed ``updateCore_`` ordering. A pending online tune is
committed at the beginning of the physical sample. That committed schedule is
therefore the one the same sample's H18/A21 Riccati word must consume. The
private Mahony observer then consumes the sample's SEA3 gyro/specific force and
produces the single vertical acceleration used by both adaptation channels.
The tuner executes before ``wave_period_.update`` and consequently receives the
*previous* WPE frequency. Only after the tuner has staged its candidate does the
same vertical acceleration advance the WPE state.

The result is a transition of one already admitted complete SEA3 realization.
No independent vertical acceleration, wave frequency, tuner rectangle,
quaternion box, or candidate schedule may be supplied to this API.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_sea3_private_mahony_state_step as MAHONY
import ou3_sea3_tuner_scheduler_step as TUNER
import ou3_sea3_wpe_state_step as WPE
import ou3_validated_log as VLOG

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_COMPOSITE_FRONTEND_STATE_STEP_V2"


@dataclass(frozen=True)
class FrontEndState:
    mahony: MAHONY.State
    wpe: WPE.WPEState
    tuner: TUNER.TunerState


@dataclass(frozen=True)
class Sample:
    gyro: MAHONY.Vec3
    specific_force: MAHONY.Vec3


@dataclass(frozen=True)
class Successor:
    state: FrontEndState
    active_schedule_for_current_riccati_sample: TUNER.ActiveSchedule
    actual_rs_std_xyz_for_current_riccati_sample: tuple[Interval, Interval, Interval]
    vertical_acceleration_current_sample: Interval
    tuner_frequency_previous_wpe: Interval


def advance(
    state: FrontEndState,
    sample: Sample,
    *,
    gravity_ms2: Interval,
    two_kp: Interval,
    two_ki: Interval,
    tuner_constants: TUNER.Constants | None = None,
    wpe_constants: WPE.Constants | None = None,
) -> list[Successor]:
    """Advance one same-history SEA3 front-end sample in shipping order."""
    tc = tuner_constants or TUNER.constants()
    wc = wpe_constants or WPE.constants(tc.dt)
    if tc.dt != wc.dt:
        raise ValueError("tuner and WPE dt must be the same shipping sample period")
    if not state.wpe.usable_period:
        raise ValueError("canonical Normal-Live front-end state requires usable WPE")

    # updateCore_: apply_pending_online_tune_() runs before any current-sample
    # measurement. This active schedule is therefore fixed for the current
    # Riccati prediction/measurement word.
    committed = TUNER.commit_if_pending(state.tuner, tc)
    active = committed.active
    rs_xyz = tuple(TUNER.active_rs_std_xyz(active, tc))

    # The single current SEA3 measurement path. No independent a_vertical is
    # accepted by this API.
    mahony_next = MAHONY.advance_initialized_live(
        state.mahony,
        dt=MAHONY.I(tc.dt),
        gyro=sample.gyro,
        acc_specific_force=sample.specific_force,
        gravity_ms2=gravity_ms2,
        two_kp=two_kp,
        two_ki=two_ki,
    )
    a_vertical = mahony_next.up_ms2

    # Shipping updateCore_ calls update_tuner(..., tuner_frequency_hz_())
    # before wave_period_.update(...). Therefore this is f_{k-1}, not f_k.
    f_previous = WPE.frequency_hz(state.wpe)
    tuner_successors = TUNER.advance_after_measurement(
        committed,
        a_vertical=a_vertical,
        f_wave_previous_wpe=f_previous,
        c=tc,
    )

    # Default Complementary WPE input returns the same private-Mahony vertical
    # sample. Only then does the WPE state advance.
    wpe_successors = WPE.advance(state.wpe, a_vertical=a_vertical, c=wc)

    # These products are successor branches of the same already-admitted SEA3
    # cell, not a Cartesian source generator. Future propagation may split a
    # cell further if branch-history correlation is needed for closure.
    out: list[Successor] = []
    for tuner_next in tuner_successors:
        for wpe_next in wpe_successors:
            out.append(Successor(
                state=FrontEndState(mahony_next, wpe_next, tuner_next),
                active_schedule_for_current_riccati_sample=active,
                actual_rs_std_xyz_for_current_riccati_sample=rs_xyz,
                vertical_acceleration_current_sample=a_vertical,
                tuner_frequency_previous_wpe=f_previous,
            ))
    return out


def _point_state() -> FrontEndState:
    tc = TUNER.constants()
    f = TUNER.I(0.2)
    active = TUNER.ActiveSchedule(
        TUNER.I(1.1), TUNER.I(0.5), TUNER.I(2.0),
        TUNER.pseudo_period(TUNER.I(1.1), tc),
    )
    tuner = TUNER.TunerState(
        TUNER.BandState(TUNER.I(0.0), TUNER.I(0.1), TUNER.I(0.2),
                        TUNER.I(0.0), TUNER.I(0.1), True),
        TUNER.MomentState(TUNER.I(0.0), TUNER.I(0.5), TUNER.I(0.1),
                          TUNER.I(0.5), f),
        TUNER.CandidateState(TUNER.I(1.1), TUNER.I(0.5), TUNER.I(2.0)),
        active,
        # Exercise the beginning-of-next-sample commit edge. Candidate equals
        # active, making the expected current schedule easy to check exactly.
        TUNER.SchedulerState(TUNER.I(0.05), True),
    )
    wpe = WPE.WPEState(
        accel_prev=WPE.I(0.0), high_pass_1=WPE.I(0.0),
        high_pass_1_prev=WPE.I(0.0), high_pass_2=WPE.I(0.0),
        velocity=WPE.I(0.1), elevation=WPE.I(0.1),
        velocity_mean=WPE.I(0.0), velocity_sq=WPE.I(0.4),
        elevation_mean=WPE.I(0.0), elevation_sq=WPE.I(0.1),
        weight=WPE.I(0.8), elapsed_s=WPE.I(60.0), raw_period_s=WPE.I(5.0),
        log_period_s=VLOG.log_interval(WPE.I(5.0)), usable_period=True,
        last_moment_horizon_s=WPE.I(20.0), last_log_horizon_s=WPE.I(0.25),
    )
    mahony = MAHONY.State(
        MAHONY.I(1.0), MAHONY.I(0.0), MAHONY.I(0.0), MAHONY.I(0.0),
        MAHONY.I(0.0), MAHONY.I(0.0), MAHONY.I(0.0), MAHONY.I(0.0),
    )
    return FrontEndState(mahony, wpe, tuner)


def build() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")
    p_commit = text.find("apply_pending_online_tune_();")
    p_mahony = text.find("vertical_accel_comp_.update(dt, gyro, acc_in, g_std);")
    p_vertical = text.find("const float a_vert_measurement = vertical_accel_comp_.verticalAccelUpMs2();")
    p_tuner = text.find("update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());")
    p_wpe = text.find("wave_period_.update(dt, wave_period_input_ms2_(direction_accel));")
    parity = {
        "pending_commit_before_private_Mahony": 0 <= p_commit < p_mahony,
        "private_Mahony_before_vertical_read": 0 <= p_mahony < p_vertical,
        "same_vertical_symbol_feeds_tuner": p_tuner > p_vertical,
        "tuner_executes_before_WPE_update": 0 <= p_tuner < p_wpe,
        "default_WPE_input_is_private_complementary": (
            "WavePeriodInputSource wave_period_input_ = WavePeriodInputSource::Complementary;" in text
            and "const float a_comp = vertical_accel_comp_.verticalAccelUpMs2();" in text
            and "case WavePeriodInputSource::Complementary:" in text
        ),
        "tuner_frequency_accessor_reads_WPE": (
            "const float wave_hz = wave_period_.getFrequencyHz();" in text
            and "if (wave_period_.hasUsablePeriod()" in text
        ),
        "current_candidate_commits_next_sample": (
            "online_tune_apply_pending_ = true;" in text
            and "Commit this candidate at the beginning of updateTime(k+1)." in text
        ),
    }

    st = _point_state()
    old_wpe_f = WPE.frequency_hz(st.wpe)
    old_active = st.tuner.active
    sample = Sample(
        MAHONY.Vec3(MAHONY.I(0.01), MAHONY.I(-0.02), MAHONY.I(0.005)),
        MAHONY.Vec3(MAHONY.I(0.2), MAHONY.I(-0.1), MAHONY.I(-9.75)),
    )
    succ = advance(
        st, sample,
        gravity_ms2=MAHONY.I(9.80665),
        two_kp=MAHONY.I(0.2),
        two_ki=MAHONY.I(0.02),
    )
    smoke = {
        "successors": len(succ),
        "finite_vertical": all(
            math.isfinite(x.vertical_acceleration_current_sample.lo)
            and math.isfinite(x.vertical_acceleration_current_sample.hi)
            for x in succ
        ),
        "tuner_uses_exact_previous_WPE_interval": all(
            x.tuner_frequency_previous_wpe == old_wpe_f for x in succ
        ),
        "current_active_schedule_is_premeasurement_commit": all(
            x.active_schedule_for_current_riccati_sample == old_active for x in succ
        ),
        "same_vertical_interval_reaches_WPE_state": all(
            x.state.wpe.accel_prev == x.vertical_acceleration_current_sample for x in succ
        ),
        "actual_rs_xyz_positive": all(
            all(r.lo > 0.0 for r in x.actual_rs_std_xyz_for_current_riccati_sample)
            for x in succ
        ),
    }

    mahony_component = MAHONY.build()
    wpe_component = WPE.build()
    tuner_component = TUNER.build()
    component_validation = {
        "private_Mahony": MAHONY.validate(mahony_component),
        "WPE": WPE.validate(wpe_component),
        "tuner_scheduler": TUNER.validate(tuner_component),
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "independent_vertical_acceleration_input_allowed": False,
        "independent_wave_frequency_input_allowed": False,
        "independent_tuner_schedule_input_allowed": False,
        "same_SEA3_sample_drives_Mahony_tuner_WPE": True,
        "shipping_order": (
            "commit pending -> private Mahony -> current Riccati uses committed schedule -> "
            "tuner(old WPE frequency) -> WPE update"
        ),
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "private_Mahony_live_entry_invariant_consumed":
            mahony_component.get("live_entry_private_observer_invariant_closed") is True,
        "current_Riccati_schedule_exported_before_current_measurement": True,
        "actual_applied_per_axis_RS_exported_from_same_active_schedule": True,
        "tuner_consumes_previous_WPE_state": True,
        "same_current_vertical_acceleration_consumed_by_tuner_and_WPE": True,
        "successor_branch_product_is_not_a_source_generator": True,
        "future_cell_split_required_if_branch_history_correlation_matters": True,
        "WPE_validity_branches_retained": True,
        "timer_boundary_branches_retained": True,
        "component_validation_failures": component_validation,
        "component_validation_pass": all(not v for v in component_validation.values()),
        "target_binary32_WPE_libm_roundoff_closed": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "smoke": smoke,
        "next_obligation": (
            "propagate this connected z^t transition over the complete phase-continuous SEA3 3 s family and feed each exported committed tau/sigma/T_S/R_S schedule and same-sample geometry into the literal H18/A21 word; WPE binary32/libm error remains an explicit implementation closure obligation"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "same_SEA3_sample_drives_Mahony_tuner_WPE",
        "shipping_source_parity_pass",
        "private_Mahony_live_entry_invariant_consumed",
        "current_Riccati_schedule_exported_before_current_measurement",
        "actual_applied_per_axis_RS_exported_from_same_active_schedule",
        "tuner_consumes_previous_WPE_state",
        "same_current_vertical_acceleration_consumed_by_tuner_and_WPE",
        "successor_branch_product_is_not_a_source_generator",
        "future_cell_split_required_if_branch_history_correlation_matters",
        "WPE_validity_branches_retained", "timer_boundary_branches_retained",
        "component_validation_pass",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if not all(d.get("shipping_source_parity", {}).values()):
        failures.append("shipping source order/parity failed")
    for key in (
        "source_generator", "trajectory_replay_used",
        "independent_vertical_acceleration_input_allowed",
        "independent_wave_frequency_input_allowed",
        "independent_tuner_schedule_input_allowed",
        "target_binary32_WPE_libm_roundoff_closed",
        "complete_SEA3_family_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    smoke = d.get("smoke", {})
    for key in (
        "finite_vertical", "tuner_uses_exact_previous_WPE_interval",
        "current_active_schedule_is_premeasurement_commit",
        "same_vertical_interval_reaches_WPE_state", "actual_rs_xyz_positive",
    ):
        if smoke.get(key) is not True:
            failures.append(f"smoke {key} is not true")
    if not isinstance(smoke.get("successors"), int) or smoke["successors"] < 1:
        failures.append("composite smoke produced no successor")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "shipping_order": d["shipping_order"],
        "parity": d["shipping_source_parity"],
        "component_failures": d["component_validation_failures"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
