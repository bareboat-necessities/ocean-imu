#!/usr/bin/env python3
"""Trusted typed execution kernel for one complete SEA3 Normal-Live word.

This module is deliberately downstream of the canonical hard-window provider.
It creates no sea/source family.  Its inputs are already provider-certified
same-history sample coordinates and the provider-certified source-reachable
word-start covariances.  The canonical JSON boundary remains gated in
``ou3_sea3_complete_window_executor`` until SEA0 closes.

For every valid IMU sample the kernel preserves shipping order:

  committed schedule -> H/A prediction -> requested covariance-dependent
  a_w floor -> every due S=0 update -> accelerometer Joseph -> exact
  measurement-only Mahony/tuner/WPE transition -> asynchronous magnetic
  Joseph events occurring after that IMU updateCore_ call.

The private Mahony consumes the raw gyro measurement.  The MEKF prediction
consumes the bias-corrected body rate.  They are distinct coordinates of the
same provider transition witness; this module never equates them or introduces
an independent gyro-bias bound.

Front-end interval branch splits are retained.  No favorable successor is
selected.  Each successor receives a deep copy of the same post-magnetometer H
and A Riccati states and continues as a separate source cell.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, IntervalMatrix
import ou3_full_process_ucc as PROCESS
import ou3_sea3_aw_covariance_floor as FLOOR
import ou3_sea3_frontend_state_step as FRONTEND
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_shipping_prediction_primitives as PRED
import ou3_sea3_tuner_scheduler_step as TUNER
import ou3_sea3_private_mahony_state_step as MAHONY

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_TRUSTED_TYPED_COMPLETE_WINDOW_EXECUTION_KERNEL_V2"


@dataclass(frozen=True)
class MagneticEvent:
    m_body: tuple[Interval, Interval, Interval]


@dataclass(frozen=True)
class SampleCoordinates:
    """One provider-certified source transition payload.

    ``gyro_measurement`` feeds only the private measurement-only Mahony.
    ``omega_body_corrected`` feeds the shipping MEKF F/Q attitude block.
    Both must descend from the same SEA3/source transition witness upstream.
    """
    gyro_measurement: MAHONY.Vec3
    omega_body_corrected: tuple[Interval, Interval, Interval]
    specific_force: MAHONY.Vec3
    f_cog_body: tuple[Interval, Interval, Interval]
    R_wb: IntervalMatrix
    due_S: bool
    aw_floor_requested: bool
    magnetometer_events_after_imu: tuple[MagneticEvent, ...] = ()


@dataclass
class ExecutionBranch:
    frontend: FRONTEND.FrontEndState
    H: WORD.LiteralWordState
    A: WORD.LiteralWordState
    source_cell_id: str


@dataclass(frozen=True)
class KernelConstants:
    h: Interval
    gyro_variance_density_xyz: tuple[Interval, Interval, Interval]
    gyro_bias_variance_density: Interval
    accel_bias_tau_s: Interval
    accel_bias_process_variance_density: Interval
    Racc: IntervalMatrix
    Rmag: IntervalMatrix
    gravity: Interval
    two_kp: Interval
    two_ki: Interval
    tuner: TUNER.Constants


def _process_constants(domain_path: Path = DEFAULT_DOMAIN) -> KernelConstants:
    process = PROCESS.build()
    pf = PROCESS.validate(process)
    if pf:
        raise RuntimeError(f"shipping process constants invalid: {pf}")
    word = WORD.build(Path(domain_path).resolve())
    wf = WORD.validate(word)
    if wf:
        raise RuntimeError(f"literal word prerequisite invalid: {wf}")
    live = LIVE.build(Path(domain_path).resolve())
    lf = LIVE.validate(live)
    if lf:
        raise RuntimeError(f"Live seed prerequisite invalid: {lf}")
    if float(live["aw_live_seed"]["S_factor"]) != 1.0:
        raise RuntimeError("typed kernel currently requires source-certified shipping S_factor=1")

    c = process["source_constants"]
    h = Interval(*process["configured_runtime"]["imu_dt_outward_interval_s"])
    qg = tuple(PRED.I(float(x) * float(x)) for x in c["gyro_noise_density_rad_sqrt_s_per_axis"])
    qb = PRED.I(float(c["gyro_bias_rw_variance_density"]))
    tau_ba = PRED.I(float(c["accel_bias_tau_s"]))
    q_ba = PRED.I(float(c["accel_bias_process_variance_density"]))
    runtime = word["measurement_runtime"]
    Racc = WORD.diagonal_R(runtime["accelerometer_std_mps2"])
    Rmag = WORD.diagonal_R(runtime["magnetometer_std_uT"])
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    gravity = PRED.I(float(domain["startup"]["gravity_mps2"]))
    return KernelConstants(
        h=h,
        gyro_variance_density_xyz=qg,
        gyro_bias_variance_density=qb,
        accel_bias_tau_s=tau_ba,
        accel_bias_process_variance_density=q_ba,
        Racc=Racc,
        Rmag=Rmag,
        gravity=gravity,
        two_kp=MAHONY.I(0.2),
        two_ki=MAHONY.I(0.02),
        tuner=TUNER.constants(),
    )


def _prediction(
    mode: str,
    sample: SampleCoordinates,
    active: TUNER.ActiveSchedule,
    c: KernelConstants,
) -> tuple[IntervalMatrix, IntervalMatrix]:
    Faa, Qaa = PRED.attitude_gyro_bias_F_Q(
        sample.omega_body_corrected,
        c.h,
        c.gyro_variance_density_xyz,
        c.gyro_bias_variance_density,
    )
    std_xyz = (active.sigma, active.sigma, active.sigma)
    Fll, Qll = PRED.translation_F_Q(active.tau, c.h, std_xyz)
    if mode == "A":
        phi_ba, Qba = PRED.active_accel_bias_F_Q(
            c.h,
            c.accel_bias_tau_s,
            c.accel_bias_process_variance_density,
        )
        return WORD.pack_prediction(
            mode, Faa, Qaa, Fll, Qll, phi_ba=phi_ba, Q_ba=Qba
        )
    return WORD.pack_prediction(mode, Faa, Qaa, Fll, Qll)


def _apply_one_mode_imu_core(
    word: WORD.LiteralWordState,
    sample: SampleCoordinates,
    active: TUNER.ActiveSchedule,
    rs_xyz: Sequence[Interval],
    c: KernelConstants,
) -> str | None:
    """Apply prediction/floor/S/accel; async mag is applied after front-end."""
    F, Q = _prediction(word.mode, sample, active, c)
    # Prediction is executed here rather than through WORD.apply_imu_sample
    # because the shipping covariance floor depends on the just-predicted P.
    WORD.BACKEND.predict(word.riccati, F, Q)
    word.event_log.append("prediction")

    floor_case: str | None = None
    if sample.aw_floor_requested:
        target = FLOOR.stationary_sigma_isotropic(active.sigma)
        Paw = FLOOR.aw_block(word.riccati.P, WORD.OFF_AW)
        Delta, floor_case = FLOOR.positive_part_enclosure(target, Paw)
        WORD.BACKEND.add_psd_floor(word.riccati, WORD.aw_floor_increment(word.mode, Delta))
        word.aw_floor_applications += 1
        word.event_log.append("aw_floor")

    if sample.due_S:
        WORD.BACKEND.joseph_measurement(
            word.riccati,
            WORD.H_S_zero(word.mode),
            WORD.R_S_zero(rs_xyz),
        )
        word.S_updates += 1
        word.event_log.append("S_zero")

    WORD.BACKEND.joseph_measurement(
        word.riccati,
        WORD.H_accelerometer(word.mode, sample.f_cog_body, sample.R_wb),
        c.Racc,
    )
    word.accel_updates += 1
    word.imu_samples += 1
    word.event_log.append("accelerometer")
    return floor_case


def advance_branch(
    branch: ExecutionBranch,
    sample: SampleCoordinates,
    *,
    constants: KernelConstants,
    next_cell_prefix: str,
) -> tuple[list[ExecutionBranch], dict]:
    """Advance one connected source cell without selecting a front-end branch."""
    committed = TUNER.commit_if_pending(branch.frontend.tuner, constants.tuner)
    active = committed.active
    rs_xyz = tuple(TUNER.active_rs_std_xyz(active, constants.tuner))

    H_post = copy.deepcopy(branch.H)
    A_post = copy.deepcopy(branch.A)
    floor_H = _apply_one_mode_imu_core(H_post, sample, active, rs_xyz, constants)
    floor_A = _apply_one_mode_imu_core(A_post, sample, active, rs_xyz, constants)

    # Shipping updateCore_ completes the private measurement-only
    # Mahony/tuner/WPE transition after the accelerometer correction.  Magnetic
    # updates are external calls, so events declared "after_imu" are applied
    # only after this front-end successor family has been formed.
    frontend_successors = FRONTEND.advance(
        branch.frontend,
        FRONTEND.Sample(sample.gyro_measurement, sample.specific_force),
        gravity_ms2=constants.gravity,
        two_kp=constants.two_kp,
        two_ki=constants.two_ki,
        tuner_constants=constants.tuner,
    )
    if not frontend_successors:
        raise RuntimeError("complete front-end transition produced no successor")
    for succ in frontend_successors:
        if succ.active_schedule_for_current_riccati_sample != active:
            raise RuntimeError("front-end current active schedule disagrees with Riccati schedule")
        if tuple(succ.actual_rs_std_xyz_for_current_riccati_sample) != rs_xyz:
            raise RuntimeError("front-end actual R_S disagrees with Riccati R_S")

    for event in sample.magnetometer_events_after_imu:
        WORD.apply_magnetometer(H_post, m_body=event.m_body, Rmag=constants.Rmag)
        WORD.apply_magnetometer(A_post, m_body=event.m_body, Rmag=constants.Rmag)

    out: list[ExecutionBranch] = []
    for i, succ in enumerate(frontend_successors):
        out.append(ExecutionBranch(
            frontend=succ.state,
            H=copy.deepcopy(H_post),
            A=copy.deepcopy(A_post),
            source_cell_id=f"{next_cell_prefix}:{i}",
        ))
    return out, {
        "frontend_successors": len(frontend_successors),
        "H_floor_case": floor_H,
        "A_floor_case": floor_A,
        "same_active_schedule_verified": True,
        "same_actual_RS_verified": True,
        "frontend_completed_before_async_mag": True,
    }


def execute_typed_window(
    *,
    frontend_entry: FRONTEND.FrontEndState,
    P0_H: IntervalMatrix,
    P0_A: IntervalMatrix,
    samples: Sequence[SampleCoordinates],
    domain_path: Path = DEFAULT_DOMAIN,
    branch_limit: int = 100000,
) -> tuple[list[ExecutionBranch], dict]:
    """Execute already-provider-certified typed coordinates through both modes.

    This function intentionally does not call the SEA0 provider gate; the
    canonical JSON executor does that before deserialization.  Direct callers
    are test/development kernels and cannot set a P3 promotion flag.
    """
    c = _process_constants(domain_path)
    branches = [ExecutionBranch(
        frontend=frontend_entry,
        H=WORD.initialize_word("H", P0_H),
        A=WORD.initialize_word("A", P0_A),
        source_cell_id="root",
    )]
    floor_cases: dict[str, int] = {}
    max_branches = 1
    all_frontend_before_mag = True
    for k, sample in enumerate(samples):
        next_branches: list[ExecutionBranch] = []
        for j, branch in enumerate(branches):
            successors, meta = advance_branch(
                branch,
                sample,
                constants=c,
                next_cell_prefix=f"k{k}:parent{j}",
            )
            next_branches.extend(successors)
            all_frontend_before_mag = all_frontend_before_mag and bool(
                meta["frontend_completed_before_async_mag"]
            )
            for key in ("H_floor_case", "A_floor_case"):
                case = meta[key]
                if case is not None:
                    floor_cases[case] = floor_cases.get(case, 0) + 1
        if len(next_branches) > branch_limit:
            raise RuntimeError(
                f"front-end branch count {len(next_branches)} exceeds validated execution limit {branch_limit}; "
                "provider must partition the source cell more tightly, not select a successor"
            )
        branches = next_branches
        max_branches = max(max_branches, len(branches))
    return branches, {
        "samples_executed": len(samples),
        "endpoint_branches": len(branches),
        "max_branch_count": max_branches,
        "floor_cases": floor_cases,
        "same_word_executed_H18_A21": True,
        "frontend_completed_before_async_mag": all_frontend_before_mag,
        "favorable_frontend_successor_selected": False,
        "kernel_self_test_only_not_P3": True,
    }


def _diag_P(n: int, value: float) -> IntervalMatrix:
    z = PRED.I(0.0)
    d = PRED.I(value)
    return [[d if i == j else z for j in range(n)] for i in range(n)]


def _smoke(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    frontend = FRONTEND._point_state()  # retained point fixture; never theorem input
    z = MAHONY.I(0.0)
    sample = SampleCoordinates(
        gyro_measurement=MAHONY.Vec3(MAHONY.I(0.01), MAHONY.I(-0.02), MAHONY.I(0.005)),
        omega_body_corrected=(PRED.I(0.01), PRED.I(-0.02), PRED.I(0.005)),
        specific_force=MAHONY.Vec3(MAHONY.I(0.2), MAHONY.I(-0.1), MAHONY.I(-9.75)),
        f_cog_body=(PRED.I(0.0), PRED.I(0.0), PRED.I(-9.80665)),
        R_wb=[
            [PRED.I(1.0), z, z],
            [z, PRED.I(1.0), z],
            [z, z, PRED.I(1.0)],
        ],
        due_S=True,
        aw_floor_requested=True,
        magnetometer_events_after_imu=(
            MagneticEvent((PRED.I(20.0), PRED.I(0.0), PRED.I(40.0))),
        ),
    )
    branches, meta = execute_typed_window(
        frontend_entry=frontend,
        P0_H=_diag_P(18, 2.0),
        P0_A=_diag_P(21, 2.0),
        samples=[sample],
        domain_path=domain_path,
        branch_limit=32,
    )
    return {
        **meta,
        "all_endpoint_H_events_present": all(
            b.H.imu_samples == 1 and b.H.accel_updates == 1 and b.H.S_updates == 1
            and b.H.mag_updates == 1 and b.H.aw_floor_applications == 1
            for b in branches
        ),
        "all_endpoint_A_events_present": all(
            b.A.imu_samples == 1 and b.A.accel_updates == 1 and b.A.S_updates == 1
            and b.A.mag_updates == 1 and b.A.aw_floor_applications == 1
            for b in branches
        ),
        "all_endpoint_event_logs_have_mag_after_accel": all(
            b.H.event_log.index("magnetometer") > b.H.event_log.index("accelerometer")
            and b.A.event_log.index("magnetometer") > b.A.event_log.index("accelerometer")
            for b in branches
        ),
        "decomposition_identity_H_enclosed": all(
            WORD.BACKEND.decomposition_identity_enclosed(b.H.riccati) for b in branches
        ),
        "decomposition_identity_A_enclosed": all(
            WORD.BACKEND.decomposition_identity_enclosed(b.A.riccati) for b in branches
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    floor = FLOOR.build()
    ff = FLOOR.validate(floor)
    if ff:
        raise RuntimeError(f"a_w floor prerequisite invalid: {ff}")
    smoke = _smoke(domain_path)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates": True,
        "current_active_schedule_derived_before_current_Riccati_sample": True,
        "current_actual_per_axis_RS_derived_from_same_committed_schedule": True,
        "covariance_floor_request_not_increment_is_source_event": True,
        "covariance_floor_increment_computed_from_current_mode_P": True,
        "H18_A21_floor_increments_not_forced_equal": True,
        "front_end_branch_splits_retained": True,
        "frontend_completed_before_async_mag": True,
        "favorable_front_end_branch_selection_allowed": False,
        "typed_execution_kernel_ready": True,
        "canonical_provider_gate_bypassed_by_this_status": False,
        "smoke": smoke,
        "P3_promoted": False,
        "next_obligation": (
            "deserialize only a canonical provider-certified 601-sample artifact into these typed coordinates, "
            "execute all retained front-end branches, and certify every endpoint H18/A21 full-matrix LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates",
        "current_active_schedule_derived_before_current_Riccati_sample",
        "current_actual_per_axis_RS_derived_from_same_committed_schedule",
        "covariance_floor_request_not_increment_is_source_event",
        "covariance_floor_increment_computed_from_current_mode_P",
        "H18_A21_floor_increments_not_forced_equal",
        "front_end_branch_splits_retained",
        "frontend_completed_before_async_mag",
        "typed_execution_kernel_ready",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator", "trajectory_replay_used",
        "favorable_front_end_branch_selection_allowed",
        "canonical_provider_gate_bypassed_by_this_status", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    s = d.get("smoke", {})
    for key in (
        "same_word_executed_H18_A21",
        "frontend_completed_before_async_mag",
        "all_endpoint_H_events_present",
        "all_endpoint_A_events_present",
        "all_endpoint_event_logs_have_mag_after_accel",
        "decomposition_identity_H_enclosed",
        "decomposition_identity_A_enclosed",
    ):
        if s.get(key) is not True:
            f.append(f"smoke lost {key}")
    if s.get("favorable_frontend_successor_selected") is not False:
        f.append("smoke selected a favorable front-end successor")
    if int(s.get("samples_executed", 0)) != 1 or int(s.get("endpoint_branches", 0)) < 1:
        f.append("typed execution smoke did not execute one complete sample")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "typed_execution_kernel_ready": d["typed_execution_kernel_ready"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
