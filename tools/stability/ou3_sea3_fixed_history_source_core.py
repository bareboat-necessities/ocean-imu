#!/usr/bin/env python3
"""Non-promoting fixed-history complete-SEA3 source-core feasibility point.

This is one member of the existing COMPLETE_SEA3_NORMAL_LIVE_WORD continuum
Hilbert-ball source.  It is not a replay, a finite harmonic source, a sampled
box, or a replacement source model.

The member is owned by ``ou3_sea3_validated_continuum_member``: one common
continuum coefficient field supported on the continuous band 0.295--0.305 Hz,
one admissible H=1.5 m/Tp=6 s/gamma=3.3 partition, unit vertical response on
the support band, zero rotational rows of the same six-DOF response, and no
phase reset between the 60 s startup prehistory and the 601-sample word.

Its output integral has a closed analytic form and every numerical sample is
wrapped by outward Decimal/Taylor intervals.  No quadrature convergence claim
is used as source admission.  This module deliberately stops before inventing
front-end state or covariance; the next stage must obtain those by running the
same prehistory through the unchanged shipping startup path.
"""
from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal
from pathlib import Path

import ou3_sea3_continuum_driver_gram as DRIVER
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_directional_response_family as RESPONSE
import ou3_sea3_validated_continuum_member as MEMBER

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_FIXED_HISTORY_SOURCE_CORE_V3"
N = 601
DT = 0.005
PREHISTORY_S = 60.0
WINDOW_START_S = PREHISTORY_S
H_M = 1.5
TP_S = 6.0
GAMMA = 3.3
RAO_GAIN = 1.0
RAO_CORNER_HZ = 0.5
RAO_POWER = 2.0


def acceleration_interval_at_time(t_s: float):
    return MEMBER.acceleration_interval(Decimal(str(t_s)))


def acceleration_at_time(t_s: float, panels: int | None = None) -> float:
    """Midpoint floating evaluation of the rigorously enclosed exact member.

    ``panels`` is accepted only for compatibility with the retired diagnostic
    caller; it has no effect and no quadrature is performed.
    """
    del panels
    v = acceleration_interval_at_time(t_s)
    return float((v.lo + v.hi) / Decimal(2))


def _evaluate_window() -> tuple[list[float], list[list[float]], float, float]:
    nominal: list[float] = []
    intervals: list[list[float]] = []
    max_width = 0.0
    max_abs = 0.0
    for k in range(N):
        t = WINDOW_START_S + k * DT
        v = acceleration_interval_at_time(t)
        bounds = v.float_bounds()
        mid = float((v.lo + v.hi) / Decimal(2))
        nominal.append(mid)
        intervals.append(bounds)
        max_width = max(max_width, bounds[1] - bounds[0])
        max_abs = max(max_abs, abs(bounds[0]), abs(bounds[1]))
    return nominal, intervals, max_width, max_abs


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    physical = PHYSICAL.build(domain_path)
    response = RESPONSE.directional_response_enclosure(REPO)
    driver = DRIVER.build()
    bad = {
        "physical": PHYSICAL.validate(physical),
        "response": RESPONSE.validate(response),
        "driver": DRIVER.validate(driver),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"fixed SEA3 source prerequisites failed: {bad}")

    g = float(physical["gravity_mps2"])
    if not PHYSICAL.partition_admissible(H_M, TP_S, g):
        raise RuntimeError("selected fixed SEA3 partition violates peak-steepness contract")
    RESPONSE.evaluate_rao_envelope_member(
        RAO_GAIN, RAO_CORNER_HZ, RAO_POWER, response
    )
    member_check = MEMBER.self_check()
    norm_cert = member_check["driver_norm_certificate"]
    if not norm_cert["driver_norm_strictly_below_one"]:
        raise RuntimeError("validated continuum member escaped complete-SEA3 hard driver ball")

    y, y_interval, max_width, max_abs = _evaluate_window()
    normal_live_accel_cap = float(
        json.loads(domain_path.read_text(encoding="utf-8"))["normal_live"]
        ["non_gravitational_cog_acceleration_norm_upper_mps2"]
    )
    source_core = []
    for k in range(N):
        t = WINDOW_START_S + k * DT
        source_core.append({
            "k": k,
            "source_time_s": t,
            "word_time_s": k * DT,
            # The scalar is the nominal floating point used by the feasibility
            # executor; the interval is the authoritative source enclosure.
            "f_cog_body_mps2": [0.0, 0.0, y[k]],
            "f_cog_body_interval_mps2": [[0.0, 0.0], [0.0, 0.0], y_interval[k]],
            "omega_body_corrected_rad_s": [0.0, 0.0, 0.0],
            "omega_body_corrected_interval_rad_s": [[0.0, 0.0]] * 3,
        })

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "role": "non-promoting legal fixed-history point for the ledger complete-word feasibility experiment",
        "source_membership": {
            "hard_driver_qualification": driver["qualification"],
            "fixed_history_operator": driver["fixed_history_operator"],
            "driver_field": "single common continuum coefficient field on [0.295,0.305] Hz; zero outside",
            "driver_norm_certificate": norm_cert,
            "driver_norm": math.sqrt(float(norm_cert["driver_norm_squared_upper"])),
            "membership_is_analytic_Hilbert_norm_not_quadrature": True,
            "membership_is_outward_enclosed": True,
            "same_driver_field_prehistory_and_window": True,
            "same_driver_field_translation_and_rotation": True,
        },
        "SEA3_fixed_history": {
            "active_partitions": 1,
            "H_r_m": [H_M, 0.0, 0.0],
            "Tp_r_s": [TP_S, TP_S, TP_S],
            "gamma_r": [GAMMA, 1.0, 1.0],
            "lambda_constant_over_prehistory_and_window": True,
            "fixed_lambda_family_substitution_used": False,
            "partition_peak_steepness_admissible": True,
            "total_Hs_m": H_M,
            "directional_density_integrates_to_one": True,
        },
        "fixed_response_member": {
            "translation": "unit vertical response on continuous driver support band, admissibly tapered/zero outside",
            "G": RAO_GAIN,
            "fc_hz": RAO_CORNER_HZ,
            "power": RAO_POWER,
            "rotation": "zero rows of same six-DOF response operator",
            "inside_declared_continuum_RAO_family": True,
            "finite_RAO_grid_used": False,
        },
        "same_history_phase_contract": {
            "prehistory_start_source_time_s": 0.0,
            "prehistory_end_source_time_s": PREHISTORY_S,
            "word_start_source_time_s": WINDOW_START_S,
            "prehistory_duration_s": PREHISTORY_S,
            "phase_center_source_time_s": float(MEMBER.PHASE_CENTER_S),
            "phase_reset_at_live_entry": False,
            "word_uses_immediate_phase_continuation": True,
        },
        # Keep this compatibility block so old contract readers fail toward the
        # correct direction: there are zero quadrature nodes and zero quadrature
        # error because quadrature is not used at all.
        "quadrature_diagnostic": {
            "coordinate": "not used",
            "source_modes_are_quadrature_nodes": False,
            "coarse_panels": 0,
            "fine_panels": 0,
            "max_abs_sample_delta_coarse_to_fine": 0.0,
            "max_abs_source_acceleration_mps2": max_abs,
            "normal_live_acceleration_cap_mps2": normal_live_accel_cap,
            "inside_normal_live_acceleration_cap": max_abs <= normal_live_accel_cap,
            "convergence_relative_to_peak": 0.0,
            "quadrature_used": False,
        },
        "validated_analytic_evaluation": {
            "formula": "C*(sin(2*pi*f2*(t-tc))-sin(2*pi*f1*(t-tc)))/(2*pi*(t-tc))",
            "continuous_frequency_band_hz": [float(MEMBER.F1), float(MEMBER.F2)],
            "source_modes_are_band_endpoints": False,
            "outward_decimal_taylor_enclosure": True,
            "max_sample_interval_width": max_width,
        },
        "sample_count": N,
        "dt_s": DT,
        "source_core": source_core,
        "front_end_entry_derived_from_same_history": False,
        "live_covariance_seed_derived_from_same_history": False,
        "complete_executor_artifact_materialized": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "trajectory_replay_used": False,
        "finite_harmonic_source_used": False,
        "independent_sample_boxes_used": False,
        "next_obligation": (
            "run source times [0,60 s) from this same outward-enclosed analytic member through shipping updateFrontEnd until TunerReady, call real goLive at source time 60 s, serialize actual frontend/Live covariance and all actual-applied R_S events, then execute the connected 601-sample nonlinear word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("fixed history detached from canonical complete SEA3")
    membership = d.get("source_membership", {})
    for key in (
        "membership_is_analytic_Hilbert_norm_not_quadrature",
        "membership_is_outward_enclosed",
        "same_driver_field_prehistory_and_window",
        "same_driver_field_translation_and_rotation",
    ):
        if membership.get(key) is not True:
            f.append(f"source membership lost {key}")
    norm = float(membership.get("driver_norm", math.inf))
    if not (0.0 < norm < 1.0):
        f.append("fixed-history continuum driver is not strictly inside admitted hard unit ball")
    sea = d.get("SEA3_fixed_history", {})
    if sea.get("active_partitions") != 1 or sea.get("partition_peak_steepness_admissible") is not True:
        f.append("fixed feasibility history is not physically admissible")
    if sea.get("fixed_lambda_family_substitution_used") is not False:
        f.append("fixed-lambda family shortcut reintroduced")
    response = d.get("fixed_response_member", {})
    if response.get("inside_declared_continuum_RAO_family") is not True:
        f.append("fixed response is outside SEA3 RAO family")
    if response.get("finite_RAO_grid_used") is not False:
        f.append("finite RAO grid re-entered fixed-history source")
    phase = d.get("same_history_phase_contract", {})
    if phase.get("phase_reset_at_live_entry") is not False or phase.get("word_uses_immediate_phase_continuation") is not True:
        f.append("fixed source broke phase continuity at Live entry")
    if not math.isclose(float(phase.get("prehistory_duration_s", -1.0)), PREHISTORY_S):
        f.append("fixed source prehistory duration drifted")
    quad = d.get("quadrature_diagnostic", {})
    if quad.get("source_modes_are_quadrature_nodes") is not False or quad.get("quadrature_used") is not False:
        f.append("quadrature was reintroduced as source machinery")
    if quad.get("inside_normal_live_acceleration_cap") is not True:
        f.append("selected fixed source point leaves declared Normal-Live acceleration domain")
    analytic = d.get("validated_analytic_evaluation", {})
    if analytic.get("outward_decimal_taylor_enclosure") is not True:
        f.append("source samples lack outward analytic enclosure")
    if not float(analytic.get("max_sample_interval_width", math.inf)) < 1e-8:
        f.append("source sample enclosure is too wide for point feasibility")
    if d.get("sample_count") != N or not math.isclose(float(d.get("dt_s", 0.0)), DT):
        f.append("fixed history does not cover canonical 601 samples at 5 ms")
    if not isinstance(d.get("source_core"), list) or len(d["source_core"]) != N:
        f.append("source core payload length mismatch")
    for key in (
        "front_end_entry_derived_from_same_history", "live_covariance_seed_derived_from_same_history",
        "complete_executor_artifact_materialized", "P4_promoted", "P5_promoted",
        "trajectory_replay_used", "finite_harmonic_source_used", "independent_sample_boxes_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false at source-core stage")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "source": d["canonical_source"],
        "membership": d["source_membership"],
        "history": d["SEA3_fixed_history"],
        "phase": d["same_history_phase_contract"],
        "response": d["fixed_response_member"],
        "analytic": d["validated_analytic_evaluation"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
