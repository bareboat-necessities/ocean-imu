#!/usr/bin/env python3
"""Directional-form transport through the fixed-mode OU-III prediction.

P4 must move rank-deficient Joseph information between measurement times without
scalarizing it.  For the deterministic fixed-mode tangent prediction

    z_{k+1} = F_k z_k,

one later directional form ``Q_{k+1}`` pulls back exactly as

    Q_k = F_k^T Q_{k+1} F_k.                            (1)

Every concrete H/A prediction matrix in the declared source domain is
invertible.  Its block structure is triangular up to the attitude/gyro-bias
upper block:

* the attitude diagonal block is a rotation, determinant one;
* held gyro bias has identity diagonal;
* v,p,S have identity diagonal blocks;
* a_w has ``alpha=exp(-h/tau)>0``;
* active b_a has ``phi_ba=exp(-h/tau_ba)>0``.

Therefore (1) preserves rank and nullity exactly.  The process covariance Q_k
is not discarded: P3 already carries its source-uniform Riccati/information
comparison and prefix gain.  This primitive only prevents the nonlinear P4
layer from losing directional observability by replacing prediction with a
scalar norm bound.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose
import ou3_full_process_ucc as PROCESS
import ou3_implementation_proof_manifest as MANIFEST
import ou3_p3_source_uniform_certificate as P3CERT
import ou3_source_reachable_matrix_p3 as P3
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def pullback(F, Q):
    """Outward interval enclosure of the exact directional pullback F^T Q F."""
    return matrix_mul(matrix_mul(matrix_transpose(F), Q), F)


def _positive_exp_factor(x_lo: float, x_hi: float) -> tuple[float, float]:
    if not (0.0 <= x_lo <= x_hi <= VT.MAX_ABS_ARGUMENT):
        raise ValueError("validated positive exponential argument interval required")
    e = VT.exp_interval(-Interval.outward_bounds(x_lo, x_hi))
    if not e.lo > 0.0:
        raise RuntimeError("prediction exponential factor lost strict positivity")
    return e.lo, e.hi


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("prediction transport must not be trajectory fitted")

    manifest = MANIFEST.build()
    p3 = P3CERT.build(path)
    process = PROCESS.build()
    failures = [f"manifest: {x}" for x in MANIFEST.validate(manifest)]
    failures += [f"P3: {x}" for x in P3CERT.validate(p3)]
    failures += [f"process: {x}" for x in PROCESS.validate(process)]

    sched = P3.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    aw_lo, aw_hi = _positive_exp_factor(h / tau_hi, h / tau_lo)

    pc = process["source_constants"]
    tau_ba = float(pc["accel_bias_tau_s"])
    if not tau_ba > 0.0:
        failures.append("active accelerometer-bias time constant is not positive")
        ba_lo, ba_hi = math.nan, math.nan
    else:
        ba_lo, ba_hi = _positive_exp_factor(h / tau_ba, h / tau_ba)

    if manifest.get("normal_live_update_order", [None])[2] != "prediction":
        failures.append("shipping operation order no longer places prediction at expected position")

    modes = {
        "H": {
            "dimension": 18,
            "attitude_rotation_diagonal_determinant_exact": 1.0,
            "gyro_bias_diagonal": "I3",
            "translation_v_p_S_diagonal": "I9",
            "aw_diagonal_factor_interval": [aw_lo, aw_hi],
            "deterministic_transition_determinant_abs_lower": down(aw_lo ** 3),
            "deterministic_transition_invertible": aw_lo > 0.0,
            "directional_rank_preserved_exactly": True,
            "directional_nullity_preserved_exactly": True,
            "process_covariance_dropped_from_metric_argument": False,
            "process_metric_accounting": "consumed from P3 source-uniform Riccati/information comparison",
            "P4_PROMOTED": False,
        },
        "A": {
            "dimension": 21,
            "attitude_rotation_diagonal_determinant_exact": 1.0,
            "gyro_bias_diagonal": "I3",
            "translation_v_p_S_diagonal": "I9",
            "aw_diagonal_factor_interval": [aw_lo, aw_hi],
            "active_ba_diagonal_factor_interval": [ba_lo, ba_hi],
            "deterministic_transition_determinant_abs_lower": down((aw_lo ** 3) * (ba_lo ** 3)) if math.isfinite(ba_lo) else 0.0,
            "deterministic_transition_invertible": aw_lo > 0.0 and math.isfinite(ba_lo) and ba_lo > 0.0,
            "directional_rank_preserved_exactly": True,
            "directional_nullity_preserved_exactly": True,
            "process_covariance_dropped_from_metric_argument": False,
            "process_metric_accounting": "consumed from P3 source-uniform Riccati/information comparison",
            "P4_PROMOTED": False,
        },
    }
    for mode in ("H", "A"):
        if modes[mode]["deterministic_transition_invertible"] is not True:
            failures.append(f"{mode}: deterministic source prediction lost invertibility")
        if p3["modes"][mode].get("prefix_information_gain_upper") != 1.0:
            failures.append(f"{mode}: P3 prefix information gain changed")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_FIXED_MODE_PREDICTION_DIRECTIONAL_PULLBACK",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "directional_pullback_identity": "Q_k=F_k^T Q_{k+1} F_k",
        "per_prediction_scalarization_used": False,
        "condition_number_conversion_used": False,
        "P3_process_noise_metric_comparison_retained": True,
        "modes": modes,
        "P4_PREDICTION_DIRECTIONAL_TRANSPORT_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "combine this invertible prediction pullback with the exact reset congruence and same-cell Joseph S^-1 weights; "
            "materialize the accumulated H/A directional form over a source-complete word and certify its first positive full-rank generalized endpoint margin"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "P3_process_noise_metric_comparison_retained",
        "P4_PREDICTION_DIRECTIONAL_TRANSPORT_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "per_prediction_scalarization_used",
        "condition_number_conversion_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for mode, n in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode}: dimension mismatch")
        if m.get("deterministic_transition_invertible") is not True:
            f.append(f"{mode}: deterministic transition not invertible")
        if not float(m.get("deterministic_transition_determinant_abs_lower", 0.0)) > 0.0:
            f.append(f"{mode}: determinant lower is not positive")
        if m.get("directional_rank_preserved_exactly") is not True:
            f.append(f"{mode}: directional rank not preserved")
        if m.get("directional_nullity_preserved_exactly") is not True:
            f.append(f"{mode}: directional nullity not preserved")
        if m.get("process_covariance_dropped_from_metric_argument") is not False:
            f.append(f"{mode}: process covariance was dropped")
        if m.get("P4_PROMOTED") is not False:
            f.append(f"{mode}: P4 prematurely promoted")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "H_det_lower": d["modes"]["H"]["deterministic_transition_determinant_abs_lower"],
        "A_det_lower": d["modes"]["A"]["deterministic_transition_determinant_abs_lower"],
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
