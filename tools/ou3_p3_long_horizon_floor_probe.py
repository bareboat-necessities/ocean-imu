#!/usr/bin/env python3
"""Falsifiable P3 probe for the zero-start segment-floor hypothesis.

PR #473 established that canonical P3 is about nine orders below the frozen
1e-18 useful gate even after replacing the isotropic identity floor by the
available anisotropic endpoint floor.  Enclosure refinements are therefore not
large enough to close the observed deficit.  The remaining concrete hypothesis
is structural: ``common_boundary_floor`` throws away all covariance older than
one 13..26-sample source segment by restarting from P=0.

This producer tests that hypothesis without promoting P3 and without tuning any
proof constant.  It keeps one frozen P2 source cell, starts from zero once, and
propagates the already validated full 4x4 Loewner lower all the way to the exact
source-uniform covariance-word sample target used by the P2-V1 history quotient.
At selected horizons it recomputes the same relative injection margin against
``HIST.endpoint_phase_upper`` used by canonical P3.

The fixed source/x interval does not change inside this experiment.  Rebuilding
its transition matrix, process covariance and physical measurement variances at
every one of hundreds of samples is therefore pure repeated arithmetic, not a
proof obligation.  The probe precomputes those immutable objects once per x
subcell, checks the first prepared step bit-for-bit against ``SEG.one_step``,
and then iterates the identical interval/Loewner algebra.  No enclosure is
cached across changing covariance arguments.

The fixed-source construction is deliberately only a diagnostic.  A favorable
result means that a source-varying full-word lower should be built over the
exact P2-V1 path language.  An unfavorable result rejects the zero-reset
explanation even before paying for that substantially larger construction.
It never sets P3/P4/P5 promotion flags.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_segment as SEG
import ou3_p3_frozen_full_matrix_translation as FROZEN
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_p3_p2_v1_stage_phase_translation as STAGE
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DEFAULT_SOURCE_NODE = 137
REFERENCE_SEGMENT_SAMPLES = 13
_COORDINATES = ("v", "p", "S", "a_w")


def _copy_matrix(P):
    return [[P[i][j] for j in range(4)] for i in range(4)]


def _same_interval_matrix(A, B) -> bool:
    if len(A) != len(B):
        return False
    return all(
        len(A[i]) == len(B[i])
        and all(
            float(A[i][j].lo) == float(B[i][j].lo)
            and float(A[i][j].hi) == float(B[i][j].hi)
            for j in range(len(A[i]))
        )
        for i in range(len(A))
    )


def _prepare_step_kernel(node: dict, x: Interval, h: float) -> dict:
    """Precompute only source/x objects that SEG.one_step recomputes verbatim."""
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    sigma2_lower = BASE.down(sigma.lo * sigma.lo)
    if not sigma2_lower > 0.0:
        raise RuntimeError("source-cell sigma lower lost positivity")
    Fm = FROZEN._transition(x)
    Ft = FROZEN._transpose(Fm)
    Qz = SEG._scale(FROZEN._scaled_Q(x), SEG._point(sigma2_lower))
    R_aw, R_S_z = SEG._physical_measurement_variances(node, h)
    return {
        "F": Fm,
        "Ft": Ft,
        "Qz": Qz,
        "R_aw": float(R_aw),
        "R_S_z": float(R_S_z),
    }


def _prepared_one_step(P, kernel: dict):
    """Exact SEG.one_step algebra with immutable source/x terms precomputed."""
    L0, _ = SEG._common_point_lower(P)
    prediction_interval = FROZEN._sym(
        FROZEN._add(
            FROZEN._mul(FROZEN._mul(kernel["F"], L0), kernel["Ft"]),
            kernel["Qz"],
        )
    )
    pred, _ = SEG._common_point_lower(prediction_interval)
    post_aw, _ = SEG._point_measurement_lower(pred, 3, kernel["R_aw"])
    post_s, _ = SEG._point_measurement_lower(post_aw, 2, kernel["R_S_z"])
    return post_s


def _physical_diagonal_ratios(Pz, h: float, upper: list[float]) -> list[float]:
    if len(upper) != 4:
        raise ValueError("four translation covariance upper entries required")
    d = (h, h * h, h * h * h, 1.0)
    out = []
    for i in range(4):
        u = float(upper[i])
        if not (math.isfinite(u) and u > 0.0):
            raise ValueError("positive finite translation covariance upper required")
        physical_lo = BASE.down(BASE.down(d[i] * d[i]) * float(Pz[i][i].lo))
        out.append(BASE.down(max(0.0, physical_lo) / u))
    return out


def _horizons(word_samples: int) -> tuple[int, ...]:
    if word_samples < REFERENCE_SEGMENT_SAMPLES:
        raise ValueError("full-word target shorter than reference segment")
    candidates = (
        REFERENCE_SEGMENT_SAMPLES,
        26,
        52,
        104,
        208,
        416,
        int(word_samples),
    )
    return tuple(sorted({k for k in candidates if 0 < k <= int(word_samples)}))


def probe_source(source_node: int = DEFAULT_SOURCE_NODE,
                 domain_path: Path = DEFAULT_DOMAIN,
                 x_subcells: int = SEG.X_SUBCELLS) -> dict:
    path = Path(domain_path).resolve()
    fr = HIST.frontier_runtime(path)
    rt = fr["rt"]
    s = int(source_node)
    if not 0 <= s < len(rt["nodes"]):
        raise IndexError("source node outside frozen P2 correlation partition")
    if s not in fr["endpoint_frontiers"]:
        raise RuntimeError("probe source is not a full-word P2 endpoint")

    target_samples = int(fr["target"]["target_samples"])
    horizons = _horizons(target_samples)
    upper_row = HIST.endpoint_phase_upper(s, 0, fr)
    upper = list(map(float, upper_row["Sigma_translation_diagonal_upper_envelope"]))
    h = float(rt["clock"]["dt_binary32_s"])
    node = rt["nodes"][s]
    xcells = SEG._node_x_subcells(node, h, int(x_subcells))
    if not xcells:
        raise RuntimeError("probe source emitted no x subcells")

    zero = FROZEN._mat_zero()
    kernels = [_prepare_step_kernel(node, x, h) for x in xcells]
    prepared_equivalence_checked = 0
    states = []
    # This exact endpoint comparison makes the cache a performance transform,
    # not a second numerical route.  The expensive source/x terms are then
    # reused; every covariance-dependent collapse/update is still recomputed.
    for x, kernel in zip(xcells, kernels):
        direct = SEG.one_step(zero, node, x, h)
        prepared = _prepared_one_step(zero, kernel)
        if not _same_interval_matrix(direct, prepared):
            raise RuntimeError("prepared long-horizon kernel changed SEG.one_step enclosure")
        prepared_equivalence_checked += 1
        states.append(prepared)

    rows = []
    horizon_set = set(horizons)
    if 1 in horizon_set:
        raise RuntimeError("probe horizons unexpectedly include the precomputed first sample")

    for sample in range(2, target_samples + 1):
        for j, kernel in enumerate(kernels):
            states[j] = _prepared_one_step(states[j], kernel)
        if sample not in horizon_set:
            continue

        deltas = [STAGE._certified_delta(P, h, upper) for P in states]
        diag_ratios = [_physical_diagonal_ratios(P, h, upper) for P in states]
        delta = min(deltas)
        worst_flat = min(
            (ratio, j, i)
            for j, row in enumerate(diag_ratios)
            for i, ratio in enumerate(row)
        )
        _, worst_x, worst_coord = worst_flat
        rows.append({
            "samples": int(sample),
            "time_s_lower": BASE.down(sample * math.nextafter(h, -math.inf)),
            "time_s_upper": BASE.up(sample * math.nextafter(h, math.inf)),
            "translation_relative_Riccati_injection_margin_lower": float(delta),
            "minimum_physical_diagonal_ratio_lower": float(worst_flat[0]),
            "minimum_physical_diagonal_ratio_coordinate": _COORDINATES[worst_coord],
            "minimum_physical_diagonal_ratio_x_subcell": int(worst_x),
            "all_x_subcells_checked": len(states),
        })

    if not rows or rows[0]["samples"] != REFERENCE_SEGMENT_SAMPLES:
        raise RuntimeError("probe did not emit the reference segment horizon")
    if rows[-1]["samples"] != target_samples:
        raise RuntimeError("probe did not reach the exact full-word target")

    baseline = float(rows[0]["translation_relative_Riccati_injection_margin_lower"])
    terminal = float(rows[-1]["translation_relative_Riccati_injection_margin_lower"])
    growth = math.inf if baseline == 0.0 and terminal > 0.0 else (
        terminal / baseline if baseline > 0.0 else 1.0
    )
    orders = math.inf if math.isinf(growth) else (
        math.log10(growth) if growth > 0.0 else -math.inf
    )

    if terminal >= BASE.MIN_USEFUL_DELTA:
        classification = "FIXED_SOURCE_FULL_WORD_REACHES_CANONICAL_GATE"
        next_obligation = (
            "replace the one-segment zero reset by a source-varying full-word Loewner lower "
            "over the exact P2-V1 path language, then rerun canonical P3 without changing its gate"
        )
    elif growth >= 1.0e6:
        classification = "FIXED_SOURCE_FULL_WORD_MOVES_BY_MANY_ORDERS"
        next_obligation = (
            "construct the source-varying full-word lower over exact P2-V1 paths and determine "
            "whether the remaining gap survives source changes; do not return to enclosure tuning"
        )
    else:
        classification = "ZERO_START_HORIZON_HYPOTHESIS_NOT_SUPPORTED_BY_FIXED_SOURCE_PROBE"
        next_obligation = (
            "reconsider the P3 theorem representation or relative covariance comparison; "
            "do not spend more work tightening the one-segment enclosure"
        )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_ZERO_START_LONG_HORIZON_STRUCTURAL_PROBE",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "canonical_gate_changed": False,
        "canonical_useful_gate": BASE.MIN_USEFUL_DELTA,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "fixed_source_diagnostic_only": True,
        "source_complete_full_word_lower_established": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "source_node": s,
        "tau_s": node["tau_s"],
        "sigma_filter_committed_mps2": node["sigma_filter_committed_mps2"],
        "R_S_filter_std": node["R_S_filter_std"],
        "x_subcells": len(xcells),
        "zero_start_used_once_at_word_start": True,
        "zero_reset_at_each_13_26_sample_segment_used": False,
        "prepared_source_x_kernel_used": True,
        "prepared_kernel_changes_covariance_enclosure": False,
        "prepared_kernel_exact_SEG_one_step_equivalence_subcells": prepared_equivalence_checked,
        "same_validated_SEG_one_step_algebra_used": True,
        "strongest_translation_measurements_every_sample_retained": True,
        "same_history_Sigma_upper_used": True,
        "Sigma_translation_diagonal_upper_envelope": upper,
        "full_word_target": dict(fr["target"]),
        "horizons": rows,
        "reference_segment_margin_lower": baseline,
        "full_word_margin_lower": terminal,
        "margin_growth_factor_lower": growth,
        "margin_growth_orders_log10": orders,
        "classification": classification,
        "next_obligation": next_obligation,
        "failures": [],
    }


def build(domain_path: Path = DEFAULT_DOMAIN,
          source_node: int = DEFAULT_SOURCE_NODE,
          x_subcells: int = SEG.X_SUBCELLS) -> dict:
    return probe_source(source_node, Path(domain_path).resolve(), int(x_subcells))


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_ZERO_START_LONG_HORIZON_STRUCTURAL_PROBE":
        f.append("wrong qualification")
    for key in (
        "source_only",
        "P2_correlation_interface_consumed",
        "fixed_source_diagnostic_only",
        "zero_start_used_once_at_word_start",
        "prepared_source_x_kernel_used",
        "same_validated_SEG_one_step_algebra_used",
        "strongest_translation_measurements_every_sample_retained",
        "same_history_Sigma_upper_used",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "canonical_gate_changed",
        "source_complete_full_word_lower_established",
        "P3_PROMOTED",
        "P4_PROMOTED",
        "P5_PROMOTED",
        "zero_reset_at_each_13_26_sample_segment_used",
        "prepared_kernel_changes_covariance_enclosure",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("prepared_kernel_exact_SEG_one_step_equivalence_subcells", 0)) != int(d.get("x_subcells", -1)):
        f.append("prepared kernel was not checked against SEG.one_step on every x subcell")
    if float(d.get("canonical_useful_gate", math.nan)) != 1.0e-18:
        f.append("canonical useful gate changed")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("probe lost frozen P2 V1 binding")
    target = d.get("full_word_target", {})
    target_samples = int(target.get("target_samples", 0) or 0)
    rows = list(d.get("horizons", []))
    if target_samples <= 0 or not rows:
        f.append("missing full-word target/horizon rows")
    elif int(rows[-1].get("samples", 0) or 0) != target_samples:
        f.append("long-horizon probe did not reach exact full-word target")
    upper = list(d.get("Sigma_translation_diagonal_upper_envelope", []))
    if len(upper) != 4 or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in upper):
        f.append("invalid same-history translation covariance upper")
    for row in rows:
        delta = float(row.get("translation_relative_Riccati_injection_margin_lower", math.nan))
        ratio = float(row.get("minimum_physical_diagonal_ratio_lower", math.nan))
        if not (math.isfinite(delta) and delta >= 0.0):
            f.append("non-finite horizon delta")
        if not (math.isfinite(ratio) and ratio >= 0.0):
            f.append("non-finite horizon diagonal ratio")
        if row.get("minimum_physical_diagonal_ratio_coordinate") not in _COORDINATES:
            f.append("invalid limiting coordinate")
    if not isinstance(d.get("classification"), str) or not d["classification"]:
        f.append("missing structural classification")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node", type=int, default=DEFAULT_SOURCE_NODE)
    ap.add_argument("--x-subcells", type=int, default=SEG.X_SUBCELLS)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.source_node, args.x_subcells)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "source_node": d["source_node"],
        "full_word_target_samples": d["full_word_target"]["target_samples"],
        "horizons": d["horizons"],
        "reference_segment_margin_lower": d["reference_segment_margin_lower"],
        "full_word_margin_lower": d["full_word_margin_lower"],
        "margin_growth_factor_lower": d["margin_growth_factor_lower"],
        "margin_growth_orders_log10": d["margin_growth_orders_log10"],
        "classification": d["classification"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
