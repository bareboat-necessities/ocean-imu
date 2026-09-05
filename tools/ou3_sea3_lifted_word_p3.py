#!/usr/bin/env python3
"""Canonical architecture for the OU-III SEA3 finite-word linear P3 proof.

This module intentionally contains no per-sample Riccati-lower propagation.
The previous one-step and commit-segment experiments failed for the same
structural reason: the integrated-OU [v,p,S,a_w] process covariance is almost
rank deficient over one 5 ms sample, so demanding a useful point lower after
small recursive steps throws away the finite-horizon controllability that P3 is
supposed to use.

The canonical proof is instead a lifted finite-dimensional Gaussian word.
Choose a finite set of independent standardized process-noise modes zeta over
one recurrent Normal-Live word.  Condition the initial state, every unselected
process input, and every nuisance noise source as known.  Extra conditioning can
only reduce conditional covariance.  For the selected modes,

    x_W = B_W zeta,
    y_W = A_W zeta + nu,          Cov(nu)=R_W,

so the exact selected-mode posterior is

    C_zeta = (I + A_W^T R_W^-1 A_W)^-1

and its endpoint covariance contribution is

    Omega_sel = B_W C_zeta B_W^T.

Therefore the full shipping word injection satisfies

    Omega_W >= Omega_sel.

The moving shipping covariance metric still supplies samplewise non-expansion:
for e+ = C_k e and P+ = C_k P C_k^T + Omega_k with Omega_k >= 0,

    e+^T P+^-1 e+ <= e^T P^-1 e.

Strictness is required only over the recurrent word.  With a validated endpoint
covariance upper Pbar, the quantitative gate is the full-matrix comparison

    Omega_sel >= delta Pbar,      delta >= 1e-18.

No determinant/trace conversion is canonical.  No scalar beta attenuation is
canonical.  No source-history graph, predecessor enumeration, endpoint tau
partition, per-sample SPD lower, or commit-aligned source word is canonical.
The only remaining numerical work is to enclose the two small lifted matrices
B_W and J_W=A_W^T R^-1 A_W directly over the SEA3 dynamic invariant, then use
validated LDLT/generalized-eigenvalue tests on the resulting full matrices.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import ou3_full_process_ucc as PROCESS
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_LIFTED_FINITE_WORD_P3_ARCHITECTURE"
USEFUL_GATE = 1.0e-18


def _positive(x: Any, label: str) -> float:
    v = float(x)
    if not (math.isfinite(v) and v > 0.0):
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return v


def _load_tube(path: Path, tube_path: Path | None) -> dict:
    tube = TUBE.build(path) if tube_path is None else json.loads(
        Path(tube_path).read_text(encoding="utf-8")
    )
    failures = TUBE.validate(tube)
    if failures:
        raise RuntimeError(f"SEA3 covariance-upper prerequisite failed: {failures}")
    return tube


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("lifted P3 architecture may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    vector = VECTOR.build()
    process = PROCESS.build()
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    prereq_failures = (
        [f"dynamic: {x}" for x in DYNAMIC.validate(dynamic)]
        + [f"vector: {x}" for x in VECTOR.validate(vector)]
        + [f"process: {x}" for x in PROCESS.validate(process)]
        + [f"translation: {x}" for x in TRANS.validate(trans)]
    )
    if prereq_failures:
        raise RuntimeError(f"lifted P3 prerequisites failed: {prereq_failures}")

    tube = _load_tube(path, tube_path)
    live = domain["normal_live"]
    pe_window = _positive(live["vector_pe_recurrence_window_s"], "PE window")
    timing = tube["covariance_memory"]
    covariance_memory = _positive(
        timing["covariance_memory_window_s_upper"], "covariance memory"
    )
    word_horizon = max(pe_window, covariance_memory)

    alpha6 = _positive(
        vector["gyro_bias_two_packet"]["alpha_6_information_lower"],
        "vector UCO alpha6",
    )
    q_trans = _positive(
        trans["process_ucc"]["Q_axis_lambda_min_lower"],
        "translation UCC floor",
    )
    q_att = _positive(
        process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"],
        "attitude/gyro-bias process floor",
    )
    q_ba = _positive(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "accelerometer-bias process floor",
    )

    modes: dict[str, dict[str, Any]] = {}
    for mode, dim, selected_dim in (("H", 18, 18), ("A", 21, 21)):
        pdiag = [float(x) for x in tube["modes"][mode]["Pbar_diagonal_variance_upper"]]
        if len(pdiag) != dim or any(not (math.isfinite(x) and x > 0.0) for x in pdiag):
            raise RuntimeError(f"{mode} endpoint covariance upper is invalid")
        modes[mode] = {
            "dimension": dim,
            "selected_process_mode_dimension": selected_dim,
            "selected_process_mode_layout": {
                "attitude_gyro_bias": {
                    "modes": 6,
                    "construction": "two independent temporal process modes per physical axis spanning theta and b_g",
                },
                "translation": {
                    "modes": 12,
                    "construction": "four independent source-independent temporal process modes per axis spanning v,p,S,a_w over the full recurrent word",
                },
                "accelerometer_bias": {
                    "modes": 0 if mode == "H" else 3,
                    "construction": "one independent active-bias process mode per axis" if mode == "A" else "inactive in H mode",
                },
            },
            "Pbar_diagonal_variance_upper": pdiag,
            "Pbar_lambda_max_trace_upper": float(tube["modes"][mode]["Pbar_lambda_max_trace_upper"]),
            "lifted_endpoint_map_B_closed": False,
            "lifted_measurement_information_J_upper_closed": False,
            "selected_mode_posterior_lower_closed": False,
            "Omega_selected_full_matrix_lower_closed": False,
            "generalized_eigenvalue_comparison_closed": False,
            "relative_Riccati_injection_margin_lower": None,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "LIFTED_FINITE_WORD_SELECTED_PROCESS_MODES",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_history_frontier_consumed": False,
        "endpoint_tau_partition_consumed": False,
        "commit_aligned_source_word_consumed": False,
        "per_sample_Riccati_lower_propagation_used": False,
        "per_sample_SPD_lower_required": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "gain_history_enumeration_used": False,
        "word_horizon_s": word_horizon,
        "samplewise_nonexpansion": {
            "closed": True,
            "identity": "P_plus = C_k P C_k^T + Omega_k, Omega_k >= 0",
            "consequence": "e_plus^T P_plus^-1 e_plus <= e^T P^-1 e",
            "strict_per_sample_margin_required": False,
        },
        "lifted_word_identity": {
            "selected_modes": "zeta ~ N(0,I)",
            "endpoint": "x_W = B_W zeta",
            "measurements": "y_W = A_W zeta + nu, Cov(nu)=R_W",
            "posterior_modes": "C_zeta = (I + A_W^T R_W^-1 A_W)^-1",
            "selected_endpoint_injection": "Omega_sel = B_W C_zeta B_W^T",
            "conditioning_monotonicity": "Omega_W >= Omega_sel after conditioning initial state, unselected process inputs, and nuisance sources as known",
            "canonical_gate": "validated LDLT of Omega_sel - delta*Pbar with delta >= 1e-18",
        },
        "source_uniform_matrix_enclosure_contract": {
            "B_W": "direct full-word enclosure over the compact SEA3 dynamic invariant; source-independent temporal modes; no path enumeration",
            "J_W": "direct full lifted information upper A_W^T R^-1 A_W including every admissible accelerometer, S=0 and magnetic packet",
            "posterior": "validated matrix inverse/factorization of I+J_upper; no scalar trace attenuation",
            "endpoint_comparison": "full H18/A21 matrix generalized-eigenvalue/LDLT comparison against endpoint Pbar",
        },
        "retained_strict_subcertificates": {
            "vector_UCO_alpha6_lower": alpha6,
            "translation_process_UCC_one_axis_strict_floor": q_trans,
            "attitude_gyro_bias_process_UCC_strict_floor": q_att,
            "active_accelerometer_bias_process_strict_floor": q_ba,
        },
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "P3_ARCHITECTURE_READY": True,
        "P3_QUANTITATIVE_WORD_MATRIX_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "build exactly two validated full-word matrix enclosures, B_W and J_W, then form Omega_sel and run the full-matrix H18/A21 generalized-eigenvalue gate; do not introduce another recursive covariance or source-word proof route"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "LIFTED_FINITE_WORD_SELECTED_PROCESS_MODES":
        f.append("wrong canonical P3 architecture")
    for key in (
        "source_generated_not_trajectory_fit",
        "SEA3_dynamic_source_consumed",
        "P3_ARCHITECTURE_READY",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "old_P2_history_frontier_consumed",
        "endpoint_tau_partition_consumed", "commit_aligned_source_word_consumed",
        "per_sample_Riccati_lower_propagation_used", "per_sample_SPD_lower_required",
        "determinant_trace_scalarization_used", "scalar_information_beta_used",
        "gain_history_enumeration_used", "P3_QUANTITATIVE_WORD_MATRIX_CLOSED",
        "P3_CANONICAL_PASS", "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("canonical useful gate changed")
    ne = d.get("samplewise_nonexpansion", {})
    if ne.get("closed") is not True or ne.get("strict_per_sample_margin_required") is not False:
        f.append("samplewise nonexpansion contract changed")
    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if row.get("dimension") != dim or row.get("selected_process_mode_dimension") != dim:
            f.append(f"{mode} selected lifted basis is not full-dimensional")
        if row.get("pass") is not False:
            f.append(f"{mode} promoted before quantitative matrix closure")
        if row.get("relative_Riccati_injection_margin_lower") is not None:
            f.append(f"{mode} emitted a margin before quantitative matrix closure")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P3_architecture"],
        "word_horizon_s": d["word_horizon_s"],
        "H_selected_modes": d["modes"]["H"]["selected_process_mode_dimension"],
        "A_selected_modes": d["modes"]["A"]["selected_process_mode_dimension"],
        "P3_ARCHITECTURE_READY": d["P3_ARCHITECTURE_READY"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
