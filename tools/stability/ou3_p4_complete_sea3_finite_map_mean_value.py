#!/usr/bin/env python3
"""Finite-map mean-value bridge for the canonical complete-SEA3 P4 theorem.

This module does not define a new Lyapunov metric and does not generate a source
word.  It supplies the theorem bridge that lets the retained outward AD event
machinery enclose the *finite physical error map* required by the paper.

Let F_W(0,zeta)=0 be the homogeneous physical true-minus-estimated complete-word
map on a star-shaped finite-error cell X.  If one SAME-HISTORY complete-SEA3
construction outward-encloses every ordinary/Clarke generalized Jacobian along
all radial segments t e, 0<=t<=1, by a matrix family J_W, then the generalized
mean-value theorem gives

    F_W(e,zeta) in J_W e,     e in X.

Thus the paper's finite quadratic endpoint condition is certified directly by
one full-matrix enclosure

    rho M_0 - J_W^T M_N J_W > 0,       0<rho<1,

where M_0 and M_N are the source-indexed endpoint information metrics from the
same SEA3 word.  This is a finite-state statement, not a differential-metric
replacement theorem.

Exactly the same construction applies to every literal prefix F_ell.  A finite
prefix gain kappa_V is certified by

    kappa_V M_0 - J_ell^T M_ell J_ell >= 0

for every 0<=ell<N_W, together with proof that each finite prefix remains in the
same certified chart/source domain.  The H18->A21 map may be rectangular; the
matrix test below therefore permits J to have a different output dimension.

The interval matrix supplied to these routines must come from one correlated
complete-SEA3 source cell.  Independently boxing P/H/R/K, tuner values, R_S,
event timing, or successive Jacobians is expressly not justified by this
bridge.  Every due S=0 event with the actual applied anisotropic SpectralMSE R_S
must already be inside the literal word whose Jacobian is enclosed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_word as DWRD

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_FINITE_MAP_MEAN_VALUE_BRIDGE_V1"


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    rows = len(A)
    cols = len(A[0]) if rows else 0
    if any(len(row) != cols for row in A):
        raise ValueError("ragged interval matrix")
    return rows, cols


def _scale(A: Sequence[Sequence[Interval]], c: float) -> IntervalMatrix:
    if not math.isfinite(float(c)):
        raise ValueError("finite scalar required")
    ci = Interval.outward_bounds(float(c), float(c))
    return [[ci * x for x in row] for row in A]


def finite_quadratic_margin(
    M0: Sequence[Sequence[Interval]],
    M1: Sequence[Sequence[Interval]],
    J: Sequence[Sequence[Interval]],
    factor: float,
) -> IntervalMatrix:
    """Return ``factor*M0 - J^T*M1*J`` with outward interval arithmetic.

    ``J`` may be rectangular, which is required for the H18->A21 hybrid.
    The caller is responsible for same-source provenance of M0/M1/J.
    """
    n0, n0b = _shape(M0)
    n1, n1b = _shape(M1)
    jr, jc = _shape(J)
    if n0 == 0 or n0 != n0b or n1 == 0 or n1 != n1b:
        raise ValueError("endpoint metrics must be nonempty square matrices")
    if (jr, jc) != (n1, n0):
        raise ValueError("Jacobian shape must map input metric dimension to output metric dimension")
    if not (math.isfinite(float(factor)) and float(factor) > 0.0):
        raise ValueError("positive finite factor required")
    pullback = matrix_mul(matrix_mul(matrix_transpose(J), M1), J)
    return matrix_symmetric_hull(matrix_sub(_scale(M0, factor), pullback))


def certify_strict_endpoint_contraction(
    M0: Sequence[Sequence[Interval]],
    MN: Sequence[Sequence[Interval]],
    JN: Sequence[Sequence[Interval]],
    rho: float,
) -> tuple[bool, list[float]]:
    """Strict full-matrix endpoint test for the finite physical map."""
    if not (0.0 < float(rho) < 1.0):
        raise ValueError("endpoint rho must lie strictly in (0,1)")
    margin = finite_quadratic_margin(M0, MN, JN, rho)
    ok, pivots = symmetric_positive_definite_ldlt(margin)
    return bool(ok), [float(x) for x in pivots]


def certify_prefix_gain(
    M0: Sequence[Sequence[Interval]],
    Mell: Sequence[Sequence[Interval]],
    Jell: Sequence[Sequence[Interval]],
    kappa_v: float,
) -> tuple[bool, list[float]]:
    """Full-matrix finite prefix-gain test for one literal prefix.

    A strict positive enclosure is used in computation.  It is stronger than
    the paper's non-strict prefix inequality and avoids an uncertified
    semidefinite numerical boundary.
    """
    if not (math.isfinite(float(kappa_v)) and float(kappa_v) >= 1.0):
        raise ValueError("prefix gain kappa_V must be finite and at least one")
    margin = finite_quadratic_margin(M0, Mell, Jell, kappa_v)
    ok, pivots = symmetric_positive_definite_ldlt(margin)
    return bool(ok), [float(x) for x in pivots]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    cayley = CAYLEY.build(path)
    events = EVENTS.build(path)
    word = DWRD.build(path)
    failures = (
        [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
        + [f"event AD: {x}" for x in EVENTS.validate(events)]
        + [f"word AD: {x}" for x in DWRD.validate(word)]
    )
    if failures:
        raise RuntimeError(f"finite-map mean-value prerequisites failed: {failures}")

    projection = bool(
        events["A21_bias_projection_generalized_Jacobian_available"]
        and events["A21_bias_projection_same_source_true_bias_required"]
        and word["A21_bias_projection_generalized_Jacobian_available"]
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "paper_finite_state_quadratic_theorem_retained": True,
        "finite_map_not_differential_metric_replacement": True,
        "generalized_mean_value_inclusion": "F(e,zeta)-F(0,zeta) in J_X(zeta)*e",
        "homogeneous_zero_error_invariance_required": True,
        "finite_error_cell_star_shaped_about_zero_required": True,
        "same_source_correlated_generalized_Jacobian_required": True,
        "same_source_endpoint_metrics_required": True,
        "same_source_prefix_metrics_required": True,
        "all_radial_segment_Jacobians_must_be_enclosed": True,
        "Clarke_generalized_Jacobian_handles_A21_projection": projection,
        "A21_projection_assumed_inactive": False,
        "rectangular_H18_to_A21_Jacobian_supported": True,
        "endpoint_matrix_test": "rho*M0-J_N^T*M_N*J_N > 0",
        "prefix_matrix_test": "kappa_V*M0-J_ell^T*M_ell*J_ell >= 0",
        "prefix_gain_required_for_every_literal_prefix": True,
        "prefix_domain_retention_required_for_every_literal_prefix": True,
        "actual_applied_RS_required_inside_same_word": True,
        "actual_RS_provenance_token": EVENTS.ACTUAL_RS_PROVENANCE,
        "independent_P_H_R_K_boxes_authorized": False,
        "independent_tuner_RS_schedule_authorized": False,
        "packetwise_remainder_sum_authorized": False,
        "finite_harmonic_or_replay_source_authorized": False,
        "source_uniform_complete_word_generalized_Jacobian_enclosed": False,
        "source_uniform_endpoint_finite_map_closed": False,
        "source_uniform_all_prefix_gains_closed": False,
        "source_uniform_all_prefix_domains_closed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "construct one source-correlated COMPLETE_SEA3_NORMAL_LIVE_WORD AD enclosure carrying state, P/H/R, "
            "committed tuner schedule, actual applied R_S, event timing and the A21 projection generalized Jacobian; "
            "use the resulting endpoint and every prefix Jacobian in the full-matrix tests above"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "paper_finite_state_quadratic_theorem_retained",
        "finite_map_not_differential_metric_replacement",
        "homogeneous_zero_error_invariance_required",
        "finite_error_cell_star_shaped_about_zero_required",
        "same_source_correlated_generalized_Jacobian_required",
        "same_source_endpoint_metrics_required", "same_source_prefix_metrics_required",
        "all_radial_segment_Jacobians_must_be_enclosed",
        "Clarke_generalized_Jacobian_handles_A21_projection",
        "rectangular_H18_to_A21_Jacobian_supported",
        "prefix_gain_required_for_every_literal_prefix",
        "prefix_domain_retention_required_for_every_literal_prefix",
        "actual_applied_RS_required_inside_same_word",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "A21_projection_assumed_inactive", "independent_P_H_R_K_boxes_authorized",
        "independent_tuner_RS_schedule_authorized", "packetwise_remainder_sum_authorized",
        "finite_harmonic_or_replay_source_authorized",
        "source_uniform_complete_word_generalized_Jacobian_enclosed",
        "source_uniform_endpoint_finite_map_closed", "source_uniform_all_prefix_gains_closed",
        "source_uniform_all_prefix_domains_closed", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("actual_RS_provenance_token") != EVENTS.ACTUAL_RS_PROVENANCE:
        f.append("actual R_S provenance changed")
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
        "qualification": d["qualification"],
        "finite_map_bridge": d["finite_map_not_differential_metric_replacement"],
        "endpoint_closed": d["source_uniform_endpoint_finite_map_closed"],
        "prefix_gains_closed": d["source_uniform_all_prefix_gains_closed"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
