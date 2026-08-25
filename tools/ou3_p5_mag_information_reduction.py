#!/usr/bin/env python3
"""Exact magnetometer Joseph-information reduction for OU-III P5.

The effective-vector-input lemma already proves that the radial component of the
configured magnetometer residual is killed by the Kalman gain.  This producer
closes the corresponding *information* identity as well.

For the shipping magnetometer

    H = -[v]_x,       R = r I,
    S = H P H^T + r I,

let ``Pi = I-vhat vhat^T`` and split the exact finite-angle residual
``y=y_T+y_R`` into the tangent and radial subspaces of ``v``.  Since

    H^T v = 0,
    H P H^T v = 0,
    S v = r v,

both ``S`` and ``R`` leave ``span(v)`` and its orthogonal complement invariant.
The tangent linear residual ``h=H c`` has zero radial component, so the radial
component of ``eta=y-h`` is exactly ``y_R``.  Therefore the radial terms cancel
*exactly* in the Joseph information identity:

    y^T S^-1 y - eta^T R^-1 eta
      = y_T^T S^-1 y_T - eta_T^T R^-1 eta_T.             (1)

Using the exact effective tangent coordinate from
``ou3_p5_effective_vector_input``, ``y_T=H d_eff`` and
``eta_T=H(d_eff-c_perp)``.  Thus

    y^T S^-1 y - eta^T R^-1 eta
      = (H d_eff)^T S^-1(H d_eff)
        -(H(d_eff-c_perp))^T R^-1 H(d_eff-c_perp).       (2)

No radial ``eta`` budget remains anywhere in the information calculation.
Moreover, on a Cayley cell ``||c||<=q``, exact vector geometry gives

    ||d_eff-c_perp|| / ||c_perp|| <= q/sqrt(4+q^2),

and because ``H`` is ``||v||`` times an isometry on the tangent plane,

    ||H(d_eff-c_perp)||_{R^-1}^2 / ||H c||_{R^-1}^2
      <= q^2/(4+q^2) < 1.                                (3)

This is strictly below one for every finite Cayley radius, unlike the retired
standalone ``||eta||^2/||y||^2=q^2/4`` diagnostic, which can exceed one on the
widened first-S chart.  Equation (3) is still only a geometry factor: the final
positive word margin must use the source-correlated interval ``P,H,R,K,S`` cell
to lower-bound the useful ``S^-1`` term in (2).  This module therefore narrows
the remaining numerical obligation but does not promote P5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_effective_vector_input as VEFF

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def tangent_penalty_ratio_upper(q_hi: float) -> float:
    """Outward upper of q^2/(4+q^2), strictly below one for finite q."""
    q = float(q_hi)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    q2 = up(q * q)
    den = down(4.0 + down(q * q))
    if not den > q2:
        # The exact denominator exceeds the exact numerator by four.  At very
        # large q binary64 subtraction could hide that fact; P5's finite chart
        # is far below that regime, so fail closed rather than return one.
        raise RuntimeError("cannot resolve strict finite-Cayley tangent ratio")
    return up(q2 / den)


def tangent_effective_gain_lower(q_hi: float) -> float:
    """Lower of 4/(4+q^2) for ||d_eff||/||c_perp||."""
    q = float(q_hi)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    den = up(4.0 + up(q * q))
    return down(4.0 / den)


def _cell(row: dict) -> dict:
    q_lo, q_hi = map(float, row["q_interval"])
    defect = float(row["mag_effective_vs_tangent_defect_ratio_upper"])
    penalty = tangent_penalty_ratio_upper(q_hi)
    # defect^2 is another outward representation of the same exact factor.
    defect_sq = up(defect * defect)
    return {
        "index": int(row["index"]),
        "q_interval": [q_lo, q_hi],
        "effective_tangent_gain_lower": tangent_effective_gain_lower(q_hi),
        "effective_tangent_gain_upper": float(row["mag_effective_tangent_coordinate_gain_upper"]),
        "effective_vs_linear_tangent_defect_norm_ratio_upper": defect,
        "effective_vs_linear_tangent_penalty_information_ratio_upper": penalty,
        "defect_ratio_squared_consistency_upper": defect_sq,
        "radial_Joseph_information_contribution_exact_zero": True,
        "strict_tangent_penalty_ratio_below_one": penalty < 1.0,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("mag-information reduction domain must not be trajectory fitted")

    veff = VEFF.build(domain_path)
    failures = [f"effective-vector-input: {x}" for x in VEFF.validate(veff)]
    mag = veff.get("magnetometer", {})
    if mag.get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("radial gain-null prerequisite missing")
    if veff.get("source_semantics", {}).get("configured_magnetometer_covariance_isotropic") is not True:
        failures.append("configured magnetometer isotropy prerequisite missing")

    cells = [_cell(row) for row in veff.get("annular_effective_input_cells", [])]
    if not cells:
        failures.append("mag-information subdivision is empty")
    for row in cells:
        if row["radial_Joseph_information_contribution_exact_zero"] is not True:
            failures.append("radial Joseph cancellation lost")
            break
        if row["strict_tangent_penalty_ratio_below_one"] is not True:
            failures.append("finite-cell tangent penalty is not strict")
            break
        if not (0.0 < float(row["effective_tangent_gain_lower"])
                <= float(row["effective_tangent_gain_upper"]) <= 1.0):
            failures.append("effective tangent gain enclosure invalid")
            break
        # The two outward forms need not be bit-identical, but the direct
        # information ratio must enclose the squared norm-ratio representation.
        if float(row["effective_vs_linear_tangent_penalty_information_ratio_upper"]) + 1e-15 \
                < float(row["defect_ratio_squared_consistency_upper"]):
            failures.append("tangent penalty ratio does not enclose defect squared")
            break

    qmax = float(cells[-1]["q_interval"][1]) if cells else math.nan
    max_penalty = max((float(x["effective_vs_linear_tangent_penalty_information_ratio_upper"])
                       for x in cells), default=math.inf)
    min_gain = min((float(x["effective_tangent_gain_lower"]) for x in cells), default=0.0)
    if cells and not (math.isfinite(qmax) and qmax > 0.0 and 0.0 <= max_penalty < 1.0 and min_gain > 0.0):
        failures.append("widened magnetometer information factors invalid")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_MAGNETOMETER_JOSEPH_RADIAL_CANCELLATION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "configured_R_isotropic": True,
        "radial_subspace_invariant_under_S_and_R": True,
        "radial_residual_equals_radial_eta": True,
        "radial_Joseph_information_cancels_exactly": True,
        "exact_reduced_Joseph_identity": (
            "y^T S^-1 y-eta^T R^-1 eta=(H d_eff)^T S^-1(H d_eff)"
            "-(H(d_eff-c_perp))^T R^-1 H(d_eff-c_perp)"
        ),
        "tangent_penalty_ratio_identity": (
            "||H(d_eff-c_perp)||_R^-1^2/||H c||_R^-1^2 <= q^2/(4+q^2) < 1"
        ),
        "standalone_radial_eta_information_budget_used": False,
        "standalone_full_vector_eta_information_budget_used": False,
        "useful_S_inverse_term_still_requires_source_correlated_interval_cell": True,
        "annular_information_cells": cells,
        "subdivision_cell_count": len(cells),
        "widened_q_upper": qmax,
        "widened_tangent_penalty_ratio_upper": max_penalty,
        "widened_effective_tangent_gain_lower": min_gain,
        "complete_word_numerical_certificate_closed_here": False,
        "P5_MAGNETOMETER_INFORMATION_REDUCTION_CERTIFICATE": "PASS" if not failures else "FAIL",
        "next_obligation": (
            "combine each strict tangent-geometry cell with the jointly reachable interval P,H,R,K,S cell and signed quaternion/reset cell over every later prefix; only the tangent deformation penalty remains"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "configured_R_isotropic",
                "radial_subspace_invariant_under_S_and_R", "radial_residual_equals_radial_eta",
                "radial_Joseph_information_cancels_exactly",
                "useful_S_inverse_term_still_requires_source_correlated_interval_cell"):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed",
                "standalone_radial_eta_information_budget_used",
                "standalone_full_vector_eta_information_budget_used",
                "complete_word_numerical_certificate_closed_here"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    cells = d.get("annular_information_cells", [])
    if not cells or len(cells) != d.get("subdivision_cell_count"):
        failures.append("magnetometer information cells missing")
    else:
        for row in cells:
            if row.get("radial_Joseph_information_contribution_exact_zero") is not True:
                failures.append("cell lost radial Joseph cancellation")
                break
            p = float(row.get("effective_vs_linear_tangent_penalty_information_ratio_upper", math.inf))
            if not 0.0 <= p < 1.0:
                failures.append("cell tangent penalty ratio is not strict")
                break
    if not (0.0 <= float(d.get("widened_tangent_penalty_ratio_upper", math.inf)) < 1.0):
        failures.append("widened tangent penalty ratio is not strict")
    if not float(d.get("widened_effective_tangent_gain_lower", 0.0)) > 0.0:
        failures.append("widened effective tangent gain lower is not positive")
    if not failures and d.get("P5_MAGNETOMETER_INFORMATION_REDUCTION_CERTIFICATE") != "PASS":
        failures.append("magnetometer information reduction did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_MAGNETOMETER_INFORMATION_REDUCTION_CERTIFICATE"],
        "q_upper": out["widened_q_upper"],
        "tangent_penalty_ratio_upper": out["widened_tangent_penalty_ratio_upper"],
        "effective_tangent_gain_lower": out["widened_effective_tangent_gain_lower"],
        "cells": out["subdivision_cell_count"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
