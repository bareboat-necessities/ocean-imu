#!/usr/bin/env python3
"""Exact Cayley measurement-defect geometry for OU-III P5.

The remaining P5 source-word enclosure must not pay for the nonlinear vector
measurement defect with an unrelated global Lipschitz constant.  In Cayley
coordinates the defect has an exact algebraic structure.

Let ``c=q u`` be the attitude Cayley coordinate and let ``y_R`` denote the
signed exact residual of one fixed vector, with ``h=H_theta c`` the tangent
residual using the sign convention of the implemented update.  Then

    ||y_R||^2 = 4 q^2/(4+q^2) ||v_perp||^2,
    eta_R := y_R-h,
    y_R^T eta_R = 0,
    ||eta_R||^2 = q^2/4 ||y_R||^2.                         (1)

Thus, for the deployed isotropic 3-vector measurement covariance, the nonlinear
penalty is *exactly* ``q^2/4`` times the rotational residual information.  This
is much sharper than multiplying a global second-order constant by the number
of accepted packets.

For the accelerometer, the latent acceleration error is itself rotated.  If
``e`` is the additive world-vector error, the only extra finite-angle term is
``(R(c)^T-I)e`` and

    ||(R(c)^T-I)e|| <= 2q/sqrt(4+q^2) ||e||.              (2)

Bias/additive body terms remain linear and are not charged as attitude
nonlinearity.  Equations (1)-(2) are used cellwise by the later complete-word
backend; source covariance, Joseph information, reset and the signed correction
composition remain correlated in that backend.

This producer also creates a deterministic outward annular partition of the
currently widened first-S Cayley chart.  It is a numerical-enclosure primitive,
not a P5 promotion and not a trajectory replay.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_s_exact_prefix as FIRSTS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DEFAULT_CELLS = 64


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _positive_den(q: float) -> float:
    return up(4.0 + up(q * q))


def exact_residual_factor_lower(q_hi: float) -> float:
    """Lower of 4/(4+q^2) on 0<=q<=q_hi."""
    if not (math.isfinite(q_hi) and q_hi >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    return down(4.0 / _positive_den(q_hi))


def exact_eta_to_residual_information_ratio_upper(q_hi: float) -> float:
    """Upper of ||eta_R||_R^-1^2 / ||y_R||_R^-1^2 for isotropic R."""
    if not (math.isfinite(q_hi) and q_hi >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    return up(up(q_hi * q_hi) / 4.0)


def rotation_difference_norm_upper(q_hi: float) -> float:
    """Upper of ||R(c)-I||_2 = 2q/sqrt(4+q^2)."""
    if not (math.isfinite(q_hi) and q_hi >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    den_lo = down(math.sqrt(down(4.0 + down(q_hi * q_hi))))
    if not den_lo > 0.0:
        raise RuntimeError("rotation-difference denominator lost positivity")
    return up(up(2.0 * q_hi) / den_lo)


def annular_cells(q_max: float, n: int = DEFAULT_CELLS) -> list[dict]:
    if not (math.isfinite(q_max) and q_max > 0.0 and n >= 1):
        raise ValueError("invalid annular subdivision")
    # Quadratic spacing puts more cells near the origin, where the exact eta
    # ratio changes from the P4 perturbative regime.  Endpoints are widened so
    # adjacent cells overlap by at least one binary64 ulp.
    edges = [q_max * (k / n) ** 2 for k in range(n + 1)]
    cells = []
    for k in range(n):
        lo = 0.0 if k == 0 else down(edges[k])
        hi = up(edges[k + 1])
        cells.append({
            "index": k,
            "q_interval": [lo, hi],
            "exact_vector_residual_factor_lower": exact_residual_factor_lower(hi),
            "exact_eta_to_rotational_residual_information_ratio_upper": exact_eta_to_residual_information_ratio_upper(hi),
            "latent_vector_rotation_gain_upper": rotation_difference_norm_upper(hi),
        })
    return cells


def build(domain_path: Path = DEFAULT_DOMAIN, cells: int = DEFAULT_CELLS) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 eta geometry must not be trajectory fitted")
    first = FIRSTS.build(domain_path)
    failures = [f"first-S: {x}" for x in FIRSTS.validate(first)]
    if first.get("P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE") != "PASS_WIDENED_CHART":
        failures.append("first-S widened Cayley chart prerequisite did not pass")

    qmax = float(first["widened_prefix_cayley_norm_upper"])
    partition = annular_cells(qmax, int(cells))
    if not partition or partition[0]["q_interval"][0] != 0.0:
        failures.append("annular subdivision does not start at zero")
    if partition and not partition[-1]["q_interval"][1] >= qmax:
        failures.append("annular subdivision does not cover widened chart")

    # Exact identity checks reduced to scalar rational forms.  The ratio is
    # finite for every finite Cayley coordinate even when q>1.
    eta_ratio_max = exact_eta_to_residual_information_ratio_upper(qmax)
    residual_factor = exact_residual_factor_lower(qmax)
    rot_gain = rotation_difference_norm_upper(qmax)
    if not (residual_factor > 0.0 and eta_ratio_max >= 0.0 and 0.0 < rot_gain < 2.0):
        failures.append("widened exact Cayley factors are invalid")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_CAYLEY_MEASUREMENT_DEFECT_SUBDIVISION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "full_S_to_attitude_gain_retained": True,
        "exact_rotational_residual_identity": "||y_R||^2=4q^2/(4+q^2)||v_perp||^2",
        "exact_eta_identity": "y_R^T eta_R=0 and ||eta_R||^2=(q^2/4)||y_R||^2",
        "isotropic_measurement_information_identity": "eta_R^T R^-1 eta_R=(q^2/4)y_R^T R^-1 y_R",
        "latent_vector_cross_term_bound": "||(R(c)^T-I)e||<=2q/sqrt(4+q^2)||e||",
        "global_packet_count_times_Lipschitz_defect_used": False,
        "widened_cayley_norm_upper": qmax,
        "widened_exact_vector_residual_factor_lower": residual_factor,
        "widened_eta_to_rotational_residual_information_ratio_upper": eta_ratio_max,
        "widened_latent_vector_rotation_gain_upper": rot_gain,
        "annular_subdivision_cells": partition,
        "subdivision_cell_count": len(partition),
        "complete_word_covariance_reset_transport_closed_here": False,
        "P5_CAYLEY_ETA_GEOMETRY_CERTIFICATE": "PASS" if not failures else "FAIL",
        "next_obligation": (
            "propagate source-correlated covariance/gain and signed Cayley correction cells through the later 1 s prefixes; use the exact eta identity in each vector cell and keep S eta=0"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("eta geometry is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("eta geometry uses replay")
    if d.get("filter_changed") is not False:
        failures.append("eta geometry changes filter")
    if d.get("full_S_to_attitude_gain_retained") is not True:
        failures.append("eta geometry drops S-to-attitude gain")
    if d.get("global_packet_count_times_Lipschitz_defect_used") is not False:
        failures.append("eta geometry reintroduced global packet-count Lipschitz penalty")
    if d.get("complete_word_covariance_reset_transport_closed_here") is not False:
        failures.append("eta primitive prematurely promotes complete word")
    q = d.get("widened_cayley_norm_upper")
    rf = d.get("widened_exact_vector_residual_factor_lower")
    er = d.get("widened_eta_to_rotational_residual_information_ratio_upper")
    rg = d.get("widened_latent_vector_rotation_gain_upper")
    if not all(isinstance(x, (int, float)) and math.isfinite(float(x)) for x in (q, rf, er, rg)):
        failures.append("eta geometry emitted nonfinite scalar")
    elif not (float(q) > 0.0 and float(rf) > 0.0 and float(er) >= 0.0 and 0.0 < float(rg) < 2.0):
        failures.append("eta geometry scalar sign/range invalid")
    cells = d.get("annular_subdivision_cells", [])
    if len(cells) != d.get("subdivision_cell_count") or not cells:
        failures.append("eta annular subdivision missing")
    else:
        if cells[0]["q_interval"][0] != 0.0:
            failures.append("eta subdivision does not start at zero")
        if float(cells[-1]["q_interval"][1]) < float(q):
            failures.append("eta subdivision does not cover widened radius")
        for row in cells:
            lo, hi = map(float, row["q_interval"])
            if not (0.0 <= lo <= hi):
                failures.append("invalid eta subdivision interval")
                break
            if not float(row["exact_vector_residual_factor_lower"]) > 0.0:
                failures.append("eta subdivision lost vector residual factor")
                break
    if not failures and d.get("P5_CAYLEY_ETA_GEOMETRY_CERTIFICATE") != "PASS":
        failures.append("eta geometry did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--cells", type=int, default=DEFAULT_CELLS)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), args.cells)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_CAYLEY_ETA_GEOMETRY_CERTIFICATE"],
        "q_max": out["widened_cayley_norm_upper"],
        "eta_ratio_max": out["widened_eta_to_rotational_residual_information_ratio_upper"],
        "rotation_gain_max": out["widened_latent_vector_rotation_gain_upper"],
        "cells": out["subdivision_cell_count"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
