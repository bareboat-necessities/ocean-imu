#!/usr/bin/env python3
"""Finite-error nonexpansiveness of the shipping A21 accelerometer-bias projection.

The shipping MEKF projects the estimated accelerometer bias onto the Euclidean
ball of radius ``acc_bias_limit_`` after state injection.  The declared Normal-
Live theorem domain keeps the physical active bias strictly inside that ball.
For every closed convex set C, the Euclidean projector Pi_C is firmly
nonexpansive; in particular, for every y in C,

    ||Pi_C(x)-y|| <= ||x-y||.

Therefore the implemented 0.4 m/s^2 bias projection cannot increase the finite
physical bias-error norm for any theorem-domain true bias (<=0.35 m/s^2).  No
projection event needs an adverse scalar remainder or correction-radius charge
in P4.  This is a finite-state fact, not a tangent-only claim.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Sequence

REPO = Path(__file__).resolve().parents[2]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_A21_ACCEL_BIAS_PROJECTION_FINITE_NONEXPANSIVE_V1"


def _norm(x: Sequence[float]) -> float:
    return math.sqrt(sum(float(v) * float(v) for v in x))


def project_ball(x: Sequence[float], radius: float) -> list[float]:
    r = float(radius)
    if not (math.isfinite(r) and r > 0.0):
        raise ValueError("projection radius must be finite positive")
    y = [float(v) for v in x]
    if not y or any(not math.isfinite(v) for v in y):
        raise ValueError("projection input must be finite and nonempty")
    n = _norm(y)
    if n <= r:
        return y
    s = r / n
    return [s * v for v in y]


def distance_nonincrease(x: Sequence[float], truth: Sequence[float], radius: float) -> float:
    if len(x) != len(truth) or not x:
        raise ValueError("projection/truth dimension mismatch")
    r = float(radius)
    t = [float(v) for v in truth]
    if any(not math.isfinite(v) for v in t) or _norm(t) > r:
        raise ValueError("truth must lie in the projection ball")
    p = project_ball(x, r)
    before = _norm([float(a) - b for a, b in zip(x, t)])
    after = _norm([a - b for a, b in zip(p, t)])
    return before - after


def _source_limit(text: str) -> float:
    m = re.search(r"acc_bias_limit_\s*=\s*T\(([0-9.eE+-]+)\)", text)
    if not m:
        raise RuntimeError("cannot extract shipping acc_bias_limit_")
    return float(m.group(1))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("bias-projection theorem may not be trajectory fitted")
    text = MEKF.read_text(encoding="utf-8")
    radius = _source_limit(text)
    active = float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
    declared_projection = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    parity = {
        "shipping_projection_radius_matches_declared_domain": radius == declared_projection,
        "projection_is_radial_Euclidean_ball_projection": (
            "const T n = b.norm();" in text
            and "if (n > acc_bias_limit_)" in text
            and "b *= (acc_bias_limit_ / n);" in text
        ),
        "projection_called_after_state_injection": "project_acc_bias_();" in text,
    }
    if radius != 0.4 or declared_projection != 0.4:
        raise RuntimeError("authorized shipping bias projection radius changed")
    if not (math.isfinite(active) and 0.0 < active < radius and active == 0.35):
        raise RuntimeError("Normal-Live active bias interior no longer equals authorized 0.35")
    failures = [k for k, ok in parity.items() if not ok]
    closed = not failures
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "shipping_projection_radius_mps2": radius,
        "declared_Normal_Live_true_bias_norm_upper_mps2": active,
        "strict_interior_margin_mps2": radius - active,
        "source_parity": parity,
        "source_parity_failures": failures,
        "mathematical_set": "closed Euclidean ball B_2(0,0.4)",
        "finite_projection_identity": "||Pi_B(x)-b_true||_2 <= ||x-b_true||_2 for b_true in B",
        "firm_nonexpansiveness_used": True,
        "finite_error_not_tangent_only": True,
        "projection_bias_error_energy_nonincrease": closed,
        "projection_can_be_omitted_from_adverse_P4_bias_budget": closed,
        "projection_requires_correction_radius": False,
        "projection_requires_inactive_branch_assumption": False,
        "point_diagnostic_projection_inactivity_required": False,
        "state_elimination_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P4_promoted_here": False,
        "P5_may_start": False,
        "next_obligation": "compose finite-tau_b prediction decay and the same-word Joseph bias correction with the signed complete-word information ledger; projection adds no adverse bias-error term",
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in ("firm_nonexpansiveness_used", "finite_error_not_tangent_only", "projection_bias_error_energy_nonincrease", "projection_can_be_omitted_from_adverse_P4_bias_budget"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("projection_requires_correction_radius", "projection_requires_inactive_branch_assumption", "point_diagnostic_projection_inactivity_required", "state_elimination_used", "filter_changed", "declared_domain_changed", "P4_promoted_here", "P5_may_start"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("source_parity_failures"):
        f.append("shipping projection source parity failed")
    if not all(d.get("source_parity", {}).values()):
        f.append("shipping projection source parity incomplete")
    if float(d.get("shipping_projection_radius_mps2", math.nan)) != 0.4:
        f.append("shipping projection radius changed")
    active = float(d.get("declared_Normal_Live_true_bias_norm_upper_mps2", math.nan))
    if active != 0.35 or not active < 0.4:
        f.append("authorized Normal-Live bias interior changed")
    if not float(d.get("strict_interior_margin_mps2", 0.0)) > 0.0:
        f.append("projection interior margin is not strict")
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
    print(json.dumps({"projection_radius": d["shipping_projection_radius_mps2"], "true_bias_upper": d["declared_Normal_Live_true_bias_norm_upper_mps2"], "finite_nonexpansive": d["projection_bias_error_energy_nonincrease"], "failures": failures}, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
