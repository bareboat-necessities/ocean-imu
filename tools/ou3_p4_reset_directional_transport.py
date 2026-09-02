#!/usr/bin/env python3
"""Exact directional-information transport through the OU-III left-error reset.

After every accepted correction the shipping filter applies

    G(d) = I + 0.5 [d]x,
    P_r  = G_e P_+ G_e^T,

where ``G_e=diag(G,I)`` on the fixed-dimensional H/A state.  This reset is a
coordinate congruence, not an independent disturbance.  For the homogeneous
post-correction tangent state ``t`` and ``z_r=G_e t``,

    z_r^T P_r^-1 z_r = t^T P_+^-1 t.                    (1)

Likewise any PSD directional form ``Q`` transports as

    Q_r = G_e^-T Q G_e^-1,                              (2)

which preserves rank and nullity exactly.  No scalar covariance condition
number belongs in either identity.

The 3x3 reset block has the exact Gram identity

    G^T G = I + 0.25 (||d||^2 I - d d^T),

so its squared singular values are ``1, 1+||d||^2/4, 1+||d||^2/4``.  Hence G
is invertible for every finite correction, ``||G^-1||_2=1``, and the extended
reset cannot destroy a certified directional rank.

This module certifies only the homogeneous reset transport.  The exact deployed
quaternion injection differs from ``G_e t`` by a finite-angle Cayley reset
defect; that defect remains an explicit downstream P4 nonlinear obligation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST

REPO = Path(__file__).resolve().parents[1]
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def reset_spectral_facts(correction_norm: float) -> dict:
    """Return outward scalar spectral facts for a fixed correction-norm bound."""
    d = float(correction_norm)
    if not (math.isfinite(d) and d >= 0.0):
        raise ValueError("finite nonnegative correction norm required")
    transverse = up(1.0 + up(0.25 * up(d * d)))
    return {
        "correction_norm_upper": d,
        "reset_Gram_lambda_min_lower": 1.0,
        "reset_Gram_lambda_max_upper": transverse,
        "reset_min_singular_value_lower": 1.0,
        "reset_max_singular_value_upper": up(math.sqrt(transverse)),
        "reset_inverse_operator_norm_upper": 1.0,
        "reset_determinant_lower": 1.0,
        "reset_determinant_upper": transverse,
    }


def _source_binding_failures() -> list[str]:
    text = CORE.read_text(encoding="utf-8")
    required = {
        "reset_helper": "apply_left_error_reset",
        "reset_matrix": "Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta)",
        "attitude_covariance_congruence": "covariance(i,j) = sum;",
        "attitude_cross_covariance_left_action": "const T new0 = G(0,0)*old0",
    }
    return [name for name, marker in required.items() if marker not in text]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("reset transport must not be trajectory fitted")

    manifest = MANIFEST.build()
    failures = [f"manifest: {x}" for x in MANIFEST.validate(manifest)]
    failures += [f"missing source semantic: {x}" for x in _source_binding_failures()]
    resets = manifest.get("same_sample_reset_policy", {})
    if resets.get("single_shared_end_of_sample_reset") is not False:
        failures.append("source manifest merged sequential resets")

    modes = {}
    for mode, n in (("H", 18), ("A", 21)):
        modes[mode] = {
            "dimension": n,
            "extended_reset": "G_e=diag(I+0.5[d]x,I_{n-3})",
            "homogeneous_information_quadratic_gain_exact": 1.0,
            "PSD_directional_form_transport": "Q_r=G_e^-T Q G_e^-1",
            "directional_rank_preserved_exactly": True,
            "directional_nullity_preserved_exactly": True,
            "condition_number_multiplier_required": False,
            "nonlinear_Cayley_injection_defect_closed_here": False,
            "P4_PROMOTED": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_RESET_DIRECTIONAL_CONGRUENCE_TRANSPORT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_reset_matrix": "G=I+0.5[d]x",
        "reset_Gram_identity": "G^T G=I+0.25*(||d||^2 I-d d^T)",
        "reset_squared_singular_values": ["1", "1+||d||^2/4", "1+||d||^2/4"],
        "reset_invertible_for_every_finite_correction": True,
        "reset_inverse_operator_norm_exact": 1.0,
        "homogeneous_information_congruence_exact": True,
        "directional_form_congruence_exact": True,
        "condition_number_conversion_used": False,
        "modes": modes,
        "P4_RESET_DIRECTIONAL_TRANSPORT_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "carry each Joseph directional form through the exact reset congruence and source prediction to a common word endpoint; "
            "keep the deployed Cayley injection mismatch as the only reset nonlinear defect, then accumulate full-rank H/A word credit before scalarization"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "reset_invertible_for_every_finite_correction",
        "homogeneous_information_congruence_exact",
        "directional_form_congruence_exact",
        "P4_RESET_DIRECTIONAL_TRANSPORT_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "condition_number_conversion_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("reset_inverse_operator_norm_exact") != 1.0:
        f.append("reset inverse norm is not exact one")
    for mode, n in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode}: dimension mismatch")
        if m.get("homogeneous_information_quadratic_gain_exact") != 1.0:
            f.append(f"{mode}: information congruence gain changed")
        if m.get("directional_rank_preserved_exactly") is not True:
            f.append(f"{mode}: directional rank is not preserved")
        if m.get("directional_nullity_preserved_exactly") is not True:
            f.append(f"{mode}: directional nullity is not preserved")
        if m.get("condition_number_multiplier_required") is not False:
            f.append(f"{mode}: condition-number reset penalty reintroduced")
        if m.get("nonlinear_Cayley_injection_defect_closed_here") is not False:
            f.append(f"{mode}: reset defect prematurely closed")
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
        "inverse_norm": d["reset_inverse_operator_norm_exact"],
        "H_rank_preserved": d["modes"]["H"]["directional_rank_preserved_exactly"],
        "A_rank_preserved": d["modes"]["A"]["directional_rank_preserved_exactly"],
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
