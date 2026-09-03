#!/usr/bin/env python3
"""Exact magnetometer radial cancellation in the OU-III Joseph identity.

For the shipping magnetometer update on the declared proof branch,

    H = -[v]_x,          R = r I,
    S = H P H' + r I,

where v is the predicted body-frame magnetic reference.  Since H'v=0,

    P H' v = 0,          S v = r v,          S^-1 v = v/r.

Write the exact finite-angle residual as y=t+a v with t perpendicular to v.
The tangent model h=H c is also perpendicular to v.  Hence eta=y-h has radial
part a v and

    y' S^-1 y - eta' R^-1 eta
      = t' S^-1 t - (t-h)' R^-1 (t-h).

The radial finite-angle residual therefore cancels *exactly in Joseph energy*;
it is not merely a zero-gain state-correction direction.

Moreover

    d = H' y / ||v||^2

satisfies H d=t, so the retained signed contribution is the tangent directional
form

    (H d)' S^-1 (H d) - (H(d-c))' R^-1 H(d-c).

For c=q u the exact Cayley formula gives

    d = A c_perp - B alpha (c_perp x v_hat),
    A=4/(4+q^2), B=2/(4+q^2), alpha=c'v_hat.

The two terms are orthogonal, so over ||c||<=q_o

    A_o ||c_perp|| <= ||d|| <= 2/sqrt(4+q_o^2) ||c_perp||,
    ||d-c_perp|| <= q_o/sqrt(4+q_o^2) ||c_perp||,

with A_o=4/(4+q_o^2).  The lower gain is what a directional P4 accumulator
needs to compare exact finite-angle tangent information to the P3 linear packet.

This primitive is source-bound and valid over the retained 0.8-rad Cayley
sector.  It supplies an exact directional operation reduction; it does not
scalarize a packet or promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SIM = REPO / "src" / "util" / "W3dSimCommon.h"
CERT_SIM = REPO / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp"
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _source_contract() -> list[str]:
    mekf = MEKF.read_text(encoding="utf-8")
    sim = SIM.read_text(encoding="utf-8")
    cert = CERT_SIM.read_text(encoding="utf-8")
    checks = (
        (mekf, "const Matrix3 J_att = -skew_symmetric_matrix(v2hat);", "mag H=-[v]x"),
        (mekf, "S_mat.noalias() += J_att * P_th_th * J_att.transpose();", "mag S=HPH'+R"),
        (mekf, "PCt.noalias() += Pext.template block<NX,3>(0,0) * J_att.transpose();", "mag PH'"),
        (sim, "const Vector3f sigma_m(sigma_m_uT, sigma_m_uT, sigma_m_uT);", "isotropic mag sigma"),
        (cert, "cfg.sigma_m = sigma_m * kSigmaMRescale;", "certificate mag rescale preserves isotropy"),
    )
    return [f"shipping source marker changed: {label}" for text, marker, label in checks if marker not in text]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("magnetometer radial Joseph proof must not be trajectory fitted")

    cayley = CAYLEY.build(path)
    vector = VECTOR.build()
    failures = [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    failures += _source_contract()

    q = float(cayley["cayley_radius_upper"])
    if not (math.isfinite(q) and 0.0 < q < 1.0):
        failures.append("invalid 0.8-rad Cayley radius")
    q2_hi = up(q*q) if math.isfinite(q) else math.inf
    den_lo = down(math.sqrt(down(4.0 + q2_hi))) if q > 0.0 else 2.0
    if not den_lo > 0.0:
        failures.append("effective tangent denominator lost positivity")
        gain_lo = 0.0
        gain_hi = math.inf
        defect = math.inf
    else:
        gain_lo = down(4.0 / up(4.0 + q2_hi))
        gain_hi = min(1.0, up(2.0 / den_lo))
        defect = up(q / den_lo)

    if not (0.0 < gain_lo <= gain_hi <= 1.0):
        failures.append("effective tangent coordinate gain bounds invalid")
    if not (0.0 <= defect < 1.0):
        failures.append("effective tangent defect ratio reached one")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_MAGNETOMETER_RADIAL_JOSEPH_CANCELLATION",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "shipping_H_theta": "H_m=-[v]_x",
        "configured_R_isotropic": True,
        "H_transpose_v_exact_zero": True,
        "PH_transpose_v_exact_zero": True,
        "S_action_on_radial_vector": "S_m v=r_m v",
        "S_inverse_action_on_radial_vector": "S_m^-1 v=v/r_m",
        "kalman_gain_radial_action_exact_zero": True,
        "radial_Joseph_positive_term_equals_R_inverse_term": True,
        "radial_Joseph_energy_cancellation_exact": True,
        "retained_tangent_signed_identity": (
            "D_m=(H_m d_m)^T S_m^-1(H_m d_m) "
            "-(H_m(d_m-c))^T R_m^-1 H_m(d_m-c)"
        ),
        "effective_tangent_coordinate": "d_m=H_m^T y_m/||v||^2",
        "effective_coordinate_formula": (
            "d_m=A c_perp-B alpha(c_perp x v_hat), A=4/(4+q^2), B=2/(4+q^2)"
        ),
        "effective_tangent_coordinate_gain_lower": gain_lo,
        "effective_tangent_coordinate_gain_upper": gain_hi,
        "effective_vs_linear_tangent_defect_ratio_upper": defect,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "outer_sector_covered": float(cayley["outer_angle_rad"]) >= 0.80,
        "standalone_radial_eta_penalty_used": False,
        "directional_form_retained_until_word_scalarization": True,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_EXACT_MAGNETOMETER_RADIAL_JOSEPH_CANCELLATION":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "configured_R_isotropic",
        "H_transpose_v_exact_zero", "PH_transpose_v_exact_zero",
        "kalman_gain_radial_action_exact_zero",
        "radial_Joseph_positive_term_equals_R_inverse_term",
        "radial_Joseph_energy_cancellation_exact", "outer_sector_covered",
        "directional_form_retained_until_word_scalarization",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "standalone_radial_eta_penalty_used", "complete_H18_A21_word_established_here",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    lo = d.get("effective_tangent_coordinate_gain_lower")
    hi = d.get("effective_tangent_coordinate_gain_upper")
    defect = d.get("effective_vs_linear_tangent_defect_ratio_upper")
    if not isinstance(lo, (int, float)) or isinstance(lo, bool) or not (0.0 < float(lo) <= 1.0):
        f.append("invalid effective tangent gain lower")
    if not isinstance(hi, (int, float)) or isinstance(hi, bool) or not (0.0 < float(hi) <= 1.0):
        f.append("invalid effective tangent gain upper")
    if isinstance(lo, (int, float)) and isinstance(hi, (int, float)) and float(lo) > float(hi):
        f.append("effective tangent gain interval inverted")
    if not isinstance(defect, (int, float)) or isinstance(defect, bool) or not (0.0 <= float(defect) < 1.0):
        f.append("invalid effective tangent defect ratio")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("magnetometer radial cancellation not attached to 0.8-rad sector")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "outer_angle_rad": d["outer_angle_rad"],
        "effective_gain_lower": d["effective_tangent_coordinate_gain_lower"],
        "effective_gain_upper": d["effective_tangent_coordinate_gain_upper"],
        "effective_defect_ratio_upper": d["effective_vs_linear_tangent_defect_ratio_upper"],
        "radial_Joseph_cancelled": d["radial_Joseph_energy_cancellation_exact"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
