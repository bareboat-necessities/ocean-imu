#!/usr/bin/env python3
"""Source-only design probe for a usable OU-III P3/P4 superword.

The old quantitative vector-UCO constant uses the two consecutive magnetic
packets *inside one PE event* to expose gyro bias.  Their separation is only
about 40 ms, so the analytical inverse pays roughly 1/Delta^2 even though the
deployment theorem already assumes a PE event in every sliding one-second
window.

This probe uses that recurrence instead of throwing it away.  In a three-second
superword choose one PE event from [0,T] and one from [a T,(a+1)T], with a=1.9
and T the declared recurrence window.  Both windows are guaranteed by exactly
the existing theorem hypothesis.  Accounting conservatively for the finite
packet span gives a positive event separation near 0.86 s while the largest
separation stays below three seconds.  The same rate-transport perturbation
lemma used by ``ou3_vector_uco_certificate`` then yields a much stronger
attitude/gyro-bias information lower bound.

No deployment assumption is narrowed.  This is initially a design probe: it
re-evaluates the reconstructed P3 limiting source corner with the stronger
source-implied UCO number and reports the resulting covariance upper/margin.
Promotion is left false until the superword language and complete H/A matrix
backend consume the result directly.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3_direct as P3D
import ou3_p4_worst_translation_cell as WORST

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
WINDOW_OFFSET = 1.9
SUPERWORD_WINDOWS = 3


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _recurrence_alpha6(live: dict, vector: dict) -> dict:
    B = P3D.BASE
    vc = vector["configured_measurement_bounds"]
    f = float(live["specific_force_norm_lower_mps2"])
    m = float(live["magnetic_vector_norm_lower_uT"])
    s = float(live["vector_sine_separation_lower"])
    rate = float(live["body_rate_norm_upper_deg_s"]) * math.pi / 180.0
    T = float(live["vector_pe_recurrence_window_s"])
    packet_span = float(vector["operating_envelope"]["packet_gap_s"][1])
    ra = up(float(vc["acc_measurement_variance_upper"]))
    rm = up(float(vc["mag_measurement_variance_upper"]))
    angular = down(s*s / up(1.0 + math.sqrt(max(0.0, 1.0-s*s))))
    mu = down(min(f*f/ra, m*m/rm) * angular)

    # Event 1 lies in [0,T], event 2 in [aT,(a+1)T].  The packet reference
    # within each certified pair may shift by one packet span, so debit/add it
    # on the separation endpoints rather than assuming coincident references.
    dmin = down((WINDOW_OFFSET - 1.0) * T - packet_span)
    dmax = up((WINDOW_OFFSET + 1.0) * T + packet_span)
    if not dmin > 0.0:
        raise RuntimeError("recurrence event windows do not have positive separation")
    omega_dmax = up(rate * dmax)
    bracket = down(1.0 - up(0.5 * omega_dmax))
    if not bracket > 0.0:
        raise RuntimeError("recurrence-separated gyro transport perturbation lost positivity")
    gamma = down(dmin * bracket / float(vector["operating_envelope"]["gyro_bias_time_scale_s"]))
    alpha = down(mu / up(1.0 + up(2.0 / down(gamma*gamma))))
    old = float(B.vector_alpha6(live, vector))
    return {
        "recurrence_window_s": T,
        "superword_horizon_s": SUPERWORD_WINDOWS * T,
        "window_offset_multiple": WINDOW_OFFSET,
        "packet_span_upper_s": packet_span,
        "selected_event_separation_lower_s": dmin,
        "selected_event_separation_upper_s": dmax,
        "omega_times_separation_upper": omega_dmax,
        "transport_bracket_lower": bracket,
        "mu_theta_lower": mu,
        "scaled_gamma_lower": gamma,
        "old_two_consecutive_packet_alpha6_lower": old,
        "recurrence_separated_alpha6_lower": alpha,
        "alpha6_widening_factor_lower": down(alpha / old),
    }


def _limiting_cell(mode: str, path: Path, alpha6: float) -> dict:
    c = WORST.build_cell(mode, path)
    row = P3D.mode_cell(
        mode, c["x"], c["rho_translation_lower"], c["sigma"], c["rs"],
        c["live"], c["vector"], c["process"], c["sched"], alpha6,
    )
    old = c["row"]
    return {
        "source_cell": WORST.serializable(c),
        "old_attitude_covariance_upper": float(old["Sigma_diagonal_upper"][0]),
        "new_attitude_covariance_upper": float(row["Sigma_diagonal_upper"][0]),
        "old_gyro_bias_covariance_upper": float(old["Sigma_diagonal_upper"][3]),
        "new_gyro_bias_covariance_upper": float(row["Sigma_diagonal_upper"][3]),
        "old_full_margin_lower": float(old["relative_Riccati_injection_margin_lower"]),
        "new_full_margin_lower": float(row["relative_Riccati_injection_margin_lower"]),
        "old_nontranslation_margin_lower": float(old["direct_nontranslation_margin_lower"]),
        "new_nontranslation_margin_lower": float(row["direct_nontranslation_margin_lower"]),
        "attitude_covariance_tightening_factor_lower": down(
            float(old["Sigma_diagonal_upper"][0]) / float(row["Sigma_diagonal_upper"][0])
        ),
        "margin_widening_factor_lower": down(
            float(row["relative_Riccati_injection_margin_lower"]) /
            float(old["relative_Riccati_injection_margin_lower"])
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("superword design must not be trajectory fitted")
    live = domain["normal_live"]
    vector = P3D.BASE.VECTOR.build()
    failures = [f"vector: {x}" for x in P3D.BASE.VECTOR.validate(vector)]
    recurrence = None
    modes = {}
    if not failures:
        try:
            recurrence = _recurrence_alpha6(live, vector)
            if recurrence["superword_horizon_s"] < recurrence["selected_event_separation_upper_s"]:
                failures.append("three-window superword does not contain selected PE events")
            for mode in ("H", "A"):
                modes[mode] = _limiting_cell(
                    mode, path, float(recurrence["recurrence_separated_alpha6_lower"])
                )
        except Exception as exc:
            failures.append(str(exc))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_RECURRENCE_SEPARATED_SUPERWORD_DESIGN",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "deployment_domain_changed": False,
        "new_PE_hypothesis_added": False,
        "existing_sliding_recurrence_hypothesis_reused": True,
        "superword_is_concatenation_of_existing_source_complete_windows": True,
        "recurrence_geometry": recurrence,
        "modes": modes,
        "theorem_promotion": "DESIGN_ONLY",
        "P3_USABLE_SUPERWORD_PROMOTED": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "existing_sliding_recurrence_hypothesis_reused",
        "superword_is_concatenation_of_existing_source_complete_windows",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployment_domain_changed",
        "new_PE_hypothesis_added", "P3_USABLE_SUPERWORD_PROMOTED",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    r = d.get("recurrence_geometry") or {}
    if not float(r.get("recurrence_separated_alpha6_lower", 0.0)) > float(r.get("old_two_consecutive_packet_alpha6_lower", math.inf)):
        f.append("recurrence separation did not improve alpha6")
    if not float(r.get("omega_times_separation_upper", math.inf)) < 2.0:
        f.append("recurrence separation left gyro transport lemma")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if not float(m.get("new_attitude_covariance_upper", math.inf)) < float(m.get("old_attitude_covariance_upper", 0.0)):
            f.append(f"{mode}: attitude covariance upper did not tighten")
        if not float(m.get("new_full_margin_lower", 0.0)) >= float(m.get("old_full_margin_lower", math.inf)):
            f.append(f"{mode}: P3 margin regressed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "recurrence": d.get("recurrence_geometry"),
        "modes": d.get("modes"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
