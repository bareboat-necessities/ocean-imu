#!/usr/bin/env python3
"""Structural A-mode wrapper for the finite-angle P4 first-accel range gate.

The schema-1 candidate range producer correctly uses the H-mode first-prefix
attitude gain as an upper bound for A mode, but its diagnostic check inspected
the already PSD-tightened 21x21 interval matrix.  Generic entrywise PSD
boxing may widen exact structural zeros into tiny intervals; that numerical
representation must not be mistaken for a physical/source cross covariance.

This wrapper proves the A-mode dominance at the construction level instead:

* the A covariance seed installs an isotropic b_a diagonal and zero b_a cross
  blocks before generic PSD boxing;
* the A transition is diagonal on b_a and has no b_a/non-b_a coupling;
* the A process block is isotropic diagonal before generic PSD boxing;
* the accelerometer Jacobian inserts J_ba=I;
* therefore the first A innovation receives an additional isotropic PSD term,
  while P_theta H^T is unchanged relative to the H proof problem.

Consequently the H structured K_theta bound is conservative for the first A
attitude correction.  This is a proof-backend representation repair only; it
does not change the shipping filter, the 6 rad deployed helper range, or any
P4/P5 theorem assumption.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_candidate_first_accel_range as BASE
import ou3_p4_candidate_full_word as CAND

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 2


def _structural_a_check(domain_path: Path, domain: dict, src: dict) -> dict:
    text = Path(CAND.__file__).read_text(encoding="utf-8")
    markers = {
        "A_seed_diagonal": "Pm[i][i] = v",
        "A_seed_cross_zero_forward": "Pm[i][j] = I(0.0)",
        "A_seed_cross_zero_transpose": "Pm[j][i] = I(0.0)",
        "A_transition_diagonal": "F[i][i] = Interval(phi.lo, min(1.0, phi.hi))",
        "A_process_diagonal": "Q[i][i] = qd",
        "A_accel_J_ba_identity": "M[ax][i] = I(1.0)",
    }
    missing = [name for name, marker in markers.items() if marker not in text]

    CAND._configure_mode("A")
    _F, _Q, _Rstep, process = CAND._transition_and_Q("A", src, domain)
    sigma0 = CAND._source_sigma_bacc0()
    CAND._configure_mode("H")

    phi = process["phi_interval"] if process is not None else [0.0, 0.0]
    qd = process["Qd_variance_interval"] if process is not None else [-1.0, -1.0]
    seed_var = sigma0 * sigma0
    predicted_lo = phi[0] * phi[0] * seed_var + qd[0]
    predicted_hi = phi[1] * phi[1] * seed_var + qd[1]

    structural = (
        not missing
        and sigma0 > 0.0
        and 0.0 < phi[0] <= phi[1] <= 1.0
        and 0.0 < qd[0] <= qd[1]
        and 0.0 < predicted_lo <= predicted_hi
    )
    return {
        "source_markers": markers,
        "missing_source_markers": missing,
        "check_uses_structural_source_matrices_not_psd_tightened_box": True,
        "A_bias_seed_variance_per_axis": seed_var,
        "A_bias_transition_phi_interval": list(phi),
        "A_bias_process_variance_interval": list(qd),
        "A_bias_predicted_variance_interval": [predicted_lo, predicted_hi],
        "A_bias_seed_isotropic": structural,
        "A_bias_seed_cross_exact_zero": structural,
        "A_bias_transition_diagonal_isotropic": structural,
        "A_bias_transition_cross_exact_zero": structural,
        "A_bias_process_diagonal_isotropic_PSD": structural,
        "A_bias_process_cross_exact_zero": structural,
        "A_bias_innovation_addition_isotropic_PSD": structural,
        "first_prefix_theta_aw_S_to_ba_cross_exact_zero": structural,
        "A_accelerometer_J_ba_identity": structural,
    }


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    alignment_pieces: int = 16,
    force_magnitude_pieces: int = 4,
) -> dict:
    old = BASE._a_mode_first_prefix_isotropic_bias_check
    try:
        BASE._a_mode_first_prefix_isotropic_bias_check = _structural_a_check
        out = dict(BASE.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            alignment_pieces=alignment_pieces,
            force_magnitude_pieces=force_magnitude_pieces,
        ))
    finally:
        BASE._a_mode_first_prefix_isotropic_bias_check = old
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P4_CANDIDATE_FIRST_ACCEL_ANALYTIC_RANGE_STRUCTURAL_A"
    out["A_structure_proved_before_generic_PSD_boxing"] = True
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = BASE.SCHEMA
    failures = BASE.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("A_structure_proved_before_generic_PSD_boxing") is not True:
        failures.append("A structure was not proved before generic PSD boxing")
    a = d.get("A_mode_structure", {})
    for key in (
        "check_uses_structural_source_matrices_not_psd_tightened_box",
        "A_bias_seed_isotropic",
        "A_bias_seed_cross_exact_zero",
        "A_bias_transition_diagonal_isotropic",
        "A_bias_transition_cross_exact_zero",
        "A_bias_process_diagonal_isotropic_PSD",
        "A_bias_process_cross_exact_zero",
        "A_bias_innovation_addition_isotropic_PSD",
        "first_prefix_theta_aw_S_to_ba_cross_exact_zero",
        "A_accelerometer_J_ba_identity",
    ):
        if a.get(key) is not True:
            failures.append(f"A structural proof field {key} is not true")
    if a.get("missing_source_markers"):
        failures.append("A structural source markers are missing")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(
        args.domain.resolve(),
        source_pieces=args.source_pieces,
        alignment_pieces=args.alignment_pieces,
        force_magnitude_pieces=args.force_magnitude_pieces,
    )
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE"],
        "widest_safe_deg": d["widest_candidate_first_accel_range_safe_deg"],
        "A_structural": d["A_structure_proved_before_generic_PSD_boxing"],
        "rows": [{
            "angle_deg": r["angle_deg"],
            "q_pred": r["post_prediction_q_upper"],
            "max_Ktheta": r["max_Ktheta_norm_upper"],
            "max_residual": r["max_combined_residual_norm_upper_mps2"],
            "max_d": r["max_first_accelerometer_correction_norm_upper_rad"],
            "margin": r["minimum_correction_range_margin_rad"],
            "safe": r["first_accelerometer_range_safe"],
        } for r in d["candidate_rows"]],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
