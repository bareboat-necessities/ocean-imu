#!/usr/bin/env python3
"""Augmented complete-word P4 design with exact S-selector cancellation.

The v1 augmented word differentiated covariance/gain/reset correctly but still
updated the three measured S states as ``z_S-dx_S``.  For the exact homogeneous
S=0 selector this destroys the same near-identity-gain cancellation diagnosed
in the non-augmented word.  This wrapper keeps all v1 Joseph/covariance AD and
finite-cell affine accounting, while replacing only the measured S-state value
and derivative by the algebraically identical

    H z_plus = R (H P H^T + R)^-1 r.

The innovation solution is already an AD quantity in v1, so this replacement
also preserves its full word-entry derivative.  Nonzero physical S_true remains
an external forcing contribution to the final affine b.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_augmented_complete_word_design as M
import ou3_p4_augmented_covariance_ad as AUG
import ou3_p4_joint_word_dissipation_design as D

SCHEMA = 2
_ORIGINAL_MEASUREMENT = M._measurement


def _measurement_v2(mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
                    residual, center_residual, eta, entry_center,
                    signed_word, beta_total, epsilon):
    # Reproduce v1 exactly so the signed form, affine beta, covariance AD and
    # reset derivatives remain on the same arithmetic path.
    cell = AUG.measurement_and_reset(
        Pad, Hm, Rm, residual, value_tighten=M.H._psd_tighten
    )
    ccell = AUG.measurement_and_reset(
        Pcenter, Hm, Rm, center_residual, value_tighten=M.H._psd_tighten
    )
    dform, beta, bmeta = M._finite_cell_terms(
        residual, eta, center_residual, entry_center,
        AUG.values(cell["Sinv"]), Rm, epsilon,
    )
    signed_word = M._add_form(signed_word, dform)
    beta_total = math.nextafter(beta_total + beta, math.inf) if beta > 0.0 else beta_total
    zout = M._state_update(z, cell["dx"])
    zcout = M._state_update(zcenter, ccell["dx"])

    selector_used = False
    if kind == "S":
        n = len(z)
        # The caller's S branch is the exact selector by construction.  Keep
        # this explicit rather than pattern-matching an arbitrary measurement.
        rplus = D._ad_matvec(Rm, cell["innovation_solution"])
        rcplus = D._ad_matvec(Rm, ccell["innovation_solution"])
        for i in range(3):
            zout[12 + i] = rplus[i]
            zcout[12 + i] = rcplus[i]
        selector_used = True

    meta = {
        "kind": kind,
        "inverse_backend": cell["inverse_meta"].get("inverse_backend"),
        "correction_theta_norm_upper": M.AD._norm_upper([q.val for q in cell["dx"][:3]]),
        "covariance_path_derivative_included": cell["covariance_path_derivative_included"],
        "reset_covariance_depends_on_correction_AD": cell["reset_covariance_depends_on_correction_AD"],
        "dependency_preserving_S_selector_identity_used": selector_used,
        "physical_S_true_forcing_must_enter_affine_b": bool(selector_used),
        "K_interval_matrix_materialized": False,
        **bmeta,
    }
    return (cell["P_accepted_post_reset"], zout,
            ccell["P_accepted_post_reset"], zcout,
            signed_word, beta_total, meta)


M._measurement = _measurement_v2


def build(domain_path: Path = M.DEFAULT_DOMAIN, *, source_node_index: int = 0,
          attitude_halfwidth_factor: float = 1.0 / 32.0,
          epsilon: float = 0.5):
    d = M.build(domain_path, source_node_index=source_node_index,
                attitude_halfwidth_factor=attitude_halfwidth_factor,
                epsilon=epsilon)
    d["schema_v2"] = SCHEMA
    d["dependency_preserving_S_selector_AD_used"] = True
    d["physical_S_true_forcing_must_enter_affine_b"] = True
    for row in d.get("modes", {}).values():
        sops = [op for op in row.get("operations", []) if op.get("kind") == "S"]
        row["all_S_operations_used_dependency_preserving_identity"] = bool(sops) and all(
            op.get("dependency_preserving_S_selector_identity_used") is True for op in sops
        )
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=M.DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--attitude-halfwidth-factor", type=float, default=1.0 / 32.0)
    ap.add_argument("--epsilon", type=float, default=0.5)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_node_index=a.source_node_index,
              attitude_halfwidth_factor=a.attitude_halfwidth_factor,
              epsilon=a.epsilon)
    failures = M.validate(d)
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode)
        if row is not None and row.get("all_S_operations_used_dependency_preserving_identity") is not True:
            failures.append(f"mode {mode} did not preserve every S selector cancellation")
    failures = list(dict.fromkeys(failures))
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "beta_measurement": d.get("modes", {}).get(m, {}).get("beta_measurement_design_upper"),
            "operations": d.get("modes", {}).get(m, {}).get("operation_count"),
            "S_identity": d.get("modes", {}).get(m, {}).get("all_S_operations_used_dependency_preserving_identity"),
        } for m in ("H", "A")},
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
