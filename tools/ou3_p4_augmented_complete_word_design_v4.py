#!/usr/bin/env python3
"""Pre-composition mean-value reconditioned augmented OU-III P4 word design.

V3 mean-value reconditions state and covariance *after* an accepted measurement,
but the deployed quaternion composition is evaluated inside the state update.
A broad natural interval for the Joseph correction can therefore hit the
quaternion/Cayley guard before the rigorous mean-value enclosure is applied.

This wrapper moves the same rigorous intersection onto the correction itself:

    dx(X) in natural(dx,X) intersect [dx(xc) + J_dx(X)(X-xc)]

before composing its attitude component with the physical error quaternion.
The AD derivative enclosure is unchanged.  No state/source assumption is
narrowed, no gain interval matrix is materialized, and the covariance path
remains the augmented Joseph/reset path from V3.  The raw covariance reset is
kept as a (possibly wider) rigorous enclosure; only the state correction value
box is tightened before the exact deployed update.

This remains a one-cell / one-source-node design gate and cannot promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_augmented_complete_word_design_v3 as V3
import ou3_p4_augmented_complete_word_design_v2 as V2

M = V3.M
REC = V3.REC
AUG = V2.AUG
D = V2.D
SCHEMA = 4


def _measurement_v4(mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
                    residual, center_residual, eta, entry_box, entry_center,
                    signed_word, beta_total, epsilon):
    # Keep the same augmented joint Joseph solve and correction-dependent
    # covariance reset as V2/V3.
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

    # Critical V4 step.  The center path and the full interval Jacobian already
    # exist before state composition, so intersect two independently rigorous
    # enclosures of the *same* Joseph correction.  REC leaves every derivative
    # interval bit-for-bit unchanged.
    dx_raw = cell["dx"]
    dx_center = ccell["dx"]
    dx_tight = REC.recondition_vector(dx_raw, dx_center, entry_box, entry_center)
    raw_theta_norm = M.AD._norm_upper([q.val for q in dx_raw[:3]])
    tight_theta_norm = M.AD._norm_upper([q.val for q in dx_tight[:3]])

    zout = M._state_update(z, dx_tight)
    zcout = M._state_update(zcenter, dx_center)

    selector_used = False
    if kind == "S":
        rplus = D._ad_matvec(Rm, cell["innovation_solution"])
        rcplus = D._ad_matvec(Rm, ccell["innovation_solution"])
        for i in range(3):
            zout[12 + i] = rplus[i]
            zcout[12 + i] = rcplus[i]
        selector_used = True

    # Preserve V3's post-operation reconditioning as a second dependency
    # control layer.  Covariance values are intersected only after the already
    # rigorous augmented Joseph/reset map has been formed.
    ztight = REC.recondition_vector(zout, zcout, entry_box, entry_center)
    Ptight = REC.recondition_matrix(
        cell["P_accepted_post_reset"], ccell["P_accepted_post_reset"],
        entry_box, entry_center,
    )

    meta = {
        "kind": kind,
        "inverse_backend": cell["inverse_meta"].get("inverse_backend"),
        "correction_theta_norm_upper_raw": raw_theta_norm,
        "correction_theta_norm_upper": tight_theta_norm,
        "correction_mean_value_reconditioning_used_before_quaternion_composition": True,
        "correction_mean_value_derivatives_left_unchanged": all(
            a.der == b.der for a, b in zip(dx_raw, dx_tight)
        ),
        "correction_value_width_ratio_after_reconditioning_max": V3._max_ratio(dx_raw, dx_tight),
        "covariance_path_derivative_included": cell["covariance_path_derivative_included"],
        "reset_covariance_depends_on_correction_AD": cell["reset_covariance_depends_on_correction_AD"],
        "dependency_preserving_S_selector_identity_used": selector_used,
        "physical_S_true_forcing_must_enter_affine_b": bool(selector_used),
        "K_interval_matrix_materialized": False,
        "mean_value_reconditioning_used": True,
        "state_value_width_ratio_after_reconditioning_max": V3._max_ratio(zout, ztight),
        "covariance_value_width_ratio_after_reconditioning_max": V3._max_ratio(
            V3._flatten(cell["P_accepted_post_reset"]),
            V3._flatten(Ptight),
        ),
        "mean_value_derivatives_left_unchanged": True,
        **bmeta,
    }
    return (Ptight, ztight,
            ccell["P_accepted_post_reset"], zcout,
            signed_word, beta_total, meta)


# V3._mode resolves this module global at call time; replacing it here leaves
# all scheduling, source, prediction, directional-form and validation logic
# unchanged while inserting the pre-composition correction enclosure.
V3._measurement_v3 = _measurement_v4


def build(domain_path: Path = M.DEFAULT_DOMAIN, *, source_node_index: int = 0,
          attitude_halfwidth_factor: float = 1.0 / 32.0,
          epsilon: float = 0.5):
    d = V3.build(
        domain_path,
        source_node_index=source_node_index,
        attitude_halfwidth_factor=attitude_halfwidth_factor,
        epsilon=epsilon,
    )
    d["schema_v4"] = SCHEMA
    d["qualification_v4"] = "OU3_P4_PRECOMPOSITION_RECONDITIONED_AUGMENTED_WORD_DESIGN"
    d["correction_reconditioned_before_quaternion_composition"] = True
    d["correction_derivative_enclosure_unchanged"] = True
    d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"] = False
    d["P4_USABLE_CERTIFICATE_PROMOTED"] = False
    return d


def validate(d):
    f = V3.validate(d)
    if d.get("schema_v4") != SCHEMA:
        f.append("schema_v4 mismatch")
    for key in (
        "correction_reconditioned_before_quaternion_composition",
        "correction_derivative_enclosure_unchanged",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode)
        if row is None:
            continue
        for op in row.get("operations", []):
            if op.get("correction_mean_value_reconditioning_used_before_quaternion_composition") is not True:
                f.append(f"mode {mode} operation lacks pre-composition correction reconditioning")
            if op.get("correction_mean_value_derivatives_left_unchanged") is not True:
                f.append(f"mode {mode} correction derivative enclosure changed")
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=M.DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--attitude-halfwidth-factor", type=float, default=1.0 / 32.0)
    ap.add_argument("--epsilon", type=float, default=0.5)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(
        a.domain,
        source_node_index=a.source_node_index,
        attitude_halfwidth_factor=a.attitude_halfwidth_factor,
        epsilon=a.epsilon,
    )
    vf = validate(d)
    d["validation_pass_v4"] = not vf
    d["validation_failures_v4"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "beta_measurement": d.get("modes", {}).get(m, {}).get("beta_measurement_design_upper"),
            "operations": d.get("modes", {}).get(m, {}).get("operation_count"),
            "max_raw_correction": max([
                float(op.get("correction_theta_norm_upper_raw", 0.0))
                for op in d.get("modes", {}).get(m, {}).get("operations", [])
            ] or [0.0]),
            "max_tight_correction": max([
                float(op.get("correction_theta_norm_upper", 0.0))
                for op in d.get("modes", {}).get(m, {}).get("operations", [])
            ] or [0.0]),
        } for m in ("H", "A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
