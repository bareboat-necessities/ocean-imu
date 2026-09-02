#!/usr/bin/env python3
"""Mean-value-reconditioned augmented complete-word OU-III P4 design.

V2 preserves the exact S-selector cancellation inside each accepted S=0 update,
but natural interval evaluation can lose the newly-created correlations again
on the following prediction.  This backend keeps the same joint Joseph,
augmented covariance AD, exact S identity and finite-cell affine ledger, and
adds rigorous mean-value reconditioning after every measurement/reset and every
prediction:

    y(X) in natural(y,X) intersection [ y(xc) + J_y(X)(X-xc) ].

The derivative enclosure is never replaced or tightened; only the value range
is intersected with a second rigorous enclosure of the same map.  The center
path is propagated through the same source cell and operation sequence.

Before the first accepted measurement the covariance is independent of the
word-entry error, so all of its AD derivatives are identically zero.  This
version therefore propagates that measurement-free prefix with the existing
validated value-only covariance map and lifts the covariance to augmented AD
only at the first measurement.  This is algebraically identical to carrying
hundreds of zero derivative vectors and makes full-cell/source sweeps practical.

This file remains a design gate for one small attitude cell / one P2 source node.
It may establish a design rho<1, but it intentionally cannot promote P4 until
the outer 0.8-rad cover, complete source/vector family, phased P2 paths and all
affine forcing terms are closed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_augmented_complete_word_design as M
import ou3_p4_augmented_complete_word_design_v2 as V2
import ou3_p4_mean_value_reconditioning as REC

SCHEMA = 3


def _flatten(A):
    return [x for row in A for x in row]


def _max_ratio(before, after):
    vals = []
    for a, b in zip(before, after):
        wa = a.val.width()
        if wa > 0.0:
            vals.append(b.val.width() / wa)
    return max(vals) if vals else 0.0


def _measurement_v3(mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
                    residual, center_residual, eta, entry_box, entry_center,
                    signed_word, beta_total, epsilon):
    raw = V2._measurement_v2(
        mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
        residual, center_residual, eta, entry_center,
        signed_word, beta_total, epsilon,
    )
    Praw, zraw, Pc, zc, signed_word, beta_total, meta = raw
    ztight = REC.recondition_vector(zraw, zc, entry_box, entry_center)
    Ptight = REC.recondition_matrix(Praw, Pc, entry_box, entry_center)
    meta = dict(meta)
    meta.update({
        "mean_value_reconditioning_used": True,
        "state_value_width_ratio_after_reconditioning_max": _max_ratio(zraw, ztight),
        "covariance_value_width_ratio_after_reconditioning_max": _max_ratio(
            _flatten(Praw), _flatten(Ptight)
        ),
        "mean_value_derivatives_left_unchanged": True,
    })
    return Ptight, ztight, Pc, zc, signed_word, beta_total, meta


def _mode(mode: str, path: Path, domain: dict, source_node_index: int,
          attitude_halfwidth_factor: float, epsilon: float):
    M.G.CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    src = M.NODES.h18_source_cell(source_node_index, M.NODES.build())
    F, Q, _meta = M.G2.corrected_zero_rate_transition_process(mode, src, domain)

    Ppre = M.G._gauged_golive_covariance(mode, src, path)
    Pentry = M._predict_value_cov(Ppre, F, Q)

    qouter = 2.0 * math.tan(0.80 / 2.0)
    hw = qouter * float(attitude_halfwidth_factor)
    cbox = [(-hw, hw), (-hw, hw), (-hw, hw)]
    zpre = M.D._initial_ad(mode, domain, cbox)
    z = M.POST._rebase_postprediction(M.D._prediction(mode, zpre, F))
    entry_box = [q.val for q in z]
    entry_center = [M._mid(q.val) for q in z]
    zcenter = M._centerize(z)

    Mupper, metric_meta = M.POST._entry_information_upper(mode)

    force, mag = M.D._canonical_vectors(domain)
    Ha, Hm, Hs = M.D._H_acc(mode, force, n), M.D._H_mag(mag, n), M.D._H_S(n)
    vc = M.VECTOR.build()["configured_measurement_bounds"]
    Racc = M.H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = M.H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = M.H._R_S(src)
    h = float(src["dt_s"])
    words = M.WORDS.build(path)
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = M.D._schedule(path, samples, h)
    event_steps = sorted(set(schedule["S_steps"] + schedule["vector_steps"]))
    if not event_steps:
        raise RuntimeError("complete word contains no accepted measurements")
    first_event = int(event_steps[0])

    signed_word = M.D._zero(n, n)
    beta_total = 0.0
    ops = []
    state_prediction_reconditioning_count = 0

    # Covariance derivatives with respect to the continuous entry error are
    # exactly zero until the first correction-dependent reset.  Carry only the
    # validated interval covariance values over this prefix.  The state still
    # carries its full entry Jacobian and is reconditioned on every step.
    Pvalue = Pentry
    for _k in range(first_event):
        Pvalue = M._predict_value_cov(Pvalue, F, Q)
        zraw = M.D._prediction(mode, z, F)
        zc = M.D._prediction(mode, zcenter, F)
        z = REC.recondition_vector(zraw, zc, entry_box, entry_center)
        zcenter = zc
        state_prediction_reconditioning_count += 1

    Pad = M.AUG.constant_matrix(Pvalue, n)
    Pcenter = M.AUG.constant_matrix(Pvalue, n)
    augmented_prediction_reconditioning_count = 0

    for k in range(first_event, samples):
        if k in schedule["S_steps"]:
            r = [z[12 + i] for i in range(3)]
            rc = [zcenter[12 + i] for i in range(3)]
            eta = [M.AD.constant(0.0, n) for _ in range(3)]
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement_v3(
                mode, "S", Pad, z, Pcenter, zcenter, Hs, RS,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        if k in schedule["vector_steps"]:
            r = M.D._exact_acc_residual(mode, z, force)
            rc = M.D._exact_acc_residual(mode, zcenter, force)
            eta = M.D._eta(r, M.D._linear_residual(Ha, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement_v3(
                mode, "acc", Pad, z, Pcenter, zcenter, Ha, Racc,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

            r = M.D._exact_mag_residual(z, mag)
            rc = M.D._exact_mag_residual(zcenter, mag)
            eta = M.D._eta(r, M.D._linear_residual(Hm, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement_v3(
                mode, "mag", Pad, z, Pcenter, zcenter, Hm, Rmag,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        Praw = M.AUG.predict(Pad, F, Q, value_tighten=M.H._psd_tighten)
        Pc = M.AUG.predict(Pcenter, F, Q, value_tighten=M.H._psd_tighten)
        zraw = M.D._prediction(mode, z, F)
        zc = M.D._prediction(mode, zcenter, F)
        Pad = REC.recondition_matrix(Praw, Pc, entry_box, entry_center)
        z = REC.recondition_vector(zraw, zc, entry_box, entry_center)
        Pcenter, zcenter = Pc, zc
        state_prediction_reconditioning_count += 1
        augmented_prediction_reconditioning_count += 1

    mu = M.D._generalized_margin(signed_word, Mupper)
    rho = math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0
    all_cov = all(op["covariance_path_derivative_included"] and
                  op["reset_covariance_depends_on_correction_AD"] for op in ops)
    no_k = all(op["K_interval_matrix_materialized"] is False for op in ops)
    all_mv = all(op.get("mean_value_reconditioning_used") is True for op in ops)
    all_s = [op for op in ops if op.get("kind") == "S"]
    return {
        "dimension": n,
        "source_node_index": int(source_node_index),
        "word_boundary": "POST_PREDICTION_PRE_MEASUREMENT",
        "attitude_design_cell_halfwidth_cayley": hw,
        "attitude_design_cell_halfwidth_factor_of_outer_q": float(attitude_halfwidth_factor),
        "outer_angle_rad_reference": 0.80,
        "samples": samples,
        "schedule": schedule,
        "first_measurement_step": first_event,
        "measurement_free_prefix_covariance_derivative_identically_zero": True,
        "measurement_free_prefix_value_only_covariance_predictions": first_event,
        "augmented_covariance_predictions": augmented_prediction_reconditioning_count,
        "operation_count": len(ops),
        "operations": ops,
        "entry_information_metric_upper": metric_meta,
        "complete_word_jacobian_covariance_paths_certified_in_design_branch": all_cov,
        "directional_forms_accumulated_before_scalarization": True,
        "finite_cell_affine_penalties_accumulated": True,
        "mean_value_reconditioning_used_after_every_operation": all_mv,
        "mean_value_reconditioning_used_after_every_state_prediction": state_prediction_reconditioning_count == samples,
        "state_prediction_reconditioning_count": state_prediction_reconditioning_count,
        "all_S_operations_used_dependency_preserving_identity": bool(all_s) and all(
            op.get("dependency_preserving_S_selector_identity_used") is True for op in all_s
        ),
        "signed_word_generalized_margin_design": mu,
        "rho_homogeneous_design_upper": rho,
        "beta_measurement_design_upper": beta_total,
        "deterministic_transport_beta_added": False,
        "stochastic_beta_added": False,
        "K_interval_matrix_materialized": not no_k,
        "design_rho_beta_pair_complete_for_P5": False,
    }


def build(domain_path: Path = M.DEFAULT_DOMAIN, *, source_node_index: int = 0,
          attitude_halfwidth_factor: float = 1.0 / 32.0,
          epsilon: float = 0.5):
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    words = M.WORDS.build(path)
    failures = [f"word: {x}" for x in M.WORDS.validate(words)]
    modes = {}
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode(
                mode, path, domain, int(source_node_index),
                float(attitude_halfwidth_factor), float(epsilon),
            )
        except Exception as exc:
            failures.append(f"{mode}: {type(exc).__name__}: {exc}")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_MEAN_VALUE_RECONDITIONED_AUGMENTED_COMPLETE_WORD_DESIGN",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "design_only_not_theorem_promotion": True,
        "joint_P_H_R_r_used": True,
        "augmented_covariance_AD_used": True,
        "state_dependent_covariance_gain_reset_paths_differentiated": True,
        "dependency_preserving_S_selector_AD_used": True,
        "mean_value_reconditioning_used": True,
        "mean_value_reconditioning_is_intersection_of_two_rigorous_enclosures": True,
        "measurement_free_covariance_AD_fast_forward_is_exact": True,
        "finite_cell_affine_Joseph_used": True,
        "directional_forms_accumulated_before_scalarization": True,
        "endpoint_scalarization_only": True,
        "K_interval_matrix_materialized": False,
        "attitude_design_cell_only": True,
        "outer_0p8rad_cover_checked_here": False,
        "phased_P2_paths_checked": False,
        "optional_accepted_branch_family_checked": False,
        "full_source_vector_family_checked": False,
        "deterministic_transport_beta_added": False,
        "stochastic_beta_added": False,
        "modes": modes,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "design_only_not_theorem_promotion",
        "joint_P_H_R_r_used", "augmented_covariance_AD_used",
        "state_dependent_covariance_gain_reset_paths_differentiated",
        "dependency_preserving_S_selector_AD_used", "mean_value_reconditioning_used",
        "mean_value_reconditioning_is_intersection_of_two_rigorous_enclosures",
        "measurement_free_covariance_AD_fast_forward_is_exact",
        "finite_cell_affine_Joseph_used", "directional_forms_accumulated_before_scalarization",
        "endpoint_scalarization_only", "attitude_design_cell_only",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "K_interval_matrix_materialized",
        "outer_0p8rad_cover_checked_here", "phased_P2_paths_checked",
        "optional_accepted_branch_family_checked", "full_source_vector_family_checked",
        "deterministic_transport_beta_added", "stochastic_beta_added",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE", "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode)
        if row is None:
            f.append(f"mode {mode} missing")
            continue
        if row.get("complete_word_jacobian_covariance_paths_certified_in_design_branch") is not True:
            f.append(f"mode {mode} covariance AD path incomplete")
        if row.get("mean_value_reconditioning_used_after_every_operation") is not True:
            f.append(f"mode {mode} operation reconditioning incomplete")
        if row.get("mean_value_reconditioning_used_after_every_state_prediction") is not True:
            f.append(f"mode {mode} state prediction reconditioning incomplete")
        if int(row.get("measurement_free_prefix_value_only_covariance_predictions", -1)) + int(row.get("augmented_covariance_predictions", -1)) != int(row.get("samples", -2)):
            f.append(f"mode {mode} covariance prediction accounting incomplete")
        if row.get("all_S_operations_used_dependency_preserving_identity") is not True:
            f.append(f"mode {mode} S selector identity incomplete")
        mu = row.get("signed_word_generalized_margin_design")
        rho = row.get("rho_homogeneous_design_upper")
        if not isinstance(mu, (int, float)) or not mu > 0.0:
            f.append(f"mode {mode} design generalized margin is not positive: {mu}")
        if not isinstance(rho, (int, float)) or not rho < 1.0:
            f.append(f"mode {mode} design rho is not below one: {rho}")
    return list(dict.fromkeys(f))


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
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "beta_measurement": d.get("modes", {}).get(m, {}).get("beta_measurement_design_upper"),
            "operations": d.get("modes", {}).get(m, {}).get("operation_count"),
            "first_measurement": d.get("modes", {}).get(m, {}).get("first_measurement_step"),
            "value_only_prefix": d.get("modes", {}).get(m, {}).get("measurement_free_prefix_value_only_covariance_predictions"),
            "augmented_predictions": d.get("modes", {}).get(m, {}).get("augmented_covariance_predictions"),
            "S_identity": d.get("modes", {}).get(m, {}).get("all_S_operations_used_dependency_preserving_identity"),
        } for m in ("H", "A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
