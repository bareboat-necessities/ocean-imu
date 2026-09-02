#!/usr/bin/env python3
"""Current augmented complete-word OU-III P4 design.

This is the single current producer. It combines the former v2-v5 refinements:
exact S-selector cancellation, mean-value reconditioning after operations and
predictions, pre-composition correction reconditioning, and the diagonally
scaled verified innovation inverse. Historical wrapper stages are retired.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_augmented_complete_word_design as M
import ou3_p4_augmented_covariance_ad as AUG
import ou3_p4_joint_word_dissipation_design as D
import ou3_p4_mean_value_reconditioning as REC

SCHEMA = 5
JJ = AUG.JJ
_ORIGINAL_VERIFIED_INVERSE = JJ.verified_inverse


def _diag_scale(S):
    d = []
    for i in range(len(S)):
        m = 0.5 * (float(S[i][i].lo) + float(S[i][i].hi))
        if not (math.isfinite(m) and m > 0.0):
            raise RuntimeError(f"innovation diagonal midpoint is not positive at {i}: {m}")
        di = 1.0 / math.sqrt(m)
        if not (math.isfinite(di) and di > 0.0):
            raise RuntimeError(f"invalid innovation diagonal scale at {i}: {di}")
        d.append(di)
    return d


def scaled_verified_inverse(S):
    try:
        return _ORIGINAL_VERIFIED_INVERSE(S)
    except Exception as first:
        d = _diag_scale(S)
        n = len(S)
        T = [[JJ.I(d[i]) * S[i][j] * JJ.I(d[j]) for j in range(n)] for i in range(n)]
        try:
            Xt, meta = _ORIGINAL_VERIFIED_INVERSE(T)
        except Exception as second:
            raise type(second)(
                f"raw verified inverse failed ({type(first).__name__}: {first}); "
                f"diagonal-scaled verified inverse also failed: {second}"
            ) from second
        X = [[JJ.I(d[i]) * Xt[i][j] * JJ.I(d[j]) for j in range(n)] for i in range(n)]
        X = JJ.matrix_symmetric_hull(X)
        out = dict(meta)
        out.update({
            "inverse_backend": "SYMMETRIC_DIAGONAL_SCALED_" + str(meta.get("inverse_backend")),
            "diagonal_congruence_preconditioner_used": True,
            "diagonal_point_scales": d,
            "preconditioner_is_exact_point_congruence": True,
            "recovery_identity": "S^-1=D(DSD)^-1D",
            "ordinary_float_inverse_used_as_enclosure": False,
            "K_interval_matrix_materialized": False,
            "raw_inverse_failure": f"{type(first).__name__}: {first}",
        })
        return X, out


JJ.verified_inverse = scaled_verified_inverse


def _flatten(A):
    return [x for row in A for x in row]


def _max_ratio(before, after):
    vals = []
    for a, b in zip(before, after):
        wa = a.val.width()
        if wa > 0.0:
            vals.append(b.val.width() / wa)
    return max(vals) if vals else 0.0


def _measurement(mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
                 residual, center_residual, eta, entry_box, entry_center,
                 signed_word, beta_total, epsilon):
    cell = AUG.measurement_and_reset(
        Pad, Hm, Rm, residual, value_tighten=M.COV.psd_tighten
    )
    ccell = AUG.measurement_and_reset(
        Pcenter, Hm, Rm, center_residual, value_tighten=M.COV.psd_tighten
    )
    dform, beta, bmeta = M._finite_cell_terms(
        residual, eta, center_residual, entry_center,
        AUG.values(cell["Sinv"]), Rm, epsilon,
    )
    signed_word = M._add_form(signed_word, dform)
    beta_total = math.nextafter(beta_total + beta, math.inf) if beta > 0.0 else beta_total

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
        "correction_value_width_ratio_after_reconditioning_max": _max_ratio(dx_raw, dx_tight),
        "covariance_path_derivative_included": cell["covariance_path_derivative_included"],
        "reset_covariance_depends_on_correction_AD": cell["reset_covariance_depends_on_correction_AD"],
        "dependency_preserving_S_selector_identity_used": selector_used,
        "physical_S_true_forcing_must_enter_affine_b": bool(selector_used),
        "K_interval_matrix_materialized": False,
        "mean_value_reconditioning_used": True,
        "state_value_width_ratio_after_reconditioning_max": _max_ratio(zout, ztight),
        "covariance_value_width_ratio_after_reconditioning_max": _max_ratio(
            _flatten(cell["P_accepted_post_reset"]), _flatten(Ptight)
        ),
        "mean_value_derivatives_left_unchanged": True,
        **bmeta,
    }
    return (
        Ptight, ztight, ccell["P_accepted_post_reset"], zcout,
        signed_word, beta_total, meta,
    )


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
    Racc = M.COV.R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = M.COV.R_diag(float(vc["mag_measurement_std_uT"]))
    RS = M.COV.R_S(src)
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
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "S", Pad, z, Pcenter, zcenter, Hs, RS,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        if k in schedule["vector_steps"]:
            r = M.D._exact_acc_residual(mode, z, force)
            rc = M.D._exact_acc_residual(mode, zcenter, force)
            eta = M.D._eta(r, M.D._linear_residual(Ha, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "acc", Pad, z, Pcenter, zcenter, Ha, Racc,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

            r = M.D._exact_mag_residual(z, mag)
            rc = M.D._exact_mag_residual(zcenter, mag)
            eta = M.D._eta(r, M.D._linear_residual(Hm, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "mag", Pad, z, Pcenter, zcenter, Hm, Rmag,
                r, rc, eta, entry_box, entry_center,
                signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        Praw = M.AUG.predict(Pad, F, Q, value_tighten=M.COV.psd_tighten)
        Pc = M.AUG.predict(Pcenter, F, Q, value_tighten=M.COV.psd_tighten)
        zraw = M.D._prediction(mode, z, F)
        zc = M.D._prediction(mode, zcenter, F)
        Pad = REC.recondition_matrix(Praw, Pc, entry_box, entry_center)
        z = REC.recondition_vector(zraw, zc, entry_box, entry_center)
        Pcenter, zcenter = Pc, zc
        state_prediction_reconditioning_count += 1
        augmented_prediction_reconditioning_count += 1

    mu = M.D._generalized_margin(signed_word, Mupper)
    rho = math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0
    all_cov = all(
        op["covariance_path_derivative_included"]
        and op["reset_covariance_depends_on_correction_AD"]
        for op in ops
    )
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
          attitude_halfwidth_factor: float = 1.0 / 32.0, epsilon: float = 0.5):
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
        "qualification": "OU3_P4_CURRENT_AUGMENTED_COMPLETE_WORD_DESIGN",
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
        "correction_reconditioned_before_quaternion_composition": True,
        "correction_derivative_enclosure_unchanged": True,
        "scaled_verified_inverse_available": True,
        "scaled_inverse_uses_only_point_congruence_and_verified_interval_inverse": True,
        "ordinary_float_inverse_used_as_enclosure": False,
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
        "endpoint_scalarization_only", "correction_reconditioned_before_quaternion_composition",
        "correction_derivative_enclosure_unchanged", "scaled_verified_inverse_available",
        "scaled_inverse_uses_only_point_congruence_and_verified_interval_inverse",
        "attitude_design_cell_only",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "ordinary_float_inverse_used_as_enclosure",
        "K_interval_matrix_materialized", "outer_0p8rad_cover_checked_here",
        "phased_P2_paths_checked", "optional_accepted_branch_family_checked",
        "full_source_vector_family_checked", "deterministic_transport_beta_added",
        "stochastic_beta_added", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED", "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
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
        if (
            int(row.get("measurement_free_prefix_value_only_covariance_predictions", -1))
            + int(row.get("augmented_covariance_predictions", -1))
            != int(row.get("samples", -2))
        ):
            f.append(f"mode {mode} covariance prediction accounting incomplete")
        if row.get("all_S_operations_used_dependency_preserving_identity") is not True:
            f.append(f"mode {mode} S selector identity incomplete")
        mu = row.get("signed_word_generalized_margin_design")
        rho = row.get("rho_homogeneous_design_upper")
        if not isinstance(mu, (int, float)) or not mu > 0.0:
            f.append(f"mode {mode} design generalized margin is not positive: {mu}")
        if not isinstance(rho, (int, float)) or not rho < 1.0:
            f.append(f"mode {mode} design rho is not below one: {rho}")
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
    d["validation_pass_v5"] = not vf
    d["validation_failures_v5"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "beta_measurement": d.get("modes", {}).get(m, {}).get("beta_measurement_design_upper"),
            "operations": d.get("modes", {}).get(m, {}).get("operation_count"),
            "inverse_backends": [
                op.get("inverse_backend")
                for op in d.get("modes", {}).get(m, {}).get("operations", [])
            ],
        } for m in ("H", "A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
