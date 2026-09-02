#!/usr/bin/env python3
"""Augmented finite-cell complete-word Joseph design for OU-III P4.

This is an end-to-end design probe for the proof calculus requested by P4.  It
combines four previously separate rigorous primitives:

* the post-prediction/process-UCC word boundary;
* the co-gauged source realization that keeps vector/covariance orientation
  correlated and makes the verified innovation solve meaningful;
* augmented interval AD of covariance, innovation inverse, Joseph posterior and
  attitude-reset covariance; and
* finite-cell affine Joseph dissipation, so an off-centre subdivision cell is
  not silently treated as a homogeneous neighbourhood of zero.

For a fixed P2 source node, zero-rate co-gauged source realization and canonical
PE vectors, each mandatory measurement contributes

    r^T S^-1 r - eta^T R^-1 eta
      >= x^T D_i x - beta_i,

where the residual Jacobian is with respect to the *word-entry* coordinates and
includes every earlier state-dependent covariance/gain/reset path.  The
positive directional matrices D_i are accumulated across the full H=18 or A=21
word and scalarized once at the endpoint against the process-UCC entry metric.
The affine penalties beta_i are accumulated separately.  Hence the diagnostic
emits the first end-to-end design quantities of the intended form

    W_next <= rho W_entry + beta_measurement,

with rho=1-mu when the accumulated generalized margin mu is positive.

Important limits: this file checks only one deliberately small attitude design
cell, one P2 source node and one canonical vector realization; it does not yet
enumerate the 0.8-rad outer cover, the phased P2 automaton, optional accepted /
rejected branches, or the full source vector family.  Deterministic transport
forcing and stochastic terms are also not yet added to beta.  It therefore
cannot promote P4 or compute P5 N_inner.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_interval_ad as AD
import ou3_p4_augmented_covariance_ad as AUG
import ou3_p4_cell_affine_joseph as AFF
import ou3_p4_joint_word_dissipation_design as D
import ou3_p4_joint_word_postprediction_design as POST
import ou3_p4_joint_word_gauge_design as G
import ou3_p4_joint_word_gauge_design_v2 as G2
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_vector_uco_certificate as VECTOR
import ou3_implementation_word_language as WORDS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _mid(x: Interval) -> float:
    return 0.5 * (float(x.lo) + float(x.hi))


def _centerize(z):
    n = len(z)
    return [AD.constant(_mid(q.val), n) for q in z]


def _state_update(z, dx):
    out = list(z)
    out[:3] = AD.deployed_correct_cayley_right(z[:3], [-q for q in dx[:3]])
    for i in range(3, len(z)):
        out[i] = z[i] - dx[i]
    return out


def _predict_value_cov(P, F, Q):
    return H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P), matrix_transpose(F)), Q))


def _add_form(A, B):
    return D._form_add(A, B)


def _finite_cell_terms(residual, eta, center_residual, entry_center,
                       Sinv_value, R, epsilon: float):
    """Directional matrix + affine beta for one actual word residual map."""
    J = AD.jacobian(residual)
    retained = math.nextafter(1.0 - float(epsilon), -math.inf)
    directional = AFF._scale_form(D._form(J, Sinv_value), retained)
    ebox = AFF._affine_intercept(AD.values(center_residual), J, entry_center)
    einfo = AFF._quad_upper(ebox, Sinv_value)
    young = math.nextafter(retained / float(epsilon), math.inf)
    beta_affine = 0.0 if einfo == 0.0 else math.nextafter(young * einfo, math.inf)
    beta_eta = AFF._eta_Rinv_upper(AD.values(eta), R)
    beta = 0.0 if beta_affine == 0.0 and beta_eta == 0.0 else math.nextafter(beta_affine + beta_eta, math.inf)
    return directional, beta, {
        "affine_intercept_information_upper": einfo,
        "beta_affine_intercept_upper": beta_affine,
        "beta_nonlinear_eta_upper": beta_eta,
        "beta_total_upper": beta,
    }


def _measurement(mode, kind, Pad, z, Pcenter, zcenter, Hm, Rm,
                 residual, center_residual, eta, entry_center,
                 signed_word, beta_total, epsilon):
    # The augmented cell differentiates P,H,S^-1, Joseph posterior and reset
    # covariance with respect to the original word-entry error.  H/R are fixed
    # within this source branch; source enumeration is a separate discrete axis.
    cell = AUG.measurement_and_reset(
        Pad, Hm, Rm, residual, value_tighten=H._psd_tighten
    )
    ccell = AUG.measurement_and_reset(
        Pcenter, Hm, Rm, center_residual, value_tighten=H._psd_tighten
    )
    dform, beta, bmeta = _finite_cell_terms(
        residual, eta, center_residual, entry_center,
        AUG.values(cell["Sinv"]), Rm, epsilon,
    )
    signed_word = _add_form(signed_word, dform)
    beta_total = math.nextafter(beta_total + beta, math.inf) if beta > 0.0 else beta_total
    zout = _state_update(z, cell["dx"])
    zcout = _state_update(zcenter, ccell["dx"])
    meta = {
        "kind": kind,
        "inverse_backend": cell["inverse_meta"].get("inverse_backend"),
        "correction_theta_norm_upper": AD._norm_upper([q.val for q in cell["dx"][:3]]),
        "covariance_path_derivative_included": cell["covariance_path_derivative_included"],
        "reset_covariance_depends_on_correction_AD": cell["reset_covariance_depends_on_correction_AD"],
        "K_interval_matrix_materialized": False,
        **bmeta,
    }
    return (cell["P_accepted_post_reset"], zout,
            ccell["P_accepted_post_reset"], zcout,
            signed_word, beta_total, meta)


def _mode(mode: str, path: Path, domain: dict, source_node_index: int,
          attitude_halfwidth_factor: float, epsilon: float):
    G.CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    src = NODES.h18_source_cell(source_node_index, NODES.build())
    F, Q, _meta = G2.corrected_zero_rate_transition_process(mode, src, domain)

    # Exact source covariance at the post-prediction/pre-measurement boundary.
    Ppre = G._gauged_golive_covariance(mode, src, path)
    Pentry = _predict_value_cov(Ppre, F, Q)

    # Deliberately small design cell.  This validates the full calculus before
    # the exact 96-cell outer cover is connected.  It is NOT a basin result.
    qouter = 2.0 * math.tan(0.80 / 2.0)
    hw = qouter * float(attitude_halfwidth_factor)
    cbox = [(-hw, hw), (-hw, hw), (-hw, hw)]
    zpre = D._initial_ad(mode, domain, cbox)
    z = POST._rebase_postprediction(D._prediction(mode, zpre, F))
    zcenter = _centerize(z)
    entry_center = [_mid(q.val) for q in z]

    Pad = AUG.constant_matrix(Pentry, n)
    Pcenter = AUG.constant_matrix(Pentry, n)
    Mupper, metric_meta = POST._entry_information_upper(mode)

    force, mag = D._canonical_vectors(domain)
    Ha, Hm, Hs = D._H_acc(mode, force, n), D._H_mag(mag, n), D._H_S(n)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    h = float(src["dt_s"])
    words = WORDS.build(path)
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = D._schedule(path, samples, h)

    signed_word = D._zero(n, n)
    beta_total = 0.0
    ops = []

    for k in range(samples):
        if k in schedule["S_steps"]:
            r = [z[12 + i] for i in range(3)]
            rc = [zcenter[12 + i] for i in range(3)]
            eta = [AD.constant(0.0, n) for _ in range(3)]
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "S", Pad, z, Pcenter, zcenter, Hs, RS,
                r, rc, eta, entry_center, signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        if k in schedule["vector_steps"]:
            r = D._exact_acc_residual(mode, z, force)
            rc = D._exact_acc_residual(mode, zcenter, force)
            eta = D._eta(r, D._linear_residual(Ha, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "acc", Pad, z, Pcenter, zcenter, Ha, Racc,
                r, rc, eta, entry_center, signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

            r = D._exact_mag_residual(z, mag)
            rc = D._exact_mag_residual(zcenter, mag)
            eta = D._eta(r, D._linear_residual(Hm, z))
            Pad, z, Pcenter, zcenter, signed_word, beta_total, meta = _measurement(
                mode, "mag", Pad, z, Pcenter, zcenter, Hm, Rmag,
                r, rc, eta, entry_center, signed_word, beta_total, epsilon,
            )
            ops.append({"step": k, **meta})

        # Meet the next tile on the same post-prediction boundary class.
        Pad = AUG.predict(Pad, F, Q, value_tighten=H._psd_tighten)
        Pcenter = AUG.predict(Pcenter, F, Q, value_tighten=H._psd_tighten)
        z = D._prediction(mode, z, F)
        zcenter = D._prediction(mode, zcenter, F)

    mu = D._generalized_margin(signed_word, Mupper)
    rho = math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0
    all_cov = all(op["covariance_path_derivative_included"] and
                  op["reset_covariance_depends_on_correction_AD"] for op in ops)
    no_k = all(op["K_interval_matrix_materialized"] is False for op in ops)
    return {
        "dimension": n,
        "source_node_index": int(source_node_index),
        "word_boundary": "POST_PREDICTION_PRE_MEASUREMENT",
        "attitude_design_cell_halfwidth_cayley": hw,
        "attitude_design_cell_halfwidth_factor_of_outer_q": float(attitude_halfwidth_factor),
        "outer_angle_rad_reference": 0.80,
        "samples": samples,
        "schedule": schedule,
        "operation_count": len(ops),
        "operations": ops,
        "entry_information_metric_upper": metric_meta,
        "complete_word_jacobian_covariance_paths_certified_in_design_branch": all_cov,
        "directional_forms_accumulated_before_scalarization": True,
        "finite_cell_affine_penalties_accumulated": True,
        "signed_word_generalized_margin_design": mu,
        "rho_homogeneous_design_upper": rho,
        "beta_measurement_design_upper": beta_total,
        "deterministic_transport_beta_added": False,
        "stochastic_beta_added": False,
        "K_interval_matrix_materialized": not no_k,
        "design_rho_beta_pair_complete_for_P5": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          attitude_halfwidth_factor: float = 1.0 / 32.0,
          epsilon: float = 0.5) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    words = WORDS.build(path)
    failures = [f"word: {x}" for x in WORDS.validate(words)]
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
        "qualification": "OU3_P4_AUGMENTED_COMPLETE_WORD_FINITE_CELL_DESIGN",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "design_only_not_theorem_promotion": True,
        "P1_changed": False,
        "P2_phased_graph_required_for_promotion": True,
        "P3_linear_foundation_retained": True,
        "P3_delta_used_as_physical_basin": False,
        "joint_P_H_R_r_used": True,
        "augmented_covariance_AD_used": True,
        "state_dependent_covariance_gain_reset_paths_differentiated": True,
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


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "design_only_not_theorem_promotion",
        "P2_phased_graph_required_for_promotion", "P3_linear_foundation_retained",
        "joint_P_H_R_r_used", "augmented_covariance_AD_used",
        "state_dependent_covariance_gain_reset_paths_differentiated",
        "finite_cell_affine_Joseph_used",
        "directional_forms_accumulated_before_scalarization", "endpoint_scalarization_only",
        "attitude_design_cell_only",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "P1_changed",
        "P3_delta_used_as_physical_basin", "K_interval_matrix_materialized",
        "outer_0p8rad_cover_checked_here", "phased_P2_paths_checked",
        "optional_accepted_branch_family_checked", "full_source_vector_family_checked",
        "deterministic_transport_beta_added", "stochastic_beta_added",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
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
            f.append(f"mode {mode} covariance-path Jacobian incomplete")
        if row.get("K_interval_matrix_materialized") is not False:
            f.append(f"mode {mode} materialized K")
        if not isinstance(row.get("beta_measurement_design_upper"), (int, float)) or not math.isfinite(float(row["beta_measurement_design_upper"])):
            f.append(f"mode {mode} finite beta missing")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--attitude-halfwidth-factor", type=float, default=1.0 / 32.0)
    ap.add_argument("--epsilon", type=float, default=0.5)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_node_index=a.source_node_index,
              attitude_halfwidth_factor=a.attitude_halfwidth_factor,
              epsilon=a.epsilon)
    failures = validate(d)
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
            "covariance_paths": d.get("modes", {}).get(m, {}).get("complete_word_jacobian_covariance_paths_certified_in_design_branch"),
        } for m in ("H", "A")},
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
