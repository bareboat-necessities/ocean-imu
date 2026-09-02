#!/usr/bin/env python3
"""Post-prediction complete-word signed-Joseph design for OU-III P4.

The first joint-word design used the broad goLive covariance *upper box* as if
it were a two-sided interval family and tried to invert it to define the word
entry metric.  That is not the right object: the box was constructed only to
upper-bound covariance and intentionally permits singular endpoint matrices.

This diagnostic moves the recurrent Live tile boundary to

    post-prediction / pre-measurement.

At every such boundary the source prediction gives

    P^- = F P F^T + Q >= Q > 0.

The full-process UCC producer certifies independent block floors

    Q_theta,bg >= q_ab I,
    Q_v,p,S,aw >= q_tr I,
    Q_ba >= q_ba I  (A mode),

so the exact entry information metric obeys

    (P^-)^-1 <= diag(q_ab^-1 I6, q_tr^-1 I12 [, q_ba^-1 I3]).

The right-hand side is therefore a rigorous *upper* comparison metric.  If the
complete signed Joseph word form D satisfies

    D >= mu M_upper,  mu>0,

then automatically D >= mu (P^-)^-1.  Scalarization occurs only once, at the
word endpoint; P3's tiny relative Riccati margin is not used.

As in ``ou3_p4_joint_word_dissipation_design.py``, this file is still a design
probe.  It checks one physical P2 node, one outer-attitude cover cell and one
canonical PE-vector realization.  Promotion remains blocked until the same
calculus is validated on the full phased P2 automaton and complete source
geometry.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_full_process_ucc as PROCESS
import ou3_implementation_word_language as WORDS
import ou3_interval_ad as AD
import ou3_p4_candidate_full_word as CAND
import ou3_p4_joint_word_dissipation_design as D
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _entry_information_upper(mode: str):
    """Block Loewner upper for the exact post-prediction information metric."""
    proc = PROCESS.build()
    vf = PROCESS.validate(proc)
    if vf:
        raise RuntimeError(f"full process UCC prerequisite failed: {vf}")
    qab = float(proc["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    qtr = float(proc["translation"]["Q_translation_lambda_min_lower"])
    qba = float(proc["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    if not (qab > 0.0 and qtr > 0.0 and qba > 0.0):
        raise RuntimeError("process block lower bound lost positivity")
    vals = [1.0 / qab] * 6 + [1.0 / qtr] * 12
    if mode == "A":
        vals += [1.0 / qba] * 3
    n = 18 if mode == "H" else 21
    if len(vals) != n:
        raise RuntimeError("entry process metric dimension mismatch")
    M = D._zero(n, n)
    for i, x in enumerate(vals):
        M[i][i] = Interval.outward_bounds(
            math.nextafter(x, -math.inf), math.nextafter(x, math.inf)
        )
    return M, {
        "attitude_gyro_bias_process_lambda_min_lower": qab,
        "translation_process_lambda_min_lower": qtr,
        "active_accel_bias_process_lambda_min_lower": qba if mode == "A" else None,
        "comparison": "Pminus_inverse <= blockdiag(1/q_ab I6,1/q_tr I12[,1/q_ba I3])",
    }


def _rebase_postprediction(z):
    """Use the post-prediction state itself as the word-entry coordinates."""
    n = len(z)
    return [AD.independent(z[i].val, i, n) for i in range(n)]


def _mode(mode: str, path: Path, domain: dict, source_node_index: int, cbox):
    CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    nodes = NODES.build()
    src = NODES.h18_source_cell(source_node_index, nodes)
    F, Q, _Rstep, _ba = CAND._transition_and_Q(mode, src, domain)

    # Build the physical source covariance family, then move the tile boundary
    # through one real prediction.  We never invert the broad goLive box.
    Ppre = CAND._initial_covariance(mode, src, path)
    Pm = H._psd_tighten(
        matrix_add(matrix_mul(matrix_mul(F, Ppre), matrix_transpose(F)), Q)
    )

    zpre = D._initial_ad(mode, domain, cbox)
    z = _rebase_postprediction(D._prediction(mode, zpre, F))

    Mupper, metric_meta = _entry_information_upper(mode)
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
    ops = []

    # The current state is sample-0 post-prediction/pre-measurement.  At the end
    # perform one extra prediction so repeated tiles meet on the same boundary
    # class.  Prediction dissipation is nonnegative and deliberately omitted
    # from signed_word, hence the accumulated measurement form remains a lower
    # bound on total word dissipation.
    for k in range(samples):
        if k in schedule["S_steps"]:
            r = [z[12 + i] for i in range(3)]
            eta = [AD.constant(0.0, n) for _ in range(3)]
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Hs, RS, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "S", **meta})

        if k in schedule["vector_steps"]:
            r = D._exact_acc_residual(mode, z, force)
            eta = D._eta(r, D._linear_residual(Ha, z))
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Ha, Racc, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "acc", **meta})

            r = D._exact_mag_residual(z, mag)
            eta = D._eta(r, D._linear_residual(Hm, z))
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Hm, Rmag, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "mag", **meta})

        # Advance to the next post-prediction boundary, including after the
        # final measurement sample.  The state derivatives remain expressed in
        # the original post-prediction word-entry coordinates.
        Pm = H._psd_tighten(
            matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q)
        )
        z = D._prediction(mode, z, F)

    mu = D._generalized_margin(signed_word, Mupper)
    return {
        "dimension": n,
        "source_node_index": source_node_index,
        "entry_cayley_box": cbox,
        "word_boundary": "POST_PREDICTION_PRE_MEASUREMENT",
        "samples": samples,
        "schedule": schedule,
        "operations": ops,
        "operation_count": len(ops),
        "entry_information_metric_upper": metric_meta,
        "goLive_covariance_inverse_used": False,
        "P3_relative_Riccati_margin_used": False,
        "prediction_dissipation_credited": False,
        "prediction_dissipation_omission_is_conservative": True,
        "signed_word_generalized_margin_design": mu,
        "rho_homogeneous_design_upper": math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0,
        "signed_word_form_interval_ldlt_positive": mu > 0.0,
        "K_interval_matrix_materialized": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          ball_inflation: float = 1.5) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    words = WORDS.build(path)
    failures = [f"word: {x}" for x in WORDS.validate(words)]
    q = 2.0 * math.tan(0.80 / 2.0)
    covers = CAND._ball_box_cover(q, max_box_norm_factor=ball_inflation)
    cbox = covers[0]
    modes = {}
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode(mode, path, domain, int(source_node_index), cbox)
        except Exception as exc:
            failures.append(f"{mode}: {type(exc).__name__}: {exc}")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_POSTPREDICTION_JOINT_JOSEPH_WORD_DISSIPATION_DESIGN",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "design_only_not_theorem_promotion": True,
        "P1_changed": False,
        "P2_phased_graph_required_for_promotion": True,
        "P3_linear_foundation_retained": True,
        "P3_delta_used_as_physical_basin": False,
        "joint_P_H_R_r_used_through_innovation_solve": True,
        "K_interval_matrix_materialized": False,
        "signed_Joseph_identity_used": "r^T S^-1 r - eta^T R^-1 eta",
        "directional_forms_accumulated_before_scalarization": True,
        "endpoint_scalarization_only": True,
        "goLive_covariance_inverse_used": False,
        "postprediction_process_metric_upper_used": True,
        "outer_angle_rad": 0.80,
        "outer_ball_cover_total": len(covers),
        "outer_ball_cells_checked": 1 if modes else 0,
        "source_node_index": int(source_node_index),
        "phased_P2_paths_checked": False,
        "optional_accepted_branch_family_checked": False,
        "full_source_vector_family_checked": False,
        "affine_b_term_established_here": False,
        "modes": modes,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d: dict):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "design_only_not_theorem_promotion",
        "P2_phased_graph_required_for_promotion", "P3_linear_foundation_retained",
        "joint_P_H_R_r_used_through_innovation_solve",
        "directional_forms_accumulated_before_scalarization", "endpoint_scalarization_only",
        "postprediction_process_metric_upper_used",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "P1_changed",
        "P3_delta_used_as_physical_basin", "K_interval_matrix_materialized",
        "goLive_covariance_inverse_used", "phased_P2_paths_checked",
        "optional_accepted_branch_family_checked", "full_source_vector_family_checked",
        "affine_b_term_established_here", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED", "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_node_index=a.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu_design": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho_design": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "ops": d.get("modes", {}).get(m, {}).get("operation_count"),
            "metric": d.get("modes", {}).get(m, {}).get("entry_information_metric_upper"),
            "last_op": (d.get("modes", {}).get(m, {}).get("operations") or [None])[-1],
        } for m in ("H", "A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
