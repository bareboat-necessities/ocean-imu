#!/usr/bin/env python3
"""Attach the canonical P3 comparison metric to the P4 source language.

P4 must not rebuild a convenient one-node covariance metric unrelated to the
source history that canonical P3 certified.  This producer joins three retained
objects without changing any theorem gate:

* the exact P2-V1 full-word Pareto history quotient;
* the source-complete 800 x 26 post-measurement translation covariance floor;
* the H/A conditional attitude/bias floor and same-history covariance ceilings
  used by the current P3 precision-block join.

The endpoint is chosen at the second accepted magnetic packet of a required
vector-PE event.  The P3 covariance word explicitly reserves one full PE
recurrence after the translation observation horizon, so such an endpoint is
reached within the certified word bound.  Translation is already post-acc/S at
that sample; magnetometer has zero translation Jacobian on the zero-lever-arm
branch.  The fresh H/A process block is then attenuated by exactly the endpoint
accelerometer/magnetometer pair as in canonical P3.

For each finite-clock source node and sample phase, the translation theorem has

    P_T|B >= rho_z D_h^2,
    D_h = diag(h, h^2, h^3, 1)

per physical axis.  The H/A producer supplies a diagonal conditional lower for
B=(theta,b_g[,b_a]).  The same precision-block theorem used by P3 gives

    P_full >= 1/2 diag(P_T|B, P_B|T).

Consequently a rigorous diagonal *upper* on the endpoint information metric is
obtained simply by reciprocal of that full covariance lower.  No matrix inverse
of an upper covariance bound is used.

The JSON representation is compressed: each of the 800 endpoint rows stores 26
translation identity floors plus two same-history upper envelopes (stage
boundary and positive phase).  A P4 word consumer can reconstruct the H=18 or
A=21 diagonal metric from the documented formulas without a 20,800-row tensor.
The absorbing frozen-clock branch uses the already certified global frozen
translation relative margin with each held node's positive-phase same-history
upper, so it is covered for arbitrary hold duration.

This is metric attachment only.  It cannot promote P4; the canonical P4 gate
still requires a source-complete nonlinear whole-word rho_H<1 and rho_A<1 and a
canonical P3 PASS at the unchanged 1e-18 gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p2_correlation_path_memory as CORR
import ou3_p3_p2_v1_full_state_join as JOIN
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_p3_p2_v1_stage_phase_translation as TRANS
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
JOIN_FACTOR = 0.5
PHASES = tuple(range(26))
STATE_GROUPS = {
    "H": ["theta", "b_g", "v", "p", "S", "a_w"],
    "A": ["theta", "b_g", "v", "p", "S", "a_w", "b_a"],
}
STATE_ORDER = {
    "H": ["theta_x", "theta_y", "theta_z", "b_gx", "b_gy", "b_gz",
          "v_x", "v_y", "v_z", "p_x", "p_y", "p_z",
          "S_x", "S_y", "S_z", "a_wx", "a_wy", "a_wz"],
    "A": ["theta_x", "theta_y", "theta_z", "b_gx", "b_gy", "b_gz",
          "v_x", "v_y", "v_z", "p_x", "p_y", "p_z",
          "S_x", "S_y", "S_z", "a_wx", "a_wy", "a_wz",
          "b_ax", "b_ay", "b_az"],
}


def _finite_positive(x) -> bool:
    return (
        not isinstance(x, bool)
        and isinstance(x, (int, float))
        and math.isfinite(float(x))
        and float(x) > 0.0
    )


def _down_mul(a: float, b: float) -> float:
    return BASE.down(float(a) * float(b))


def _up_recip(x: float) -> float:
    if not _finite_positive(x):
        raise ValueError("strict positive covariance lower required before reciprocal")
    return BASE.up(1.0 / float(x))


def _translation_physical_lower_groups(rho_z: float, h: float) -> dict:
    """Map Pz>=rho I from z=[v/h,p/h^2,S/h^3,a_w] to physical units."""
    if not (_finite_positive(rho_z) and _finite_positive(h)):
        raise ValueError("positive finite scaled translation floor and dt required")
    h2 = BASE.down(h * h)
    h4 = BASE.down(h2 * h2)
    h6 = BASE.down(h4 * h2)
    return {
        "v": _down_mul(rho_z, h2),
        "p": _down_mul(rho_z, h4),
        "S": _down_mul(rho_z, h6),
        "a_w": BASE.down(rho_z),
    }


def _bias_conditional_lower_groups(mode: str, cond: dict) -> dict:
    out = {
        "theta": float(cond["attitude_conditional_posterior_lower"]),
        "b_g": float(cond["gyro_bias_conditional_posterior_lower"]),
    }
    if mode == "A":
        out["b_a"] = float(cond["accel_bias_conditional_posterior_lower"])
    if any(not _finite_positive(x) for x in out.values()):
        raise RuntimeError(f"{mode} conditional bias floor lost strict positivity")
    return out


def _finite_mode_metric(mode: str, rho_z: float, h: float, cond: dict) -> dict:
    """Full endpoint covariance lower and information upper by P3 block join."""
    groups = {}
    groups.update(_bias_conditional_lower_groups(mode, cond))
    groups.update(_translation_physical_lower_groups(rho_z, h))
    ordered = {
        name: BASE.down(JOIN_FACTOR * float(groups[name]))
        for name in STATE_GROUPS[mode]
    }
    if any(not _finite_positive(x) for x in ordered.values()):
        raise RuntimeError(f"{mode} full covariance lower lost strict positivity")
    info = {name: _up_recip(ordered[name]) for name in STATE_GROUPS[mode]}
    return {
        "covariance_lower_group_diagonal": ordered,
        "information_metric_upper_group_diagonal": info,
        "each_group_repeated_coordinates": 3,
        "dimension": 18 if mode == "H" else 21,
    }


def _frozen_mode_metric(mode: str, relative_margin: float,
                        translation_upper: list[float], cond: dict) -> dict:
    """Arbitrary held-source tail metric from P_T >= delta diag(U_T)."""
    if not _finite_positive(relative_margin):
        raise ValueError("strict frozen translation relative margin required")
    if len(translation_upper) != 4 or any(not _finite_positive(x) for x in translation_upper):
        raise ValueError("positive four-group translation upper required")
    trans = {
        name: BASE.down(float(relative_margin) * float(value))
        for name, value in zip(("v", "p", "S", "a_w"), translation_upper)
    }
    groups = {}
    groups.update(_bias_conditional_lower_groups(mode, cond))
    groups.update(trans)
    ordered = {
        name: BASE.down(JOIN_FACTOR * float(groups[name]))
        for name in STATE_GROUPS[mode]
    }
    info = {name: _up_recip(ordered[name]) for name in STATE_GROUPS[mode]}
    return {
        "covariance_lower_group_diagonal": ordered,
        "information_metric_upper_group_diagonal": info,
        "each_group_repeated_coordinates": 3,
        "dimension": 18 if mode == "H" else 21,
    }


def _bias_upper_envelope(fr: dict, endpoint: int, positive_phase: bool,
                         mode: str, domain: dict, blocks: dict) -> dict:
    labels = JOIN._mapped_labels(fr, int(endpoint), bool(positive_phase))
    if not labels:
        raise RuntimeError("P4 metric attachment lost same-history P2 labels")
    max_theta = 0.0
    max_bg = 0.0
    max_ba = 0.0
    limiting = {"theta": None, "b_g": None, "b_a": None}
    for label in labels:
        summary = HIST.label_summary(label, fr)
        row = JOIN._history_bias_upper(summary, mode, domain, blocks)
        if float(row["theta_covariance_upper"]) > max_theta:
            max_theta = float(row["theta_covariance_upper"])
            limiting["theta"] = list(label)
        if float(row["gyro_bias_covariance_upper"]) > max_bg:
            max_bg = float(row["gyro_bias_covariance_upper"])
            limiting["b_g"] = list(label)
        if mode == "A" and float(row["accel_bias_covariance_upper"]) > max_ba:
            max_ba = float(row["accel_bias_covariance_upper"])
            limiting["b_a"] = list(label)
    out = {
        "phase_pareto_labels": len(labels),
        "theta_covariance_upper": max_theta,
        "gyro_bias_covariance_upper": max_bg,
        "accel_bias_covariance_upper": None if mode == "H" else max_ba,
        "limiting_history_labels": limiting,
        "same_history_rows_evaluated_before_envelope": True,
        "raw_tuner_cartesian_extrema_used": False,
    }
    vals = [max_theta, max_bg] + ([] if mode == "H" else [max_ba])
    if any(not _finite_positive(x) for x in vals):
        raise RuntimeError(f"{mode} same-history bias upper lost positivity")
    return out


def _translation_rows(translation: dict) -> dict[tuple[int, int], dict]:
    rows = translation.get("endpoint_phase_rows", [])
    out = {}
    for row in rows:
        key = (int(row["source_node"]), int(row["phase_samples_after_stage_boundary"]))
        if key in out:
            raise RuntimeError(f"duplicate translation endpoint phase {key}")
        out[key] = row
    expected = 800 * len(PHASES)
    if len(out) != expected:
        raise RuntimeError(f"translation candidate has {len(out)} endpoint phases, expected {expected}")
    return out


def _endpoint_envelope(endpoint: int, phase_class: int, fr: dict, row: dict,
                       domain: dict, blocks: dict) -> dict:
    positive = bool(phase_class)
    history = HIST.endpoint_phase_upper(int(endpoint), 1 if positive else 0, fr)
    upper = list(map(float, row["Sigma_translation_diagonal_upper_envelope"]))
    expected_upper = list(map(float, history["Sigma_translation_diagonal_upper_envelope"]))
    if upper != expected_upper:
        raise RuntimeError("translation artifact differs from same-history endpoint envelope")
    if int(row["same_history_covariance_labels"]) != int(history["phase_pareto_labels"]):
        raise RuntimeError("translation artifact label count differs from P2 history quotient")
    return {
        "phase_class": "positive_1_25" if positive else "stage_boundary_0",
        "phase_pareto_labels": int(history["phase_pareto_labels"]),
        "translation_covariance_upper_groups": {
            name: value for name, value in zip(("v", "p", "S", "a_w"), upper)
        },
        "H_bias_covariance_upper": _bias_upper_envelope(
            fr, endpoint, positive, "H", domain, blocks
        ),
        "A_bias_covariance_upper": _bias_upper_envelope(
            fr, endpoint, positive, "A", domain, blocks
        ),
        "same_history_envelope_before_source_endpoint_hull": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, translation_candidate: dict | None = None,
          full_state_candidate: dict | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 P3-metric attachment must not be trajectory fitted")
    if float(BASE.MIN_USEFUL_DELTA) != 1.0e-18:
        raise RuntimeError("canonical P3 useful gate changed")

    translation = TRANS.build(path) if translation_candidate is None else translation_candidate
    tf = TRANS.validate(translation)
    if tf:
        raise RuntimeError(f"translation P3 candidate invalid: {tf}")
    full = JOIN.build(path, translation) if full_state_candidate is None else full_state_candidate
    jf = JOIN.validate(full)
    if jf:
        raise RuntimeError(f"full-state P3 candidate invalid: {jf}")
    if full.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("full-state P3 candidate lost frozen P2-V1 binding")
    if translation.get("full_word_frontier_digest_sha256") is None:
        raise RuntimeError("translation candidate lost P2 history frontier digest")

    fr = HIST.frontier_runtime(path)
    if fr["frontier_digest_sha256"] != translation["full_word_frontier_digest_sha256"]:
        raise RuntimeError("P4 metric attachment history frontier differs from canonical P3")
    rows = _translation_rows(translation)
    blocks = JOIN._common_blocks(domain)
    cond = {m: JOIN._conditional_bias_floor(m, domain, blocks) for m in ("H", "A")}
    h = float(fr["rt"]["clock"]["dt_binary32_s"])

    endpoints = []
    finite_min = {m: {g: math.inf for g in STATE_GROUPS[m]} for m in ("H", "A")}
    for endpoint in fr["endpoint_nodes"]:
        t = int(endpoint)
        boundary_row = rows[(t, 0)]
        positive_row = rows[(t, 1)]
        boundary_env = _endpoint_envelope(t, 0, fr, boundary_row, domain, blocks)
        positive_env = _endpoint_envelope(t, 1, fr, positive_row, domain, blocks)

        phase_rhos = []
        phase_metrics = {"H": [], "A": []}
        for phase in PHASES:
            row = rows[(t, int(phase))]
            rho = float(row["phase_z_identity_floor_lower"])
            if not _finite_positive(rho):
                raise RuntimeError(f"translation metric floor is nonpositive at source {t}, phase {phase}")
            expected_env = boundary_env if phase == 0 else positive_env
            row_upper = {
                name: float(value)
                for name, value in zip(
                    ("v", "p", "S", "a_w"),
                    row["Sigma_translation_diagonal_upper_envelope"],
                )
            }
            if row_upper != expected_env["translation_covariance_upper_groups"]:
                raise RuntimeError("positive source phases do not share the retained same-history upper")
            phase_rhos.append(rho)
            for mode in ("H", "A"):
                metric = _finite_mode_metric(mode, rho, h, cond[mode])
                phase_metrics[mode].append(metric["information_metric_upper_group_diagonal"])
                for group, value in metric["covariance_lower_group_diagonal"].items():
                    finite_min[mode][group] = min(finite_min[mode][group], float(value))

        endpoints.append({
            "source_node": t,
            "boundary_history_envelope": boundary_env,
            "positive_phase_history_envelope": positive_env,
            "phase_z_identity_floor_lower": phase_rhos,
            "finite_phase_information_metric_upper_group_diagonal": phase_metrics,
        })

    frozen = translation.get("frozen_clock_branch", {})
    frozen_delta = float(frozen.get("worst_margin_lower", 0.0))
    if frozen.get("included") is not True or frozen.get("clock_stagnation_verified") is not True:
        raise RuntimeError("translation candidate lost frozen-clock theorem branch")
    if not _finite_positive(frozen_delta):
        raise RuntimeError("frozen-clock translation branch has no strict relative floor")

    # Each held source keeps the positive-phase same-history covariance upper.
    # Since every frozen row has delta_t >= global frozen_delta, P_T is at least
    # frozen_delta times that endpoint diagonal upper for arbitrary hold time.
    frozen_rows = []
    frozen_min = {m: {g: math.inf for g in STATE_GROUPS[m]} for m in ("H", "A")}
    for endpoint in endpoints:
        upper_groups = endpoint["positive_phase_history_envelope"]["translation_covariance_upper_groups"]
        upper = [upper_groups[k] for k in ("v", "p", "S", "a_w")]
        modes = {}
        for mode in ("H", "A"):
            metric = _frozen_mode_metric(mode, frozen_delta, upper, cond[mode])
            modes[mode] = metric["information_metric_upper_group_diagonal"]
            for group, value in metric["covariance_lower_group_diagonal"].items():
                frozen_min[mode][group] = min(frozen_min[mode][group], float(value))
        frozen_rows.append({
            "source_node": endpoint["source_node"],
            "information_metric_upper_group_diagonal": modes,
        })

    global_lower = {}
    global_info = {}
    for mode in ("H", "A"):
        lower = {
            g: BASE.down(min(float(finite_min[mode][g]), float(frozen_min[mode][g])))
            for g in STATE_GROUPS[mode]
        }
        if any(not _finite_positive(x) for x in lower.values()):
            raise RuntimeError(f"{mode} global source/phase covariance lower lost positivity")
        global_lower[mode] = lower
        global_info[mode] = {g: _up_recip(lower[g]) for g in STATE_GROUPS[mode]}

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CANONICAL_P3_SOURCE_PHASE_METRIC_ATTACHMENT",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "canonical_P3_useful_gate": BASE.MIN_USEFUL_DELTA,
        "canonical_P3_candidate_qualification": full.get("qualification"),
        "canonical_P3_candidate_numeric_pass_observed": full.get("P3_PRODUCER_NUMERIC_PASS") is True,
        "canonical_P3_pass_required_before_P4_theorem_consumption": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "full_word_frontier_digest_sha256": fr["frontier_digest_sha256"],
        "same_history_P3_frontier_consumed": True,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "metric_boundary": "SECOND_ACCEPTED_MAG_PACKET_OF_REQUIRED_VECTOR_PE_EVENT_POST_MEASUREMENT",
        "endpoint_vector_PE_recurrence_reserved_by_P3_word": True,
        "translation_post_acc_S_at_metric_boundary": True,
        "magnetometer_translation_jacobian_zero_on_declared_branch": True,
        "H_A_fresh_process_floor_at_same_endpoint_vector_packet": True,
        "precision_block_join_factor": JOIN_FACTOR,
        "precision_block_theorem_same_as_canonical_P3": True,
        "information_upper_derived_only_from_covariance_lower": True,
        "upper_covariance_never_inverted_to_claim_information_upper": True,
        "state_group_order": STATE_GROUPS,
        "state_order": STATE_ORDER,
        "translation_scaling": "z=[v/h,p/h^2,S/h^3,a_w]",
        "configured_dt_s": h,
        "finite_source_nodes": len(endpoints),
        "finite_phase_count_per_source": len(PHASES),
        "finite_source_phase_classes": len(endpoints) * len(PHASES),
        "endpoint_rows": endpoints,
        "conditional_bias_floor": cond,
        "frozen_clock": {
            "included": True,
            "absorbing_hold_arbitrary_duration_covered": True,
            "translation_relative_margin_lower_global": frozen_delta,
            "source_rows": frozen_rows,
        },
        "global_source_phase_covariance_lower_group_diagonal": global_lower,
        "global_source_phase_information_metric_upper_group_diagonal": global_info,
        "metric_attachment_structurally_complete": True,
        "complete_nonlinear_word_evaluated_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "propagate the full H=18/A=21 nonlinear word on this exact endpoint/source-phase metric family, including source-edge metric transport and the 0.8-rad signed vector remainder; submit rho_H/rho_A to the frozen canonical P4 gate"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_CANONICAL_P3_SOURCE_PHASE_METRIC_ATTACHMENT":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch",
        "canonical_P3_pass_required_before_P4_theorem_consumption",
        "same_history_P3_frontier_consumed",
        "endpoint_vector_PE_recurrence_reserved_by_P3_word",
        "translation_post_acc_S_at_metric_boundary",
        "magnetometer_translation_jacobian_zero_on_declared_branch",
        "H_A_fresh_process_floor_at_same_endpoint_vector_packet",
        "precision_block_theorem_same_as_canonical_P3",
        "information_upper_derived_only_from_covariance_lower",
        "upper_covariance_never_inverted_to_claim_information_upper",
        "metric_attachment_structurally_complete",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "independent_cartesian_tau_sigma_R_S_extrema_used",
        "complete_nonlinear_word_evaluated_here", "P4_USABLE_CERTIFICATE_PROMOTED",
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("P4 metric attachment lost frozen P2-V1 binding")
    if float(d.get("canonical_P3_useful_gate", math.nan)) != 1.0e-18:
        f.append("canonical P3 useful gate changed at P4 metric attachment")
    if float(d.get("precision_block_join_factor", math.nan)) != JOIN_FACTOR:
        f.append("P3 precision-block join factor changed")
    if d.get("state_group_order") != STATE_GROUPS or d.get("state_order") != STATE_ORDER:
        f.append("H/A state order changed")
    if int(d.get("finite_source_nodes", 0)) != 800:
        f.append("metric attachment does not cover 800 P2 source nodes")
    if int(d.get("finite_phase_count_per_source", 0)) != 26:
        f.append("metric attachment does not cover phases 0..25")
    if int(d.get("finite_source_phase_classes", 0)) != 800 * 26:
        f.append("metric attachment does not cover 800x26 finite source phases")
    endpoints = d.get("endpoint_rows", [])
    if len(endpoints) != 800:
        f.append("endpoint metric row count is not 800")
    else:
        seen = set()
        for row in endpoints:
            source = row.get("source_node")
            if not isinstance(source, int) or source in seen:
                f.append("endpoint metric source ids are invalid or duplicated")
                break
            seen.add(source)
            rhos = row.get("phase_z_identity_floor_lower", [])
            if len(rhos) != 26 or any(not _finite_positive(x) for x in rhos):
                f.append(f"source {source}: invalid finite phase covariance floors")
                break
            for env_key in ("boundary_history_envelope", "positive_phase_history_envelope"):
                env = row.get(env_key, {})
                u = env.get("translation_covariance_upper_groups", {})
                if list(u.keys()) != ["v", "p", "S", "a_w"] or any(not _finite_positive(x) for x in u.values()):
                    f.append(f"source {source}: invalid {env_key} translation upper")
                    break
    frozen = d.get("frozen_clock", {})
    if frozen.get("included") is not True or frozen.get("absorbing_hold_arbitrary_duration_covered") is not True:
        f.append("frozen-clock metric branch is incomplete")
    if not _finite_positive(frozen.get("translation_relative_margin_lower_global")):
        f.append("frozen-clock translation metric floor is not strict")
    if len(frozen.get("source_rows", [])) != 800:
        f.append("frozen-clock metric does not cover all 800 held sources")
    for mode in ("H", "A"):
        lower = d.get("global_source_phase_covariance_lower_group_diagonal", {}).get(mode, {})
        info = d.get("global_source_phase_information_metric_upper_group_diagonal", {}).get(mode, {})
        if list(lower.keys()) != STATE_GROUPS[mode] or any(not _finite_positive(x) for x in lower.values()):
            f.append(f"{mode}: invalid global covariance lower metric")
        if list(info.keys()) != STATE_GROUPS[mode] or any(not _finite_positive(x) for x in info.values()):
            f.append(f"{mode}: invalid global information upper metric")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--translation-candidate", type=Path)
    ap.add_argument("--full-state-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    tc = json.loads(a.translation_candidate.read_text(encoding="utf-8")) if a.translation_candidate else None
    fc = json.loads(a.full_state_candidate.read_text(encoding="utf-8")) if a.full_state_candidate else None
    d = build(a.domain, translation_candidate=tc, full_state_candidate=fc)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "finite_source_phase_classes": d["finite_source_phase_classes"],
        "frontier_digest": d["full_word_frontier_digest_sha256"],
        "P3_numeric_pass_observed": d["canonical_P3_candidate_numeric_pass_observed"],
        "global_information_metric_upper": d["global_source_phase_information_metric_upper_group_diagonal"],
        "frozen_translation_margin": d["frozen_clock"]["translation_relative_margin_lower_global"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
