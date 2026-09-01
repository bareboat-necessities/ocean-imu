#!/usr/bin/env python3
"""Path-dependent complete-word differential-contraction route for OU-III P4.

This is an independent alternative to operation-by-operation nonlinear budget
accounting.  It treats one source-complete fixed-mode Live word as the atomic
map and asks for a direct differential contraction certificate

    J_w(x)^T M_h J_w(x) <= gamma^2 M_g,   gamma < 1,

for every reachable source-word edge g->h and every state in the declared
finite-angle domain.  M_g and M_h are source covariance-information metrics
with one mode-global normalization; the endpoint metric is selected by the
actual reachable source node rather than by an independent worst-case choice.

The route deliberately does not require individual packet full rank,
per-operation sector invariance, an N-times-global Lipschitz remainder, or a
translation/nontranslation Schur split.  The tiny source-uniform P3 delta is a
linear-origin strictness diagnostic only, never a nonlinear radius.

A future numerical producer must outward-enclose the full H=18/A=21
Jacobian/generalized Jacobian of the exact shipping complete-word return map,
including sequential quaternion resets and the A-mode accelerometer-bias
projection.  Candidate metadata is accepted only when it is bound to the exact
proof-operating-domain file and the exact 0.80-rad outer geometry.  This route
contract never promotes P4 by accepting metadata alone.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path

import ou3_p4_exact_word_map as WORDMAP
import ou3_p4_node_metrics as METRIC
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_path_reachability as PATH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DESIGN_OUTER_RAD = 0.80
P5_ENTRANCE_DEG = 45.0
JACOBIAN_QUALIFICATION = "OU3_P4_OUTWARD_COMPLETE_WORD_GENERALIZED_JACOBIAN"


def _decimal_linear_gap(delta: float) -> dict:
    """Return a high-precision view of 1-sqrt(1-delta)."""
    d = Decimal(str(float(delta)))
    if not (Decimal(0) < d < Decimal(1)):
        raise ValueError("strict delta in (0,1) required")
    getcontext().prec = 90
    gamma = (Decimal(1) - d).sqrt()
    gap = Decimal(1) - gamma
    return {
        "delta": str(d),
        "sqrt_one_minus_delta": format(gamma, ".80E"),
        "one_minus_sqrt_one_minus_delta": format(gap, ".80E"),
    }


def _mode_route(mode: str, metric: dict) -> dict:
    """Build the fixed-mode metric and direct Jacobian acceptance contract."""
    m = metric["modes"][mode]
    dim = int(m["dimension"])
    delta = float(m["P3_word_endpoint_margin_lower"])
    return {
        "dimension": dim,
        "start_metric": "M_g=s_mode*Sigma_KF(g)^-1 in exact Cayley/linear coordinates",
        "endpoint_metric": "M_h=s_mode*Sigma_KF(h)^-1 at the actual reachable endpoint source node",
        "same_mode_global_scale_on_all_nodes": m["same_scale_on_every_source_node_in_mode"],
        "full_attitude_linear_cross_terms_retained": m["full_attitude_linear_cross_terms_retained"],
        "block_diagonal_metric_used": m["block_diagonal_metric_used"],
        "P3_linear_origin_margin_lower": delta,
        "P3_linear_origin_contraction_gap_diagnostic": _decimal_linear_gap(delta),
        "P3_delta_used_as_nonlinear_radius": False,
        "P3_delta_used_as_whole_word_jacobian_acceptance_threshold": False,
        "acceptance_matrix_inequality": "J_w(x)^T M_h J_w(x) <= gamma_mode^2 M_g",
        "equivalent_whitened_test": "||L_h J_w(x) L_g^-1||_2 <= gamma_mode < 1, L_g^T L_g=M_g",
        "required_gamma_predicate": "< 1",
        "complete_full_state_matrix_test": True,
        "translation_nontranslation_schur_split_required": False,
        "per_packet_full_rank_required": False,
        "per_operation_contraction_required": False,
        "finite_angle_whole_word_generalized_jacobian_required": True,
        "A_mode_projection_generalized_jacobian_required": mode == "A",
    }


def _candidate_status(candidate: dict | None, route: dict) -> dict:
    """Validate only the metadata contract of a future rigorous Jacobian result.

    The actual matrices must be owned and promoted by their numerical producer.
    This consumer may return CONTRACT_ACCEPTED but never establishes P4.
    """
    if candidate is None:
        return {
            "provided": False,
            "contract_accepted": False,
            "all_modes_strict": False,
            "reasons": ["no outward complete-word generalized-Jacobian certificate supplied"],
        }

    reasons: list[str] = []
    if candidate.get("qualification") != JACOBIAN_QUALIFICATION:
        reasons.append("wrong Jacobian certificate qualification")
    if candidate.get("source_only") is not True:
        reasons.append("Jacobian certificate is not source-only")
    if candidate.get("trajectory_replay_used") is not False:
        reasons.append("Jacobian certificate uses trajectory replay")
    if candidate.get("outward_validated") is not True:
        reasons.append("Jacobian certificate is not outward validated")
    if candidate.get("exact_shipping_complete_word_map") is not True:
        reasons.append("Jacobian certificate is not bound to exact shipping complete-word map")
    if candidate.get("per_operation_sector_invariance_required") is not False:
        reasons.append("Jacobian certificate reintroduced per-operation sector invariance")
    if candidate.get("N_times_global_defect_used") is not False:
        reasons.append("Jacobian certificate reintroduced N-times-global-defect accounting")
    if candidate.get("P3_delta_used_as_nonlinear_radius") is not False:
        reasons.append("Jacobian certificate reintroduced P3 delta as nonlinear radius")

    # Bind an accepted candidate to the exact physical/full-state theorem domain.
    # The content hash covers every startup/live/full-state bound in the declared
    # operating-domain JSON; the explicit angle/q fields make the finite-angle
    # coverage human-auditable and reject a rigorous result on a smaller ball.
    if candidate.get("proof_operating_domain_sha256") != route.get("proof_operating_domain_sha256"):
        reasons.append("candidate proof-operating-domain hash does not match route domain")
    candidate_angle = candidate.get("certified_outer_angle_rad")
    if isinstance(candidate_angle, bool) or not isinstance(candidate_angle, (int, float)) \
            or not math.isfinite(float(candidate_angle)) \
            or float(candidate_angle) != float(route["outer_geometry_angle_rad"]):
        reasons.append("candidate outer-angle domain is not exactly the required 0.80 rad")
    candidate_q = candidate.get("certified_outer_cayley_norm_upper")
    if isinstance(candidate_q, bool) or not isinstance(candidate_q, (int, float)) \
            or not math.isfinite(float(candidate_q)) \
            or float(candidate_q) < float(route["outer_geometry_cayley_norm_upper"]):
        reasons.append("candidate Cayley domain does not cover the required outer ball")
    if candidate.get("full_state_domain_exact_match") is not True:
        reasons.append("candidate does not declare an exact full-state theorem-domain match")

    if int(candidate.get("source_partition_state_count", -1)) != int(route["source_graph"]["partition_state_count"]):
        reasons.append("source partition state count mismatch")
    if int(candidate.get("source_transition_edge_count", -1)) != int(route["source_graph"]["transition_edge_count"]):
        reasons.append("source transition edge count mismatch")
    if int(candidate.get("recurrent_state_count", -1)) != int(route["source_graph"]["recurrent_state_count"]):
        reasons.append("recurrent source state count mismatch")
    if candidate.get("all_required_reachable_word_edges_checked") is not True:
        reasons.append("not all required reachable complete-word edges were checked")
    if candidate.get("sequential_reset_jacobians_included") is not True:
        reasons.append("sequential reset Jacobians were not included")
    if candidate.get("accepted_rejected_not_due_branches_covered") is not True:
        reasons.append("accepted/rejected/not-due branches are incomplete")

    strict = True
    for mode in ("H", "A"):
        row = candidate.get("modes", {}).get(mode, {})
        expected_dim = int(route["modes"][mode]["dimension"])
        if int(row.get("dimension", -1)) != expected_dim:
            reasons.append(f"{mode}: wrong Jacobian dimension")
            strict = False
        if row.get("endpoint_metric_uses_actual_reachable_source_node") is not True:
            reasons.append(f"{mode}: endpoint metric is not source-node correlated")
            strict = False
        if row.get("full_state_cross_terms_retained") is not True:
            reasons.append(f"{mode}: full-state Jacobian cross terms were dropped")
            strict = False
        gamma = row.get("max_whitened_generalized_jacobian_norm_upper")
        if isinstance(gamma, bool) or not isinstance(gamma, (int, float)) \
                or not math.isfinite(float(gamma)) or not (0.0 <= float(gamma) < 1.0):
            reasons.append(f"{mode}: strict whole-word contraction gamma is not certified below one")
            strict = False
        if mode == "A" and row.get("accelerometer_bias_projection_generalized_jacobian_included") is not True:
            reasons.append("A: accelerometer-bias projection generalized Jacobian missing")
            strict = False

    accepted = not reasons and strict
    return {
        "provided": True,
        "contract_accepted": accepted,
        "all_modes_strict": strict,
        "reasons": reasons,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, jacobian_candidate: dict | None = None) -> dict:
    """Build the independent whole-word differential-contraction route."""
    path = Path(domain_path).resolve()
    raw_domain = path.read_bytes()
    domain = json.loads(raw_domain.decode("utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("path-dependent P4 route must not be trajectory fitted")

    words = WORDMAP.build(path)
    metric = METRIC.build(path)
    sector = SECTOR.build(path)
    graph = PATH.build(path)

    failures = [f"word-map: {x}" for x in WORDMAP.validate(words)]
    failures += [f"metric: {x}" for x in METRIC.validate(metric)]
    failures += [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"source-path: {x}" for x in PATH.validate(graph)]

    outer = float(sector["design_full_attitude_angle_rad"])
    q_outer = float(sector["design_cayley_norm_upper"])
    if outer != DESIGN_OUTER_RAD:
        failures.append("whole-word route requires the exact 0.80-rad outer geometry")
    if not q_outer < 1.0:
        failures.append("whole-word route left the q<1 Cayley chart")

    entrance = domain.get("initial_filter_entrance", {}).get("attitude", {})
    entrance_deg = float(entrance.get("full_attitude_error_upper_deg", math.nan))
    if entrance_deg != P5_ENTRANCE_DEG:
        failures.append("whole-word route must preserve the 45-deg P5 entrance")

    modes = {mode: _mode_route(mode, metric) for mode in ("H", "A")}
    if any(modes[m]["block_diagonal_metric_used"] is not False for m in modes):
        failures.append("path-dependent route requires full source information metrics")

    route = {
        "schema": SCHEMA,
        "qualification": "OU3_P4_PATH_DEPENDENT_COMPLETE_WORD_DIFFERENTIAL_CONTRACTION_ROUTE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "proof_operating_domain_sha256": hashlib.sha256(raw_domain).hexdigest(),
        "candidate_domain_contract": {
            "proof_operating_domain_hash_exact_match_required": True,
            "full_state_domain_exact_match_required": True,
            "certified_outer_angle_rad_required": DESIGN_OUTER_RAD,
            "certified_outer_cayley_norm_must_cover": q_outer,
        },
        "P5_entrance_angle_deg": entrance_deg,
        "P5_45deg_entrance_preserved": entrance_deg == P5_ENTRANCE_DEG,
        "outer_geometry_angle_rad": outer,
        "outer_geometry_cayley_norm_upper": q_outer,
        "outer_geometry_exact_0p80_rad_required": True,
        "word_horizon_s": float(words["source_word_horizon_s"]),
        "word_samples_upper": int(words["word_samples_upper"]),
        "shipping_operation_order": words["shipping_operation_order"],
        "source_graph": {
            "partition_state_count": int(graph["partition"]["states"]),
            "transition_edge_count": int(graph["transition_edges"]),
            "strongly_connected_components": int(graph["strongly_connected_components"]),
            "recurrent_state_count": int(graph["recurrent_states"]),
            "path_graph_ready": graph["path_graph_ready"],
            "raw_tuner_sigma_subfloor_states_included": graph["raw_tuner_sigma_subfloor_states_included"],
            "RS_target_powf_tightening_used": graph["RS_target_powf_tightening_used"],
        },
        "modes": modes,
        "route_distinctions": {
            "atomic_object": "COMPLETE_SOURCE_WORD_RETURN_MAP",
            "individual_packet_rank_may_be_singular": True,
            "per_packet_full_rank_required": False,
            "per_operation_sector_invariance_required": False,
            "N_times_global_lipschitz_defect_used": False,
            "translation_nontranslation_schur_split_required": False,
            "full_18_21_state_jacobian_test_required": True,
            "source_node_dependent_metric_switching": True,
            "arbitrary_node_rescaling_allowed": False,
            "endpoint_metric_tied_to_actual_reachable_source_node": True,
            "P3_delta_used_as_nonlinear_radius": False,
            "P3_delta_used_only_as_linear_origin_strictness_fallback": True,
        },
        "finite_angle_acceptance": {
            "criterion": "for every reachable word edge g->h and every x in the finite-angle/full-state domain: J_w(x)^T M_h J_w(x) <= gamma^2 M_g with gamma<1",
            "outward_interval_or_validated_Taylor_model_required": True,
            "complete_word_composed_before_scalarization": True,
            "sequential_quaternion_resets_must_be_differentiated": True,
            "A_mode_bias_projection_requires_generalized_jacobian": True,
            "rejected_and_not_due_identity_branches_must_be_covered": True,
            "active_vibration_guard_covered": False,
            "active_vibration_guard_remains_separate_hybrid_source_obligation": True,
        },
        "P4_COMPLETE_WORD_DIFFERENTIAL_CONTRACTION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }
    route["jacobian_candidate"] = _candidate_status(jacobian_candidate, route)
    route["route_contract_pass"] = not failures
    route["next_obligation"] = (
        "build the full source-reachable 18/21-state generalized-Jacobian enclosure of the exact one-word return map, "
        "whiten each edge with the actual source-node information metrics, and directly certify the worst edge norm below one over the exact 0.80-rad/full-state domain; "
        "do not scalarize packet defects before composing the word"
    )
    return route


def validate(d: dict) -> list[str]:
    """Validate fail-closed route semantics without promoting an input candidate."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "P5_45deg_entrance_preserved",
        "outer_geometry_exact_0p80_rad_required", "route_contract_pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "P4_COMPLETE_WORD_DIFFERENTIAL_CONTRACTION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED", "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    domain_hash = d.get("proof_operating_domain_sha256")
    if not isinstance(domain_hash, str) or len(domain_hash) != 64:
        f.append("proof operating-domain hash is missing")
    contract = d.get("candidate_domain_contract", {})
    if contract.get("proof_operating_domain_hash_exact_match_required") is not True:
        f.append("candidate domain hash is not required")
    if contract.get("full_state_domain_exact_match_required") is not True:
        f.append("candidate full-state domain exact match is not required")
    if float(d.get("P5_entrance_angle_deg", 0.0)) != P5_ENTRANCE_DEG:
        f.append("P5 entrance angle changed")
    if float(d.get("outer_geometry_angle_rad", math.nan)) != DESIGN_OUTER_RAD:
        f.append("outer geometry is not exactly 0.80 rad")
    if not float(d.get("outer_geometry_cayley_norm_upper", math.inf)) < 1.0:
        f.append("outer geometry left q<1 chart")

    graph = d.get("source_graph", {})
    if graph.get("path_graph_ready") is not True:
        f.append("source graph is not ready")
    if int(graph.get("partition_state_count", 0)) <= 0 or int(graph.get("transition_edge_count", 0)) <= 0:
        f.append("source graph is empty")
    if int(graph.get("recurrent_state_count", 0)) <= 0:
        f.append("source graph has no recurrent states")
    if graph.get("raw_tuner_sigma_subfloor_states_included") is not True:
        f.append("source graph omitted raw tuner sub-floor sigma states")
    if graph.get("RS_target_powf_tightening_used") is not False:
        f.append("route relies on unqualified powf/sqrtf path tightening")

    distinctions = d.get("route_distinctions", {})
    for key in (
        "individual_packet_rank_may_be_singular", "full_18_21_state_jacobian_test_required",
        "source_node_dependent_metric_switching", "endpoint_metric_tied_to_actual_reachable_source_node",
        "P3_delta_used_only_as_linear_origin_strictness_fallback",
    ):
        if distinctions.get(key) is not True:
            f.append(f"route distinction {key} is not true")
    for key in (
        "per_packet_full_rank_required", "per_operation_sector_invariance_required",
        "N_times_global_lipschitz_defect_used", "translation_nontranslation_schur_split_required",
        "arbitrary_node_rescaling_allowed", "P3_delta_used_as_nonlinear_radius",
    ):
        if distinctions.get(key) is not False:
            f.append(f"route distinction {key} is not false")

    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if int(row.get("dimension", -1)) != dim:
            f.append(f"{mode}: wrong dimension")
        if row.get("same_mode_global_scale_on_all_nodes") is not True:
            f.append(f"{mode}: node-dependent rescaling could manufacture contraction")
        if row.get("full_attitude_linear_cross_terms_retained") is not True:
            f.append(f"{mode}: source information cross terms were dropped")
        if row.get("block_diagonal_metric_used") is not False:
            f.append(f"{mode}: block-diagonal metric reintroduced")
        if row.get("P3_delta_used_as_nonlinear_radius") is not False:
            f.append(f"{mode}: P3 delta reintroduced as nonlinear radius")
        if row.get("translation_nontranslation_schur_split_required") is not False:
            f.append(f"{mode}: route reintroduced Schur split")
        delta = row.get("P3_linear_origin_margin_lower")
        if not isinstance(delta, (int, float)) or isinstance(delta, bool) or not (0.0 < float(delta) < 1.0):
            f.append(f"{mode}: missing strict P3 origin margin")

    finite = d.get("finite_angle_acceptance", {})
    if finite.get("complete_word_composed_before_scalarization") is not True:
        f.append("finite-angle route scalarizes before complete-word composition")
    if finite.get("A_mode_bias_projection_requires_generalized_jacobian") is not True:
        f.append("A-mode projection generalized Jacobian omitted")
    if finite.get("active_vibration_guard_covered") is not False:
        f.append("active vibration guard was silently admitted")
    if finite.get("active_vibration_guard_remains_separate_hybrid_source_obligation") is not True:
        f.append("active vibration guard is not fail-closed")

    if d.get("P4_USABLE_CERTIFICATE_PROMOTED") is not False:
        f.append("route contract prematurely promoted P4")
    return list(dict.fromkeys(f))


def main() -> int:
    """CLI entry point for route diagnostics and optional candidate-schema check."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--jacobian-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    candidate = None
    if a.jacobian_candidate:
        candidate = json.loads(a.jacobian_candidate.read_text(encoding="utf-8"))
    d = build(a.domain.resolve(), candidate)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "route_pass": d["route_contract_pass"],
        "P5_entrance_deg": d["P5_entrance_angle_deg"],
        "outer_rad": d["outer_geometry_angle_rad"],
        "source_graph": d["source_graph"],
        "modes": {
            mode: {
                "dimension": d["modes"][mode]["dimension"],
                "P3_delta": d["modes"][mode]["P3_linear_origin_margin_lower"],
                "P3_gap_diagnostic": d["modes"][mode]["P3_linear_origin_contraction_gap_diagnostic"],
            } for mode in ("H", "A")
        },
        "jacobian_candidate": d["jacobian_candidate"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
