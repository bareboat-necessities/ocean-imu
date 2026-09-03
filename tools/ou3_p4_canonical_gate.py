#!/usr/bin/env python3
"""Single canonical P4 theorem-interface gate for OU-III.

Canonical P4 PASS is deliberately stronger than any retained sector/source
prerequisite.  It requires:

* the unique canonical P3 theorem artifact and unchanged 1e-18 H/A gate;
* the source-generated 0.8-rad Cayley geometry and nonlinear vector-remainder
  sector certificates;
* the source-complete P4 word timing decomposition;
* the exact 800-node P2 source partition and finite-speed sample-clock
  transition refinement, including the absorbing frozen-clock branch; and
* one complete H=18/A=21 nonlinear whole-word dissipation candidate that is
  bound to those exact inputs and proves strict rho_H<1 and rho_A<1.

This module does not create a contraction estimate.  It freezes the theorem
boundary so that sector primitives, source diagnostics, or a linear P3 PASS can
never be relabeled as P4.  A future numerical/analytic complete-word producer
must satisfy this interface without changing the definition.

Strictness representation.  ``rho<1`` is a semantic obligation, not a claim
that ``rho`` is representable below one in binary64.  The certified source-word
margin is far below binary64 epsilon, so ``1-delta/2`` rounds back to exactly
one.  A producer may therefore carry the strict gap exactly in
``one_minus_rho_{H,A}_lower``; the gate then requires that gap to be finite
positive and requires the reported ``rho_{H,A}_upper`` to be consistent with
it.  A producer that omits the gap must still show strictness inside the float
itself, exactly as before.  Either way the gate accepts only a strictly
dissipative complete word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

SCHEMA = 1
QUALIFICATION = "OU3_P4_CANONICAL_COMPLETE_WORD_DISSIPATION"
P3_QUALIFICATION = "OU3_P3_CANONICAL_THEOREM_INTERFACE"
CANDIDATE_QUALIFICATION = "OU3_P4_COMPLETE_NONLINEAR_WORD_DISSIPATION_V1"
USEFUL_P3_GATE = 1.0e-18
REQUIRED_OUTER_ANGLE_RAD = 0.80
H_DIM = 18
A_DIM = 21


def _finite_number(value) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _true(payload: dict, key: str, reasons: list[str], label: str) -> None:
    if payload.get(key) is not True:
        reasons.append(f"{label}: {key} is not true")


def _false(payload: dict, key: str, reasons: list[str], label: str) -> None:
    if payload.get(key) is not False:
        reasons.append(f"{label}: {key} is not false")


def _validated(payload: dict, reasons: list[str], label: str) -> None:
    if payload.get("validation_pass") is not True:
        reasons.append(f"{label}: validation_pass is not true")
    failures = payload.get("validation_failures", [])
    if not isinstance(failures, list) or failures:
        reasons.append(f"{label}: validation failures are nonempty or malformed")


def _check_p3(p3: dict, reasons: list[str]) -> None:
    label = "P3"
    _validated(p3, reasons, label)
    if p3.get("qualification") != P3_QUALIFICATION:
        reasons.append("P3: wrong canonical qualification")
    if p3.get("canonical_definition_frozen") is not True:
        reasons.append("P3: canonical definition is not frozen")
    if p3.get("only_this_module_may_promote_P3_for_P4") is not True:
        reasons.append("P3: promotion authority is not unique")
    if p3.get("useful_gate") != USEFUL_P3_GATE:
        reasons.append("P3: unchanged 1e-18 useful gate is not present")
    if p3.get("P3_CANONICAL_PASS") is not True:
        reasons.append("P3: canonical P3 has not passed")
    if p3.get("P4_MAY_CONSUME_P3") is not True:
        reasons.append("P3: P4 consumption has not been authorized")
    worst = p3.get("worst_H_A_margin")
    if not _finite_number(worst) or float(worst) < USEFUL_P3_GATE:
        reasons.append("P3: worst H/A margin is missing or below 1e-18")
    margins = p3.get("mode_margins", {})
    for mode in ("H", "A"):
        x = margins.get(mode) if isinstance(margins, dict) else None
        if not _finite_number(x) or float(x) < USEFUL_P3_GATE:
            reasons.append(f"P3: {mode} margin is missing or below 1e-18")


def _check_cayley(d: dict, reasons: list[str]) -> None:
    label = "Cayley"
    _validated(d, reasons, label)
    if d.get("qualification") != "OU3_P4_GLOBAL_CAYLEY_SECTOR_GEOMETRY":
        reasons.append("Cayley: wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "declared_filter_entrance_covered",
        "chart_antipode_excluded",
        "usable_sector_geometry_pass",
        "pass",
    ):
        _true(d, key, reasons, label)
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "full_18_21_state_Joseph_word_established_here",
        "signed_EKF_remainder_charged_here",
    ):
        _false(d, key, reasons, label)
    theta = d.get("outer_angle_rad")
    if not _finite_number(theta) or float(theta) < REQUIRED_OUTER_ANGLE_RAD:
        reasons.append("Cayley: certified sector does not cover 0.8 rad")


def _check_remainder(d: dict, reasons: list[str]) -> None:
    label = "Remainder"
    _validated(d, reasons, label)
    if d.get("qualification") != "OU3_P4_GLOBAL_VECTOR_NONLINEAR_REMAINDER_SECTOR":
        reasons.append("Remainder: wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "declared_filter_entrance_covered",
        "accelerometer_bias_cancels_exactly_from_eta",
        "penalties_are_homogeneous_quadratic_not_affine_beta",
        "measurement_covariance_isotropy_required",
        "usable_sector_remainder_primitive_pass",
        "pass",
    ):
        _true(d, key, reasons, label)
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "complete_Joseph_word_established_here",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        _false(d, key, reasons, label)
    theta = d.get("outer_angle_rad")
    if not _finite_number(theta) or float(theta) < REQUIRED_OUTER_ANGLE_RAD:
        reasons.append("Remainder: certified sector does not cover 0.8 rad")


def _check_timing(d: dict, reasons: list[str]) -> None:
    label = "Word timing"
    _validated(d, reasons, label)
    if d.get("qualification") != "OU3_P4_SOURCE_COMPLETE_WORD_TIMING_DECOMPOSITION":
        reasons.append("Word timing: wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "S_firing_times_are_source_intervals_not_fixed_samples",
        "S_residual_exactly_linear_selector",
        "S_nonlinear_eta_identically_zero",
        "S_timing_consumed_by_linear_P3_translation_UCO",
        "nonlinear_timing_obligations_reduce_to_vector_measurements",
        "ready_for_source_complete_nonlinear_remainder_composition",
    ):
        _true(d, key, reasons, label)
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "fixed_minimum_gap_S_schedule_is_source_complete",
        "old_terminal_192_201_cluster_required_for_promotion",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
    ):
        _false(d, key, reasons, label)


def _check_path(d: dict, reasons: list[str]) -> None:
    label = "P2 path"
    _validated(d, reasons, label)
    if d.get("P2_SOURCE_PATH_CERTIFICATE") != "PASS":
        reasons.append("P2 path: source path certificate did not pass")
    _true(d, "path_graph_ready", reasons, label)
    _false(d, "usable_P4_promoted", reasons, label)
    if int(d.get("partition", {}).get("states", 0) or 0) != 800:
        reasons.append("P2 path: physical state partition is not 800 nodes")


def _check_nodes(d: dict, reasons: list[str]) -> None:
    label = "Source nodes"
    _validated(d, reasons, label)
    if d.get("qualification") != "OU3_P2_SOURCE_NODE_CELL_MATERIALIZATION":
        reasons.append("Source nodes: wrong qualification")
    _true(d, "source_only", reasons, label)
    _true(d, "state_order_matches_P2_nested_loops", reasons, label)
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "source_graph_rebuilt_or_pruned_here",
        "P4_metric_attached_here",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        _false(d, key, reasons, label)
    if int(d.get("partition", {}).get("states", 0) or 0) != 800:
        reasons.append("Source nodes: physical state partition is not 800 nodes")


def _check_clock(d: dict, reasons: list[str]) -> None:
    label = "Sample clock"
    _validated(d, reasons, label)
    if d.get("qualification") != "OU3_P2_SAMPLE_CLOCK_COMMIT_REACHABILITY_REFINEMENT":
        reasons.append("Sample clock: wrong qualification")
    if d.get("P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE") != "PASS":
        reasons.append("Sample clock: finite-speed source refinement did not pass")
    for key in (
        "source_only",
        "same_physical_partition_as_P2",
        "EMA_updated_every_valid_sample",
        "EMA_composed_sample_by_sample",
        "sample_varying_target_and_horizon_boxes_admitted",
        "commit_only_stages_current_smoothed_candidate",
        "pending_candidate_applied_before_next_sample",
        "arbitrary_late_commit_jump_removed",
        "frozen_clock_self_loop_included",
    ):
        _true(d, key, reasons, label)
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "source_graph_all_to_all",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        _false(d, key, reasons, label)
    if int(d.get("partition", {}).get("states", 0) or 0) != 800:
        reasons.append("Sample clock: physical state partition is not 800 nodes")
    clock = d.get("clock", {})
    if not isinstance(clock, dict) or clock.get("floating_clock_stagnation_verified") is not True:
        reasons.append("Sample clock: absorbing frozen-clock semantics are not certified")
    lo = d.get("finite_stage_gap_lower_samples")
    hi = d.get("finite_stage_gap_upper_samples")
    if not (isinstance(lo, int) and isinstance(hi, int) and 0 < lo <= hi):
        reasons.append("Sample clock: finite stage-gap range is invalid")


def _check_candidate(
    cand: dict | None,
    p3: dict,
    cayley: dict,
    remainder: dict,
    timing: dict,
    clock: dict,
    reasons: list[str],
) -> dict:
    if cand is None:
        reasons.append("complete H=18/A=21 nonlinear whole-word dissipation candidate is missing")
        return {
            "rho_upper": {"H": None, "A": None},
            "one_minus_rho_lower": {"H": None, "A": None},
        }

    label = "Complete-word candidate"
    if cand.get("qualification") != CANDIDATE_QUALIFICATION:
        reasons.append(f"{label}: wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "canonical_P3_artifact_consumed",
        "same_source_history_for_metric_and_nonlinear_word",
        "implemented_prediction_measurement_Joseph_order_covered",
        "source_complete_vector_packet_language_covered",
        "S_linear_timing_discharged_by_canonical_P3",
        "Cayley_exact_geometry_consumed",
        "homogeneous_vector_remainder_consumed",
        "finite_speed_sample_clock_graph_consumed",
        "frozen_clock_absorbing_hold_branch_covered",
        "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch",
        "full_H18_state_word_covered",
        "full_A21_state_word_covered",
        "signed_nonlinear_remainder_charged",
        "complete_word_generalized_Jacobian_or_equivalent_bound",
    ):
        _true(cand, key, reasons, label)
    for key in ("trajectory_replay_used", "filter_changed"):
        _false(cand, key, reasons, label)

    if cand.get("H_dimension") != H_DIM:
        reasons.append(f"{label}: H dimension is not 18")
    if cand.get("A_dimension") != A_DIM:
        reasons.append(f"{label}: A dimension is not 21")

    theta = cand.get("outer_angle_rad")
    ctheta = cayley.get("outer_angle_rad")
    rtheta = remainder.get("outer_angle_rad")
    if not all(_finite_number(x) for x in (theta, ctheta, rtheta)):
        reasons.append(f"{label}: outer-angle binding is missing")
    elif not (
        float(theta) >= REQUIRED_OUTER_ANGLE_RAD
        and float(theta) == float(ctheta)
        and float(theta) == float(rtheta)
    ):
        reasons.append(f"{label}: candidate is not bound to the exact 0.8-rad prerequisite sector")

    p3_worst = p3.get("worst_H_A_margin")
    consumed = cand.get("canonical_P3_worst_H_A_margin_consumed")
    if not (_finite_number(p3_worst) and _finite_number(consumed) and float(consumed) == float(p3_worst)):
        reasons.append(f"{label}: candidate is not bound to the exact canonical P3 metric margin")

    if cand.get("source_word_horizon_s") != timing.get("word_horizon_s"):
        reasons.append(f"{label}: source word horizon differs from the certified timing contract")
    if cand.get("sample_clock_transition_edges") != clock.get("transition_edges"):
        reasons.append(f"{label}: sample-clock edge family differs from the certified source graph")

    rhos: dict[str, float | None] = {}
    gaps: dict[str, float | None] = {}
    for mode in ("H", "A"):
        key = f"rho_{mode}_upper"
        gap_key = f"one_minus_rho_{mode}_lower"
        rho = cand.get(key)
        gap = cand.get(gap_key)
        rhos[mode] = float(rho) if _finite_number(rho) else None
        gaps[mode] = None
        if gap is None:
            # No exact strict gap supplied: strictness must be visible in the
            # reported float itself, exactly as for a representable rho.
            if rhos[mode] is None or not 0.0 <= rhos[mode] < 1.0:
                reasons.append(f"{label}: {key} is missing or is not strictly below 1")
            else:
                gaps[mode] = 1.0 - rhos[mode]
        else:
            if rhos[mode] is None or not 0.0 <= rhos[mode] <= 1.0:
                reasons.append(f"{label}: {key} is missing or is outside [0,1]")
            if not _finite_number(gap) or not 0.0 < float(gap) <= 1.0:
                reasons.append(f"{label}: {gap_key} is not a finite positive strict gap")
            else:
                gaps[mode] = float(gap)
                if rhos[mode] is not None and rhos[mode] > math.nextafter(1.0 - float(gap), math.inf):
                    reasons.append(f"{label}: {key} is inconsistent with the exact strict gap")
        margin_key = f"strict_dissipation_margin_{mode}_lower"
        margin = cand.get(margin_key)
        if not _finite_number(margin) or float(margin) <= 0.0:
            reasons.append(f"{label}: {margin_key} is not finite positive")
        elif gaps[mode] is not None and float(margin) > gaps[mode]:
            reasons.append(f"{label}: {margin_key} exceeds 1-rho_{mode}")

    if cand.get("P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED") is not True:
        reasons.append(f"{label}: producer does not explicitly establish complete-word dissipation")
    return {"rho_upper": rhos, "one_minus_rho_lower": gaps}


def build(
    p3: dict,
    cayley: dict,
    remainder: dict,
    timing: dict,
    path: dict,
    nodes: dict,
    clock: dict,
    candidate: dict | None = None,
) -> dict:
    reasons: list[str] = []
    _check_p3(p3, reasons)
    _check_cayley(cayley, reasons)
    _check_remainder(remainder, reasons)
    _check_timing(timing, reasons)
    _check_path(path, reasons)
    _check_nodes(nodes, reasons)
    _check_clock(clock, reasons)
    contraction = _check_candidate(candidate, p3, cayley, remainder, timing, clock, reasons)

    passed = not reasons
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_definition_frozen": True,
        "only_this_module_may_promote_P4_for_P5": True,
        "P3_useful_gate": USEFUL_P3_GATE,
        "required_outer_angle_rad": REQUIRED_OUTER_ANGLE_RAD,
        "required_dimensions": {"H": H_DIM, "A": A_DIM},
        "P3_canonical_qualification": p3.get("qualification"),
        "P3_canonical_pass_consumed": p3.get("P3_CANONICAL_PASS") is True,
        "P3_worst_H_A_margin": p3.get("worst_H_A_margin"),
        "source_partition_states": nodes.get("partition", {}).get("states"),
        "sample_clock_transition_edges": clock.get("transition_edges"),
        "source_word_horizon_s": timing.get("word_horizon_s"),
        "candidate_qualification": None if candidate is None else candidate.get("qualification"),
        "rho_upper": contraction["rho_upper"],
        "one_minus_rho_lower": contraction["one_minus_rho_lower"],
        "strict_gap_may_be_below_binary64_epsilon": True,
        "P4_CANONICAL_PASS": passed,
        "P4_CANONICAL_FAIL_REASONS": reasons,
        "P5_MAY_CONSUME_P4": passed,
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "freeze P4 and prove finite startup/outer-sector capture into its robust inner funnel"
            if passed
            else "construct or tighten the source-complete H=18/A=21 nonlinear whole-word dissipation candidate under this frozen P4 interface; do not promote sector/source prerequisites"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    if d.get("canonical_definition_frozen") is not True:
        f.append("canonical P4 definition is not frozen")
    if d.get("only_this_module_may_promote_P4_for_P5") is not True:
        f.append("canonical P4 promotion authority is not unique")
    if d.get("P3_useful_gate") != USEFUL_P3_GATE:
        f.append("P3 useful gate changed at the P4 interface")
    if d.get("required_outer_angle_rad") != REQUIRED_OUTER_ANGLE_RAD:
        f.append("required P4 outer angle changed")
    if d.get("required_dimensions") != {"H": H_DIM, "A": A_DIM}:
        f.append("required H/A dimensions changed")
    reasons = d.get("P4_CANONICAL_FAIL_REASONS")
    if not isinstance(reasons, list):
        f.append("P4 fail reasons are malformed")
        reasons = []
    expected = len(reasons) == 0
    if d.get("P4_CANONICAL_PASS") is not expected:
        f.append("P4 canonical pass flag does not match obligations")
    if d.get("P5_MAY_CONSUME_P4") is not d.get("P4_CANONICAL_PASS"):
        f.append("P5 consumption gate differs from canonical P4 result")
    if d.get("P5_FINITE_CAPTURE_ESTABLISHED_HERE") is not False:
        f.append("P4 gate prematurely established P5 finite capture")
    if d.get("strict_gap_may_be_below_binary64_epsilon") is not True:
        f.append("P4 strict-gap representation contract changed")
    gaps = d.get("one_minus_rho_lower")
    if not isinstance(gaps, dict) or set(gaps) != {"H", "A"}:
        f.append("P4 strict gaps are malformed")
    elif expected and any(not _finite_number(x) or float(x) <= 0.0 for x in gaps.values()):
        f.append("P4 passed without a finite positive strict gap in both modes")
    return list(dict.fromkeys(f))


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p3", type=Path, required=True)
    ap.add_argument("--cayley", type=Path, required=True)
    ap.add_argument("--remainder", type=Path, required=True)
    ap.add_argument("--timing", type=Path, required=True)
    ap.add_argument("--path", type=Path, required=True)
    ap.add_argument("--nodes", type=Path, required=True)
    ap.add_argument("--clock", type=Path, required=True)
    ap.add_argument("--candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()

    d = build(
        _load(a.p3),
        _load(a.cayley),
        _load(a.remainder),
        _load(a.timing),
        _load(a.path),
        _load(a.nodes),
        _load(a.clock),
        None if a.candidate is None else _load(a.candidate),
    )
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P3_canonical_pass_consumed": d["P3_canonical_pass_consumed"],
        "P3_worst_H_A_margin": d["P3_worst_H_A_margin"],
        "candidate_qualification": d["candidate_qualification"],
        "rho_upper": d["rho_upper"],
        "one_minus_rho_lower": d["one_minus_rho_lower"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "P4_CANONICAL_FAIL_REASONS": d["P4_CANONICAL_FAIL_REASONS"],
        "P5_MAY_CONSUME_P4": d["P5_MAY_CONSUME_P4"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
