#!/usr/bin/env python3
"""Structural P3 probe for a horizon-matched whole-word translation floor.

PR #475 corrected the P3 diagnosis: canonical ``delta`` compares a zero-start
13--26-sample translation floor against a 3.02--3.17 s whole-word covariance
ceiling.  Tightening the ceiling alone is not enough.  This producer asks the
single falsifiable question recorded by the research ledger: can a rigorous
translation lower be retained across a complete covariance word without
changing the filter, domain, source contract, or the frozen 1e-18 gate?

The probe deliberately starts with a more conservative source projection than
an exact P2 path frontier.  For each of the ten tau cells, node ``tau*80`` is an
*actual* P2 physical cell with that same tau interval and with the smallest
filter-side sigma and R_S intervals in the tau cell.  Riccati covariance is
Loewner monotone increasing in process covariance Q and in measurement R, so
replacing any physical node by this same-tau node can only lower covariance.
We then allow any tau cell and any certified 13..26-sample gap after every
segment.  That only ADDS source paths.  Consequently a lower that survives this
projection is valid for every legal P2 V1 history, while never constructing an
independent Cartesian tau/sigma/R_S tuple.

To keep the path merge aligned with the actual theorem quantity, every complete
segment is projected to one scalar in the *fixed canonical translation metric*

    D_h P_z D_h' >= delta * diag(Sigma_upper_global),

where ``Sigma_upper_global`` is formed only AFTER each legal P2 history label is
mapped to its same-history covariance upper.  The corresponding diagonal
``P_z`` lower is fed into the next segment.  This loses cross-correlation but
preserves the v/p/S/a_w anisotropy that the old ``rho I`` merge discarded.
No interval covariance boxes are recursively hulled across source changes.
Within one tau cell the retained, already-validated segment algebra is reused.

A 635-sample word can lose at most one <26-sample partial source segment at each
end.  Therefore it contains at least 23 complete 13..26-sample segments.  The
scalar transfer is monotone in its incoming lower; starting at zero makes the
certified sequence nondecreasing.  After 23 complete segments the probe then
checks all 0..25 end phases under the same source over-approximation.

The strongest accelerometer and S measurements are still applied every sample,
as in the existing P3 lower.  This is intentionally stronger than the deployed
pseudo cadence and therefore remains a valid covariance lower.  If even this
conservative whole-word floor clears the unchanged gate against the *existing*
whole-word upper, cadence-crediting the ceiling is no longer needed to establish
that the horizon mismatch was the limiter (it remains a later tightening inside
the same research question).  If it does not clear, this projection is recorded
as falsified; do not tune x subdivisions or theorem thresholds around it.

This is a non-promoting research producer.  It cannot set P3/P4/P5 PASS.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_segment as SEG
import ou3_p3_frozen_full_matrix_translation as FROZEN
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_p3_p2_v1_stage_phase_translation as STAGE
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
TAU_CELLS = 10
SIGMA_CELLS = 8
RS_CELLS = 10
TAU_STRIDE = SIGMA_CELLS * RS_CELLS
PHASE_MAX = 25


def _copy_matrix(A):
    return [[A[i][j] for j in range(len(A[i]))] for i in range(len(A))]


def _same_interval_matrix(A, B) -> bool:
    if len(A) != len(B):
        return False
    for i in range(len(A)):
        if len(A[i]) != len(B[i]):
            return False
        for j in range(len(A[i])):
            if float(A[i][j].lo) != float(B[i][j].lo):
                return False
            if float(A[i][j].hi) != float(B[i][j].hi):
                return False
    return True


def guaranteed_complete_segments(target_samples: int,
                                 gap_max: int = 26,
                                 phase_max: int = PHASE_MAX) -> int:
    """Lower-bound complete source segments in an arbitrary target window.

    At most ``phase_max`` samples can be lost before the first complete segment
    and at most ``phase_max`` after the last one.  Every complete segment has at
    most ``gap_max`` samples.  The ceiling division is therefore the minimum
    number needed to cover the remaining samples.
    """
    N = int(target_samples)
    g = int(gap_max)
    p = int(phase_max)
    if N <= 2 * p or g <= 0 or p < 0:
        raise ValueError("invalid whole-word segment-count inputs")
    return int(math.ceil((N - 2 * p) / g))


def _dominator_nodes(rt: dict) -> list[dict]:
    """Return and verify the ten same-tau minimum-(sigma,R_S) physical nodes."""
    nodes = rt["nodes"]
    if len(nodes) != TAU_CELLS * TAU_STRIDE:
        raise RuntimeError("P2 physical source partition is no longer 10x8x10")
    out = []
    for ti in range(TAU_CELLS):
        idx = ti * TAU_STRIDE
        d = nodes[idx]
        if int(d.get("tau_index", -1)) != ti:
            raise RuntimeError("same-tau dominator index lost tau ordering")
        if int(d.get("sigma_raw_index", -1)) != 0 or int(d.get("R_S_index", -1)) != 0:
            raise RuntimeError("same-tau dominator is not the physical (sigma0,R_S0) node")
        dtau = tuple(map(float, d["tau_s"]))
        dsig = tuple(map(float, d["sigma_filter_committed_mps2"]))
        drs = tuple(map(float, d["R_S_filter_std"]))
        for n in nodes[idx:idx + TAU_STRIDE]:
            if tuple(map(float, n["tau_s"])) != dtau:
                raise RuntimeError("tau block no longer has a common tau cell")
            nsig = tuple(map(float, n["sigma_filter_committed_mps2"]))
            nrs = tuple(map(float, n["R_S_filter_std"]))
            if dsig[0] > nsig[0] or dsig[1] > nsig[1]:
                raise RuntimeError("chosen same-tau sigma node does not Loewner-dominate process lower")
            if drs[0] > nrs[0] or drs[1] > nrs[1]:
                raise RuntimeError("chosen same-tau R_S node is not the strongest measurement cell")
        out.append(d)
    return out


def _global_same_history_upper(fr: dict) -> tuple[list[float], dict]:
    """Envelope legal-history covariance results only after same-history mapping."""
    env = [0.0] * 4
    rows = 0
    phase0_rows = 0
    positive_phase_rows = 0
    for t in fr["endpoint_nodes"]:
        for phase in (0, 1):
            row = HIST.endpoint_phase_upper(int(t), phase, fr)
            if row.get("same_history_upper_evaluated_before_endpoint_envelope") is not True:
                raise RuntimeError("history upper lost same-history-before-envelope ordering")
            if row.get("raw_tuner_cartesian_extrema_used") is not False:
                raise RuntimeError("history upper returned to Cartesian tuner extrema")
            u = list(map(float, row["Sigma_translation_diagonal_upper_envelope"]))
            if len(u) != 4 or any(not (math.isfinite(x) and x > 0.0) for x in u):
                raise RuntimeError("invalid same-history translation upper")
            env = [max(env[i], u[i]) for i in range(4)]
            rows += 1
            if phase == 0:
                phase0_rows += 1
            else:
                positive_phase_rows += 1
    if any(not (math.isfinite(x) and x > 0.0) for x in env):
        raise RuntimeError("global same-history covariance envelope is empty")
    return env, {
        "endpoint_phase_classes_enveloped": rows,
        "phase0_endpoint_classes": phase0_rows,
        "positive_phase_endpoint_classes": positive_phase_rows,
        "same_history_upper_evaluated_before_global_envelope": True,
        "raw_tuner_cartesian_extrema_used": False,
    }


def _metric_lower(delta: float, h: float, upper: list[float]):
    """Map a certified physical metric margin back to a diagonal z-state lower."""
    x = float(delta)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("nonnegative finite metric margin required")
    if x == 0.0:
        return FROZEN._mat_zero()
    d = (h, h * h, h * h * h, 1.0)
    out = [[Interval.point(0.0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        d2_hi = BASE.up(d[i] * d[i])
        num_lo = BASE.down(x * float(upper[i]))
        v = BASE.down(num_lo / d2_hi)
        if not (math.isfinite(v) and v > 0.0):
            raise RuntimeError("metric lower underflowed or lost positivity")
        out[i][i] = Interval.point(v)
    return out


def _physical_diagonal_ratios(Pz, h: float, upper: list[float]) -> list[float]:
    d = (h, h * h, h * h * h, 1.0)
    ans = []
    for i in range(4):
        physical = BASE.down(BASE.down(d[i] * d[i]) * float(Pz[i][i].lo))
        ans.append(BASE.down(max(0.0, physical) / float(upper[i])))
    return ans


def _prepare_kernel(node: dict, x: Interval, h: float) -> dict:
    """Precompute only source/x terms that are invariant across covariance inputs."""
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    sigma2_lower = BASE.down(sigma.lo * sigma.lo)
    if not sigma2_lower > 0.0:
        raise RuntimeError("same-tau dominator sigma lower lost positivity")
    Fm = FROZEN._transition(x)
    Ft = FROZEN._transpose(Fm)
    Qz = SEG._scale(FROZEN._scaled_Q(x), SEG._point(sigma2_lower))
    R_aw, R_S_z = SEG._physical_measurement_variances(node, h)
    return {"F": Fm, "Ft": Ft, "Qz": Qz, "R_aw": float(R_aw), "R_S_z": float(R_S_z)}


def _prepared_one_step(P, kernel: dict):
    """Exact SEG.one_step algebra with immutable source/x objects cached."""
    L0, _ = SEG._common_point_lower(P)
    prediction_interval = FROZEN._sym(
        FROZEN._add(
            FROZEN._mul(FROZEN._mul(kernel["F"], L0), kernel["Ft"]),
            kernel["Qz"],
        )
    )
    pred, _ = SEG._common_point_lower(prediction_interval)
    post_aw, _ = SEG._point_measurement_lower(pred, 3, kernel["R_aw"])
    post_s, _ = SEG._point_measurement_lower(post_aw, 2, kernel["R_S_z"])
    return post_s


def _prepare_leaf(node: dict, x: Interval, h: float, depth: int = 0) -> list[dict]:
    try:
        kernel = _prepare_kernel(node, x, h)
    except (RuntimeError, ValueError) as exc:
        if not SEG._split_request(exc) or depth >= SEG.MAX_ADAPTIVE_X_DEPTH:
            raise
        halves = SEG._split_x(x)
        if halves is None:
            raise RuntimeError(f"cannot split whole-word x cell {x.as_list()}") from exc
        return _prepare_leaf(node, halves[0], h, depth + 1) + _prepare_leaf(node, halves[1], h, depth + 1)

    zero = FROZEN._mat_zero()
    direct = SEG.one_step(zero, node, x, h)
    prepared = _prepared_one_step(zero, kernel)
    if not _same_interval_matrix(direct, prepared):
        raise RuntimeError("prepared kernel changed retained SEG.one_step interval endpoints")
    return [{"x": x, "kernel": kernel}]


def _prepare_sources(rt: dict, dominators: list[dict]) -> tuple[list[dict], int]:
    h = float(rt["clock"]["dt_binary32_s"])
    sources = []
    checked = 0
    for node in dominators:
        leaves = []
        for x in SEG._node_x_subcells(node, h, SEG.X_SUBCELLS):
            leaves.extend(_prepare_leaf(node, x, h, 0))
        checked += len(leaves)
        sources.append({"node": node, "leaves": leaves})
    return sources, checked


def _candidate_delta(Pz, h: float, upper: list[float], current_best: float | None) -> float | None:
    """Return candidate delta only when it can lower the running minimum.

    Once one candidate has supplied ``current_best``, a single validated LDL
    test of ``P-current_best*U`` is enough to discard every non-limiting
    candidate.  Only candidates that fail that test pay the full 52-step
    certified bisection used by canonical P3.
    """
    if current_best is not None and current_best > 0.0:
        if symmetric_positive_definite_ldlt(
            STAGE._physical_gate_matrix(Pz, h, upper, current_best)
        )[0]:
            return None
    return STAGE._certified_delta(Pz, h, upper)


def _scan_complete_segment(P0, sources: list[dict], gaps: tuple[int, ...],
                           h: float, upper: list[float]) -> tuple[float, dict]:
    best = None
    limiting = None
    coordinates = ("v", "p", "S", "a_w")
    # Long tau and short gaps have historically been the difficult directions;
    # scan them first so the running-minimum shortcut becomes effective early.
    ordered_sources = list(reversed(sources))
    ordered_gaps = tuple(sorted(gaps))
    gap_set = set(ordered_gaps)
    max_gap = max(ordered_gaps)

    for src in ordered_sources:
        node = src["node"]
        for leaf in src["leaves"]:
            P = _copy_matrix(P0)
            for sample in range(1, max_gap + 1):
                P = _prepared_one_step(P, leaf["kernel"])
                if sample not in gap_set:
                    continue
                cand = _candidate_delta(P, h, upper, best)
                if cand is None:
                    continue
                if best is None or cand < best:
                    ratios = _physical_diagonal_ratios(P, h, upper)
                    ci = min(range(4), key=lambda i: ratios[i])
                    best = float(cand)
                    limiting = {
                        "tau_index": int(node["tau_index"]),
                        "physical_source_node": int(node["index"]),
                        "gap_samples": int(sample),
                        "x_h_over_tau": leaf["x"].as_list(),
                        "translation_relative_Riccati_injection_margin_lower": float(cand),
                        "minimum_physical_diagonal_ratio_lower": float(ratios[ci]),
                        "minimum_physical_diagonal_ratio_coordinate": coordinates[ci],
                        "physical_diagonal_ratio_lower": ratios,
                    }
    if best is None or limiting is None or not (math.isfinite(best) and best > 0.0):
        raise RuntimeError("whole-word complete-segment transfer lost a strict metric floor")
    return BASE.down(best), limiting


def _scan_end_phase(P0, sources: list[dict], h: float,
                    upper: list[float]) -> tuple[float, dict]:
    best = STAGE._certified_delta(P0, h, upper)
    if not (math.isfinite(best) and best > 0.0):
        raise RuntimeError("whole-word boundary metric lower is not strict")
    limiting = {
        "phase_samples": 0,
        "translation_relative_Riccati_injection_margin_lower": float(best),
        "kind": "stage_boundary",
    }
    coordinates = ("v", "p", "S", "a_w")
    for src in reversed(sources):
        node = src["node"]
        for leaf in src["leaves"]:
            P = _copy_matrix(P0)
            for sample in range(1, PHASE_MAX + 1):
                P = _prepared_one_step(P, leaf["kernel"])
                cand = _candidate_delta(P, h, upper, best)
                if cand is None:
                    continue
                if cand < best:
                    ratios = _physical_diagonal_ratios(P, h, upper)
                    ci = min(range(4), key=lambda i: ratios[i])
                    best = float(cand)
                    limiting = {
                        "kind": "positive_end_phase",
                        "phase_samples": int(sample),
                        "tau_index": int(node["tau_index"]),
                        "physical_source_node": int(node["index"]),
                        "x_h_over_tau": leaf["x"].as_list(),
                        "translation_relative_Riccati_injection_margin_lower": float(cand),
                        "minimum_physical_diagonal_ratio_lower": float(ratios[ci]),
                        "minimum_physical_diagonal_ratio_coordinate": coordinates[ci],
                        "physical_diagonal_ratio_lower": ratios,
                    }
    return BASE.down(best), limiting


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("whole-word P3 probe must not be trajectory fitted")
    runtime_cfg = domain.get("configured_runtime", {})
    if runtime_cfg.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("whole-word P3 probe requires declared zero lever arm")
    if runtime_cfg.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        raise RuntimeError("whole-word P3 probe requires dormant transparent vibration guard")
    if float(BASE.MIN_USEFUL_DELTA) != 1.0e-18:
        raise RuntimeError("canonical P3 useful gate changed")

    fr = HIST.frontier_runtime(path)
    rt = fr["rt"]
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 correlation interface failed: {cf}")
    if corr.get("interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("P2 correlation interface version changed")

    upper, upper_meta = _global_same_history_upper(fr)
    dominators = _dominator_nodes(rt)
    sources, equivalence_checks = _prepare_sources(rt, dominators)
    h = float(rt["clock"]["dt_binary32_s"])
    gaps = tuple(map(int, rt["gaps"]))
    if gaps != tuple(range(13, 27)):
        raise RuntimeError("P2 finite source gap alphabet changed")

    target_samples = int(fr["target"]["target_samples"])
    complete = guaranteed_complete_segments(target_samples, max(gaps), PHASE_MAX)
    if complete < 1:
        raise RuntimeError("whole-word target contains no guaranteed complete source segment")

    P = FROZEN._mat_zero()
    margins = []
    previous = 0.0
    monotone = True
    for k in range(1, complete + 1):
        delta, limiting = _scan_complete_segment(P, sources, gaps, h, upper)
        if delta < previous:
            monotone = False
            raise RuntimeError(
                f"certified whole-word metric recurrence decreased at segment {k}: {delta} < {previous}"
            )
        margins.append({
            "complete_segments": k,
            "translation_relative_Riccati_injection_margin_lower": delta,
            "limiting_transfer": limiting,
        })
        P = _metric_lower(delta, h, upper)
        previous = delta

    boundary_delta = previous
    final_delta, final_limiting = _scan_end_phase(P, sources, h, upper)
    useful = final_delta >= BASE.MIN_USEFUL_DELTA
    classification = (
        "WHOLE_WORD_HORIZON_MATCHED_FLOOR_CLEARS_CANONICAL_GATE"
        if useful else
        "CONSERVATIVE_WHOLE_WORD_SOURCE_PROJECTION_STILL_BELOW_CANONICAL_GATE"
    )
    next_obligation = (
        "replace the one-segment canonical P3 translation floor by the horizon-matched construction, "
        "then add the separately rigorous cell-local S-cadence ceiling tightening before H/A promotion"
        if useful else
        "perform the required critic pass on the arbitrary-switch/same-tau-dominator projection; "
        "do not tune x subdivisions or the 1e-18 gate"
    )

    physical_metric_lower = [BASE.down(final_delta * x) for x in upper]
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_WHOLE_WORD_HORIZON_MATCHED_TRANSLATION_FLOOR_PROBE",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "canonical_gate_changed": False,
        "canonical_useful_gate": BASE.MIN_USEFUL_DELTA,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "same_history_upper_evaluated_before_global_envelope": True,
        "raw_tuner_cartesian_extrema_used": False,
        "lower_projection_kind": "same_tau_physical_min_sigma_min_R_S_node_then_arbitrary_tau_segment_switching",
        "lower_projection_uses_actual_physical_nodes": True,
        "lower_projection_can_only_add_source_paths": True,
        "lower_projection_can_only_reduce_covariance": True,
        "independent_tau_sigma_R_S_coordinate_extremization_used": False,
        "same_tau_dominator_nodes": [int(x["index"]) for x in dominators],
        "same_tau_dominator_count": len(dominators),
        "prepared_kernel_exact_SEG_one_step_equivalence_subcells": equivalence_checks,
        "strongest_accelerometer_measurement_every_sample_retained": True,
        "strongest_S_measurement_every_sample_retained": True,
        "deployed_S_cadence_used_in_lower": False,
        "lower_is_conservative_against_deployed_S_cadence": True,
        "global_same_history_Sigma_translation_diagonal_upper_envelope": upper,
        "upper_envelope_metadata": upper_meta,
        "whole_word_target": dict(fr["target"]),
        "phase_prefix_suffix_max_samples": PHASE_MAX,
        "guaranteed_complete_segments": complete,
        "complete_segment_gap_alphabet_samples": list(gaps),
        "whole_word_metric_recurrence_monotone": monotone,
        "complete_segment_margin_sequence": margins,
        "whole_word_boundary_margin_lower": boundary_delta,
        "whole_word_arbitrary_phase_margin_lower": final_delta,
        "whole_word_arbitrary_phase_limiting_transfer": final_limiting,
        "whole_word_physical_metric_diagonal_lower": physical_metric_lower,
        "structural_probe_useful_at_frozen_gate": useful,
        "classification": classification,
        "source_complete_whole_word_floor_established_by_probe": True,
        "canonical_P3_translation_route_replaced": False,
        "cadence_aware_whole_word_ceiling_established": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "next_obligation": next_obligation,
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_WHOLE_WORD_HORIZON_MATCHED_TRANSLATION_FLOOR_PROBE":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "same_history_upper_evaluated_before_global_envelope",
        "lower_projection_uses_actual_physical_nodes",
        "lower_projection_can_only_add_source_paths",
        "lower_projection_can_only_reduce_covariance",
        "strongest_accelerometer_measurement_every_sample_retained",
        "strongest_S_measurement_every_sample_retained",
        "lower_is_conservative_against_deployed_S_cadence",
        "whole_word_metric_recurrence_monotone",
        "source_complete_whole_word_floor_established_by_probe",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "canonical_gate_changed", "raw_tuner_cartesian_extrema_used",
        "independent_tau_sigma_R_S_coordinate_extremization_used",
        "deployed_S_cadence_used_in_lower", "canonical_P3_translation_route_replaced",
        "cadence_aware_whole_word_ceiling_established", "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("canonical_useful_gate", math.nan)) != 1.0e-18:
        f.append("canonical useful gate changed")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("whole-word probe lost frozen P2 V1 binding")
    if d.get("same_tau_dominator_nodes") != [i * TAU_STRIDE for i in range(TAU_CELLS)]:
        f.append("same-tau physical dominator node list changed")
    if int(d.get("same_tau_dominator_count", 0)) != TAU_CELLS:
        f.append("whole-word lower does not cover all ten tau cells")
    target = d.get("whole_word_target", {})
    N = int(target.get("target_samples", 0) or 0)
    expected = guaranteed_complete_segments(N, 26, PHASE_MAX) if N > 2 * PHASE_MAX else 0
    if int(d.get("guaranteed_complete_segments", -1)) != expected or expected <= 0:
        f.append("whole-word guaranteed complete-segment count is wrong")
    seq = list(d.get("complete_segment_margin_sequence", []))
    if len(seq) != expected:
        f.append("whole-word recurrence did not cover every guaranteed complete segment")
    prev = 0.0
    for i, row in enumerate(seq, 1):
        if int(row.get("complete_segments", -1)) != i:
            f.append("whole-word recurrence segment numbering changed")
            break
        x = float(row.get("translation_relative_Riccati_injection_margin_lower", math.nan))
        if not (math.isfinite(x) and x > 0.0 and x >= prev):
            f.append("whole-word recurrence margin is not positive nondecreasing")
            break
        prev = x
    final = float(d.get("whole_word_arbitrary_phase_margin_lower", math.nan))
    if not (math.isfinite(final) and final > 0.0):
        f.append("whole-word arbitrary-phase margin is not strict")
    upper = list(d.get("global_same_history_Sigma_translation_diagonal_upper_envelope", []))
    if len(upper) != 4 or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in upper):
        f.append("invalid legal-history covariance envelope")
    useful = final >= BASE.MIN_USEFUL_DELTA if math.isfinite(final) else False
    if d.get("structural_probe_useful_at_frozen_gate") is not useful:
        f.append("whole-word structural verdict does not match frozen gate")
    if not isinstance(d.get("classification"), str) or not d["classification"]:
        f.append("missing whole-word structural classification")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "whole_word_target": d["whole_word_target"],
        "same_tau_dominator_nodes": d["same_tau_dominator_nodes"],
        "guaranteed_complete_segments": d["guaranteed_complete_segments"],
        "global_upper": d["global_same_history_Sigma_translation_diagonal_upper_envelope"],
        "margin_sequence": [
            [r["complete_segments"], r["translation_relative_Riccati_injection_margin_lower"]]
            for r in d["complete_segment_margin_sequence"]
        ],
        "boundary_margin_lower": d["whole_word_boundary_margin_lower"],
        "arbitrary_phase_margin_lower": d["whole_word_arbitrary_phase_margin_lower"],
        "limiting_transfer": d["whole_word_arbitrary_phase_limiting_transfer"],
        "classification": d["classification"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
