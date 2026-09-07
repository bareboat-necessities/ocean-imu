#!/usr/bin/env python3
"""Correlated finite-window source cells for complete SEA3 P4 covering.

This module supplies the *cover data model* used after the connected selector
execution exists.  It does not manufacture the still-open SEA0 hard realization
oracle.  A ``Sea3WindowCell`` is an intersection cell inside the one canonical
complete-SEA3 window set: every child retains the parent's common hard
realization, lambda-coupling, joint-response and source-history constraints.

Refinement is therefore allowed only on genuine source/generator coordinates.
Derived shipping quantities (tau, sigma_aw, R_S, T_S, P, H, K, F, Q, gains,
front-end outputs, etc.) are never split independently.  They must be recomputed
by the trusted same-history executor from each refined source cell.

A binary split means

    C = (C intersect {g <= m}) union (C intersect {g >= m}),

with *all* non-Cartesian SEA3 constraints retained on both children.  The bound
stored for ``g`` is only a search hull; membership remains constrained by the
shared SEA3 witness and coupled constraints.  Thus splitting one H_r search
hull does not discard H_s^2=sum H_r^2 or the active-partition steepness
constraint.

The second half of the module attaches one such window cell to every prefix of
one #500 selector lineage.  A21 absolute physical-bias cells can be attached to
the same exact prefixes, but source-uniform bias closure remains false until a
hard pathwise bias coordinate is supplied by SEA0; covariance/stationary
variance is not accepted as such a bound.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from ou3_interval import Interval
import ou3_p4_complete_sea3_same_history_prefix_selectors as SELECTORS
import ou3_sea3_hard_finite_window_source as HARD
import ou3_sea3_rlambda_transition as RLAMBDA
import ou3_mems_bias_contract as BIAS

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_CORRELATED_FINITE_WINDOW_SOURCE_CELLS_V1"
CANONICAL_SOURCE = "COMPLETE_SEA3_NORMAL_LIVE_WORD"

_ALLOWED_SOURCE_PREFIXES = (
    "lambda.",
    "xs.",
    "hard_realization.",
    "response.",
    "physical_bias.",
)
_FORBIDDEN_DERIVED_PREFIXES = (
    "tau",
    "sigma",
    "r_s",
    "rs",
    "t_s",
    "ts",
    "p_",
    "riccati",
    "kalman",
    "gain",
    "f_",
    "q_",
    "frontend.",
    "tuner.",
    "scheduler.",
)


@dataclass(frozen=True)
class SourceBound:
    coordinate: str
    interval: Interval
    provenance: str

    def __post_init__(self) -> None:
        name = self.coordinate.strip().lower()
        if not any(name.startswith(prefix) for prefix in _ALLOWED_SOURCE_PREFIXES):
            raise ValueError(f"not a refinable complete-SEA3 source coordinate: {self.coordinate}")
        if any(name.startswith(prefix) for prefix in _FORBIDDEN_DERIVED_PREFIXES):
            raise ValueError(f"derived shipping coordinate cannot define source refinement: {self.coordinate}")
        if not isinstance(self.interval, Interval):
            raise TypeError("source search hull must be an outward Interval")
        if not isinstance(self.provenance, str) or not self.provenance:
            raise ValueError("source bound requires provenance")


@dataclass(frozen=True)
class Sea3WindowCell:
    cell_id: str
    parent_cell_id: str | None
    depth: int
    source_identity: str
    hard_realization_witness_id: str
    joint_response_witness_id: str
    source_bounds: tuple[SourceBound, ...]
    shared_constraint_tokens: tuple[str, ...]
    bias_path_witness_id: str | None = None
    window_horizon_s: float = 3.0
    sample_period_s: float = 0.005
    complete_window_samples: int = 601

    def bound(self, coordinate: str) -> SourceBound:
        matches = [b for b in self.source_bounds if b.coordinate == coordinate]
        if len(matches) != 1:
            raise KeyError(f"source coordinate {coordinate!r} is not unique in cell")
        return matches[0]


@dataclass(frozen=True)
class AttachedPrefix:
    selector_source_cell_id: str
    selector_parent_source_cell_id: str
    prefix_length: int
    window_source_cell_id: str
    canonical_source_identity: str
    absolute_bias_true: tuple[Interval, Interval, Interval] | None
    absolute_bias_same_history: bool


def _validate_cell(cell: Sea3WindowCell) -> list[str]:
    f: list[str] = []
    if not cell.cell_id:
        f.append("empty window source cell id")
    if cell.depth < 0:
        f.append("negative source-cell depth")
    if cell.source_identity != CANONICAL_SOURCE:
        f.append("source cell detached from canonical complete SEA3")
    if not cell.hard_realization_witness_id:
        f.append("source cell lost hard-realization witness")
    if not cell.joint_response_witness_id:
        f.append("source cell lost joint-response witness")
    if not cell.shared_constraint_tokens:
        f.append("source cell lost shared non-Cartesian constraints")
    if not math.isclose(cell.window_horizon_s, 3.0, rel_tol=0.0, abs_tol=1e-12):
        f.append("source cell changed 3 s horizon")
    if not math.isclose(cell.sample_period_s, 0.005, rel_tol=0.0, abs_tol=1e-12):
        f.append("source cell changed 5 ms sample period")
    if cell.complete_window_samples != 601:
        f.append("source cell does not cover all 601 source transitions")
    names = [b.coordinate for b in cell.source_bounds]
    if len(set(names)) != len(names):
        f.append("duplicate source search coordinate")
    return f


def split_source_cell(
    cell: Sea3WindowCell,
    coordinate: str,
    *,
    split_value: float | None = None,
) -> tuple[Sea3WindowCell, Sea3WindowCell]:
    """Bisect one source-search hull while retaining every coupled constraint."""
    failures = _validate_cell(cell)
    if failures:
        raise ValueError(f"cannot split invalid SEA3 source cell: {failures}")
    target = cell.bound(coordinate)
    lo = float(target.interval.lo)
    hi = float(target.interval.hi)
    if not lo < hi:
        raise ValueError("cannot split a point source coordinate")
    mid = 0.5 * (lo + hi) if split_value is None else float(split_value)
    if not (math.isfinite(mid) and lo < mid < hi):
        raise ValueError("source split must lie strictly inside current search hull")

    def child(side: str, interval: Interval) -> Sea3WindowCell:
        bounds = tuple(
            SourceBound(b.coordinate, interval, b.provenance)
            if b.coordinate == coordinate
            else b
            for b in cell.source_bounds
        )
        relation = f"split:{coordinate}:{side}:{mid:.17g}"
        return Sea3WindowCell(
            cell_id=f"{cell.cell_id}/{coordinate}:{side}",
            parent_cell_id=cell.cell_id,
            depth=cell.depth + 1,
            source_identity=cell.source_identity,
            hard_realization_witness_id=cell.hard_realization_witness_id,
            joint_response_witness_id=cell.joint_response_witness_id,
            bias_path_witness_id=cell.bias_path_witness_id,
            source_bounds=bounds,
            shared_constraint_tokens=cell.shared_constraint_tokens + (relation,),
            window_horizon_s=cell.window_horizon_s,
            sample_period_s=cell.sample_period_s,
            complete_window_samples=cell.complete_window_samples,
        )

    left = child("lo", Interval.outward_bounds(lo, mid))
    right = child("hi", Interval.outward_bounds(mid, hi))
    return left, right


def validate_binary_split(
    parent: Sea3WindowCell,
    children: Sequence[Sea3WindowCell],
    coordinate: str,
) -> list[str]:
    f: list[str] = []
    if len(children) != 2:
        return ["binary source refinement must have exactly two children"]
    try:
        pb = parent.bound(coordinate)
        cb = [child.bound(coordinate) for child in children]
    except KeyError as exc:
        return [str(exc)]
    for child in children:
        f.extend(_validate_cell(child))
        if child.parent_cell_id != parent.cell_id:
            f.append("source child parent id changed")
        if child.depth != parent.depth + 1:
            f.append("source child depth changed")
        if child.source_identity != parent.source_identity:
            f.append("source child identity changed")
        if child.hard_realization_witness_id != parent.hard_realization_witness_id:
            f.append("source child changed hard-realization witness")
        if child.joint_response_witness_id != parent.joint_response_witness_id:
            f.append("source child changed joint-response witness")
        if child.bias_path_witness_id != parent.bias_path_witness_id:
            f.append("source child changed bias-path witness")
        if child.shared_constraint_tokens[: len(parent.shared_constraint_tokens)] != parent.shared_constraint_tokens:
            f.append("source child discarded a coupled parent constraint")
        for b in parent.source_bounds:
            if b.coordinate == coordinate:
                continue
            if child.bound(b.coordinate) != b:
                f.append(f"source split changed unrelated coordinate {b.coordinate}")
    lo = min(cb[0].interval.lo, cb[1].interval.lo)
    hi = max(cb[0].interval.hi, cb[1].interval.hi)
    if lo != pb.interval.lo or hi != pb.interval.hi:
        f.append("source split children do not cover parent hull")
    intersection_lo = max(cb[0].interval.lo, cb[1].interval.lo)
    intersection_hi = min(cb[0].interval.hi, cb[1].interval.hi)
    if intersection_lo != intersection_hi:
        f.append("source split children overlap away from one shared boundary")
    return list(dict.fromkeys(f))


def attach_window_cell_to_endpoint_lineage(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
    window_cell: Sea3WindowCell,
    *,
    absolute_bias_true_by_selector_cell: Mapping[
        str, tuple[Interval, Interval, Interval]
    ] | None = None,
) -> tuple[AttachedPrefix, ...]:
    """Attach one correlated window cell and optional A21 bias path to every prefix."""
    failures = _validate_cell(window_cell)
    if failures:
        raise ValueError(f"invalid complete-SEA3 window source cell: {failures}")
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint_source_cell_id)
    ids = {s.source_cell_id for s in lineage}
    if absolute_bias_true_by_selector_cell is not None:
        if set(absolute_bias_true_by_selector_cell) != ids:
            raise RuntimeError("absolute physical-bias attachment must cover every and only selector prefix")
    attached: list[AttachedPrefix] = []
    for selector in lineage:
        bias = None
        same_history = False
        if absolute_bias_true_by_selector_cell is not None:
            bias = absolute_bias_true_by_selector_cell[selector.source_cell_id]
            if len(bias) != 3 or any(not isinstance(x, Interval) for x in bias):
                raise TypeError("absolute physical-bias source cell must contain three intervals")
            # This flag establishes ID attachment only. BIAS1 dynamics are
            # checked by the nonlinear lineage's retained root/GM graph.
            same_history = True
        attached.append(
            AttachedPrefix(
                selector_source_cell_id=selector.source_cell_id,
                selector_parent_source_cell_id=selector.parent_source_cell_id,
                prefix_length=selector.prefix_length,
                window_source_cell_id=window_cell.cell_id,
                canonical_source_identity=window_cell.source_identity,
                absolute_bias_true=bias,
                absolute_bias_same_history=same_history,
            )
        )
    return tuple(attached)


def _root_smoke_cell() -> Sea3WindowCell:
    constraints = (
        "H_s^2=sum_r H_r^2",
        "active partition peak-steepness",
        "lambda_{k+1} in coupled Rhat_lambda(lambda_k)",
        "one common continuum phase history with no reseed",
        "one common joint translational/rotational response witness",
        "B^601_SEA3 common-witness membership",
    )
    return Sea3WindowCell(
        cell_id="sea3-root",
        parent_cell_id=None,
        depth=0,
        source_identity=CANONICAL_SOURCE,
        hard_realization_witness_id="X^s_SEA3(lambda_0:600)",
        joint_response_witness_id="G_imu:common-window",
        bias_path_witness_id=None,
        source_bounds=(
            SourceBound("lambda.H1_fraction", Interval.outward_bounds(0.0, 1.0), "coupled Lambda_SEA3 search hull"),
            SourceBound("lambda.nu1", Interval.outward_bounds(0.0, 1.0), "compact R_lambda coordinate"),
            SourceBound("xs.support_coordinate", Interval.outward_bounds(-1.0, 1.0), "symbolic hard-realization search coordinate"),
            SourceBound("response.gain_fraction", Interval.outward_bounds(0.0, 1.0), "joint response-family search hull"),
        ),
        shared_constraint_tokens=constraints,
    )


def _smoke(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    hard = HARD.build(domain_path)
    rlambda = RLAMBDA.build(domain_path)
    root = _root_smoke_cell()
    left, right = split_source_cell(root, "lambda.H1_fraction")
    split_failures = validate_binary_split(root, (left, right), "lambda.H1_fraction")

    frontend = SELECTORS.FRONTEND._point_state()
    P0_H, P0_A, _ = SELECTORS._live_structured_point_covariance_fixture(frontend, domain_path)
    sample = SELECTORS._point_sample()
    endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
        frontend_entry=frontend,
        P0_H=P0_H,
        P0_A=P0_A,
        samples=[sample, sample],
        domain_path=domain_path,
        branch_limit=128,
    )
    endpoint = endpoints[0].source_cell_id
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint)
    z = (Interval.point(0.0), Interval.point(0.0), Interval.point(0.0))
    bias_map = {s.source_cell_id: z for s in lineage}
    attached = attach_window_cell_to_endpoint_lineage(
        selectors,
        endpoint,
        left,
        absolute_bias_true_by_selector_cell=bias_map,
    )
    return {
        "root_valid": not _validate_cell(root),
        "binary_split_valid": not split_failures,
        "binary_split_failures": split_failures,
        "children_retain_all_parent_constraints": all(
            child.shared_constraint_tokens[: len(root.shared_constraint_tokens)] == root.shared_constraint_tokens
            for child in (left, right)
        ),
        "children_retain_same_hard_realization_witness": all(
            child.hard_realization_witness_id == root.hard_realization_witness_id
            for child in (left, right)
        ),
        "children_retain_same_joint_response_witness": all(
            child.joint_response_witness_id == root.joint_response_witness_id
            for child in (left, right)
        ),
        "selector_prefixes_attached": len(attached),
        "selector_lineage_length": len(lineage),
        "every_selector_prefix_attached_to_same_window_cell": bool(attached)
        and all(a.window_source_cell_id == left.cell_id for a in attached),
        "point_bias_attached_same_history_every_prefix": bool(attached)
        and all(a.absolute_bias_same_history for a in attached),
        "point_bias_history_is_source_uniform": False,
        "hard_provider_implementation_closed": bool(hard["provider_implementation_closed"]),
        "hard_window_realization_closed": bool(hard["finite_window_realization_certificate_closed"]),
        "R_lambda_closed": bool(rlambda["machine_readable_R_lambda_closed"]),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    hard = HARD.build(domain_path)
    hf = HARD.validate_status(hard)
    rlambda = RLAMBDA.build(domain_path)
    rf = RLAMBDA.validate(rlambda)
    selectors = SELECTORS.build(domain_path)
    sf = SELECTORS.validate(selectors)
    bad = {"hard": hf, "R_lambda": rf, "selectors": sf}
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"correlated source-cell prerequisites failed: {bad}")
    smoke = _smoke(domain_path)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "mems_bias_preconditions": BIAS.build(domain_path),
        "P3_delta_preserved": 1.0e-18,
        "complete_window_samples": 601,
        "correlated_window_source_cell_data_model_available": True,
        "binary_refinement_only_intersects_existing_source_cell": True,
        "all_coupled_parent_constraints_retained_on_children": True,
        "derived_shipping_coordinates_are_not_refinement_coordinates": True,
        "selector_lineage_window_cell_attachment_available": True,
        "absolute_A21_bias_lineage_attachment_available": True,
        "independent_true_bias_event_boxes_allowed": False,
        "source_cell_refinement_recomputes_derived_shipping_history_required": True,
        "machine_readable_R_lambda_consumed": True,
        "hard_realization_oracle_closed_here": False,
        "joint_source_output_map_closed_here": False,
        "source_uniform_absolute_bias_path_closed_here": False,
        "source_uniform_complete_window_cover_closed_here": False,
        "source_uniform_P4_graph_closed_here": False,
        "trajectory_replay_used": False,
        "finite_frequency_grid_used": False,
        "independent_per_sample_boxes_used": False,
        "filter_changed": False,
        "quality_gates_changed": False,
        "declared_domain_changed": False,
        "P4_promoted_here": False,
        "smoke": smoke,
        "next_obligation": (
            "supply a validated common-witness hard-realization/joint-output oracle for each Sea3WindowCell and a legitimate hard same-history physical A21 bias path; only then may the refinement cover be executed source-uniformly"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    f.extend(f"MEMS bias: {x}" for x in BIAS.validate(d.get("mems_bias_preconditions", {})))
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("canonical source changed")
    if float(d.get("P3_delta_preserved", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "correlated_window_source_cell_data_model_available",
        "binary_refinement_only_intersects_existing_source_cell",
        "all_coupled_parent_constraints_retained_on_children",
        "derived_shipping_coordinates_are_not_refinement_coordinates",
        "selector_lineage_window_cell_attachment_available",
        "absolute_A21_bias_lineage_attachment_available",
        "source_cell_refinement_recomputes_derived_shipping_history_required",
        "machine_readable_R_lambda_consumed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "independent_true_bias_event_boxes_allowed",
        "hard_realization_oracle_closed_here",
        "joint_source_output_map_closed_here",
        "source_uniform_absolute_bias_path_closed_here",
        "source_uniform_complete_window_cover_closed_here",
        "source_uniform_P4_graph_closed_here",
        "trajectory_replay_used",
        "finite_frequency_grid_used",
        "independent_per_sample_boxes_used",
        "filter_changed",
        "quality_gates_changed",
        "declared_domain_changed",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    smoke = d.get("smoke", {})
    for key in (
        "root_valid",
        "binary_split_valid",
        "children_retain_all_parent_constraints",
        "children_retain_same_hard_realization_witness",
        "children_retain_same_joint_response_witness",
        "every_selector_prefix_attached_to_same_window_cell",
        "point_bias_attached_same_history_every_prefix",
        "R_lambda_closed",
    ):
        if smoke.get(key) is not True:
            f.append(f"smoke lost {key}")
    if smoke.get("binary_split_failures") != []:
        f.append("binary source-cell smoke reported split failures")
    if smoke.get("selector_prefixes_attached") != smoke.get("selector_lineage_length"):
        f.append("source cell did not attach to every selector prefix")
    if smoke.get("point_bias_history_is_source_uniform") is not False:
        f.append("point A21 bias history falsely claimed source uniformity")
    if smoke.get("hard_provider_implementation_closed") is not False:
        f.append("source-cell data model falsely closed SEA0 provider")
    if smoke.get("hard_window_realization_closed") is not False:
        f.append("source-cell data model falsely closed hard window realization")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "window_cells": d["correlated_window_source_cell_data_model_available"],
        "selector_attachment": d["selector_lineage_window_cell_attachment_available"],
        "hard_oracle_closed": d["hard_realization_oracle_closed_here"],
        "source_uniform_cover": d["source_uniform_complete_window_cover_closed_here"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
