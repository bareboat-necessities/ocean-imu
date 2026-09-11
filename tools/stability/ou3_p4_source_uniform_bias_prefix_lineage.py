#!/usr/bin/env python3
"""Source-uniform absolute true-bias prefix cells for BIAS0/BIAS1/BIAS2.

Projection events need the absolute shipping-centered physical bias at the same
prefix as the nonlinear error state.  The existing joint ISS supply correctly
retains the shared driver column, but its norm bound alone is not an absolute
bias cell that the exact deployed projection can consume.

For each mandatory family this module uses that family's authoritative hard
recurrence

    b[k+1] = phi_true[k] b[k] + w[k]

with the declared component root envelope, factor interval, driver increment
bound and absolute component envelope.  It propagates one recursive interval
state per axis and intersects each successor with the family's already-declared
absolute envelope.  Thus a prefix cell is a consequence of one common root and
successive admitted recurrence steps; independent per-sample bias slots are not
introduced.  The interval recurrence is a conservative outer relation, so it
may contain histories not in the physical family, which is safe for P4.

The three families remain separate.  Their shared-w coupling to bias-error
state remains supplied by ``ou3_p4_bias_family_joint_iss_supply``; this module
never replaces that augmented-master supply by independent disturbances.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2
import ou3_p4_bias_family_joint_iss_supply as SUPPLY
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_ABSOLUTE_BIAS_PREFIX_LINEAGE_V1"
FAMILIES = ("BIAS0", "BIAS1", "BIAS2")
P3_DELTA = 1.0e-18


@dataclass(frozen=True)
class BiasPrefixLineage:
    family: str
    endpoint_source_cell_id: str
    prefix_boxes: dict[str, tuple[Interval, Interval, Interval]]
    common_root_box: tuple[Interval, Interval, Interval]
    phi_true: Interval
    driver_component: Interval
    absolute_component_cap: float
    recurrence_outer_relation: bool = True
    source_uniform: bool = True

    def at(self, source_cell_id: str) -> tuple[Interval, Interval, Interval]:
        try:
            return self.prefix_boxes[source_cell_id]
        except KeyError as exc:
            raise RuntimeError(f"bias prefix missing at {source_cell_id}") from exc


def _module(family: str):
    if family == "BIAS0":
        return BIAS0
    if family == "BIAS1":
        return BIAS1
    if family == "BIAS2":
        return BIAS2
    raise ValueError("bias family must be BIAS0/BIAS1/BIAS2")


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise RuntimeError("bias recurrence outer cell became disjoint from declared family envelope")
    return Interval(lo, hi)


def family_parameters(family: str) -> dict:
    mod = _module(family)
    d = mod.build()
    failures = mod.validate(d)
    if failures:
        raise RuntimeError(f"{family} prerequisite failed: {failures!r}")
    phi_lo, phi_hi = map(float, d["phi_true_interval"])
    cap = float(d["true_bias_component_abs_upper_mps2"])
    w = float(d["driver_increment_component_abs_upper_mps2"])
    if not (0.0 < phi_lo <= phi_hi <= 1.0 and cap > 0.0 and w >= 0.0):
        raise RuntimeError(f"{family} recurrence parameters invalid")
    return {
        "qualification": d["qualification"],
        "phi": Interval(phi_lo, phi_hi),
        "cap": cap,
        "driver": Interval(-w, w),
    }


def prefix_component_boxes(family: str, prefix_count: int) -> list[Interval]:
    """Return b[1]...b[N] outer cells from one common b[0] root box."""
    if not isinstance(prefix_count, int) or prefix_count < 1:
        raise ValueError("positive prefix count required")
    p = family_parameters(family)
    cap_box = Interval(-p["cap"], p["cap"])
    current = cap_box
    out: list[Interval] = []
    for _ in range(prefix_count):
        successor = p["phi"] * current + p["driver"]
        current = _intersect(successor, cap_box)
        out.append(current)
    return out


def attach_to_selector_lineage(
    family: str,
    lineage: Sequence[SELECTORS.PrefixSelector],
) -> BiasPrefixLineage:
    """Attach recurrence-generated absolute true-bias cells to one selector lineage."""
    if not lineage:
        raise ValueError("nonempty selector lineage required")
    for expected, selector in enumerate(lineage, start=1):
        if selector.prefix_length != expected:
            raise ValueError("selector lineage skipped or duplicated a prefix")
        if expected > 1 and selector.parent_source_cell_id != lineage[expected - 2].source_cell_id:
            raise ValueError("selector lineage ancestry broken")
    p = family_parameters(family)
    boxes = prefix_component_boxes(family, len(lineage))
    mapping = {
        selector.source_cell_id: (box, box, box)
        for selector, box in zip(lineage, boxes)
    }
    root = Interval(-p["cap"], p["cap"])
    return BiasPrefixLineage(
        family=family,
        endpoint_source_cell_id=lineage[-1].source_cell_id,
        prefix_boxes=mapping,
        common_root_box=(root, root, root),
        phi_true=p["phi"],
        driver_component=p["driver"],
        absolute_component_cap=p["cap"],
    )


def build() -> dict:
    supply = SUPPLY.build()
    sf = SUPPLY.validate(supply)
    if sf:
        raise RuntimeError(f"joint bias supply prerequisite failed: {sf!r}")
    family_rows = {}
    for family in FAMILIES:
        p = family_parameters(family)
        cells = prefix_component_boxes(family, 600)
        family_rows[family] = {
            "family_qualification": p["qualification"],
            "phi_true_interval": [p["phi"].lo, p["phi"].hi],
            "driver_component_interval_mps2": [p["driver"].lo, p["driver"].hi],
            "absolute_component_cap_mps2": p["cap"],
            "word_prefix_count": len(cells),
            "all_prefixes_inside_declared_absolute_cap": all(
                c.lo >= -p["cap"] and c.hi <= p["cap"] for c in cells
            ),
            "one_recursive_state_not_independent_slots": True,
        }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "required_bias_families": list(FAMILIES),
        "family_prefix_relations": family_rows,
        "absolute_true_bias_prefix_cells_available_for_projection": True,
        "common_root_then_recursive_successors_required": True,
        "family_absolute_cap_intersection_is_existing_hypothesis_not_new_cap": True,
        "recurrence_is_conservative_outer_relation": True,
        "independent_per_sample_absolute_bias_boxes_used": False,
        "three_bias_families_collapsed_into_generic_box": False,
        "joint_shared_w_supply_still_required": True,
        "joint_shared_w_supply_consumed": True,
        "source_uniform_projection_bias_coordinate_materialized": True,
        "production_selector_lineages_attached_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "attach each family recurrence to every retained selector lineage, feed its same-prefix absolute b_true cell into each A21 projection event, and retain the family-specific shared-w/tau-mismatch supply in the augmented LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if d.get("required_bias_families") != list(FAMILIES):
        f.append("mandatory bias family set changed")
    for k in (
        "absolute_true_bias_prefix_cells_available_for_projection",
        "common_root_then_recursive_successors_required",
        "family_absolute_cap_intersection_is_existing_hypothesis_not_new_cap",
        "recurrence_is_conservative_outer_relation",
        "joint_shared_w_supply_still_required",
        "joint_shared_w_supply_consumed",
        "source_uniform_projection_bias_coordinate_materialized",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "independent_per_sample_absolute_bias_boxes_used",
        "three_bias_families_collapsed_into_generic_box",
        "production_selector_lineages_attached_here",
        "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here",
        "P4_MOTION_PASS",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    rows = d.get("family_prefix_relations", {})
    if set(rows) != set(FAMILIES):
        f.append("family prefix relation table incomplete")
    for family in FAMILIES:
        row = rows.get(family, {})
        if row.get("word_prefix_count") != 600:
            f.append(f"{family} prefix count changed")
        if row.get("all_prefixes_inside_declared_absolute_cap") is not True:
            f.append(f"{family} recurrence escaped declared cap")
        if row.get("one_recursive_state_not_independent_slots") is not True:
            f.append(f"{family} recurrence lost common-state semantics")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "projection_bias_prefix_cells": d["source_uniform_projection_bias_coordinate_materialized"],
        "production_attached": d["production_selector_lineages_attached_here"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
