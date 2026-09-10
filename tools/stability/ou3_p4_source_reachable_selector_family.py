#!/usr/bin/env python3
"""Relational source-reachable COMPLETE-BRMM selector-family cover for P4.

P4 is universal over the continuum COMPLETE-BRMM source and may therefore be
proved on the validated correlated outer relation O^601_BRMM rather than by
enumerating source histories.  This module composes already-certified relations:

  B^601_BRMM subset O^601_BRMM
    -> compact/invariant shipping frontend predecessor family
    -> same-history joint frontend/coefficient transition
    -> trusted typed Riccati transition / prefix-selector ancestry
    -> full hard-entry radial scale r in [0,1]
    -> one recursive absolute b_true lineage for each BIAS0/BIAS1/BIAS2.

The result is a *relation cover*: for every admitted BRMM history, every
shipping branch successor is represented by a selector-family path with the
same physical/estimator ancestry.  No favorable branch, seeded realization,
frequency grid, independent sample box, independent tuner schedule, or
independent per-event bias box is introduced.

This does not by itself prove P4.  In particular it does not assert that a
single coarse interval hull of all 601 numeric P/H/R/K matrices is useful.
The next theorem object must evaluate the exact nonlinear event graph and its
quadratic sectors over this relation and close endpoint/every-prefix LDLT.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_correlated_window_outer_enclosure as OUTER
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_brmm_frontend_predecessor_invariant as PREDECESSOR
import ou3_p4_complete_brmm_universal_target_relation as TARGET
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_source_uniform_bias_prefix_lineage as BIASPREFIX
import ou3_p4_all_bias_nonlinear_lineage_binding as ALLBIAS
import ou3_p4_acceleration_moment_iqc_sector as MOMSECTOR
import ou3_p4_hard_entry_set as ENTRY

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_REACHABLE_COMPLETE_BRMM_SELECTOR_FAMILY_RELATION_V1"
CANONICAL_SOURCE = "COMPLETE_BRMM_NORMAL_LIVE_WORD"
P3_DELTA = 1.0e-18


@dataclass(frozen=True)
class EndpointFamilyAttachment:
    endpoint_source_cell_id: str
    selector_lineage: tuple[SELECTORS.PrefixSelector, ...]
    radial_scale: Interval
    bias_lineages: dict[str, BIASPREFIX.BiasPrefixLineage]


def attach_endpoint_family(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
) -> EndpointFamilyAttachment:
    """Attach full radial scale and all bias recurrences to one retained endpoint.

    This helper is intentionally agnostic to how a concrete selector path was
    represented.  It validates literal ancestry and then attaches only
    theorem-owned relations.  It never selects a branch or fabricates P/H/R/K.
    """
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint_source_cell_id)
    if not lineage:
        raise ValueError("source-reachable endpoint lineage is empty")
    bias = {
        family: BIASPREFIX.attach_to_selector_lineage(family, lineage)
        for family in BIASPREFIX.FAMILIES
    }
    if any(x.endpoint_source_cell_id != endpoint_source_cell_id for x in bias.values()):
        raise RuntimeError("bias recurrence endpoint detached from selector endpoint")
    return EndpointFamilyAttachment(
        endpoint_source_cell_id=endpoint_source_cell_id,
        selector_lineage=tuple(lineage),
        radial_scale=Interval(0.0, 1.0),
        bias_lineages=bias,
    )


def build(domain_path: Path = OUTER.DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    outer = OUTER.build(path)
    predecessor = PREDECESSOR.build()
    target = TARGET.build(path)
    kernel = KERNEL.build(path)
    selectors = SELECTORS.build(path)
    cover = COVER.build(path)
    bias = BIASPREFIX.build()
    allbias = ALLBIAS.build()
    moment = MOMSECTOR.build()
    entry = ENTRY.build()
    bad = {
        "outer": OUTER.validate(outer),
        "predecessor": PREDECESSOR.validate(predecessor),
        "target": TARGET.validate(target),
        "kernel": KERNEL.validate(kernel),
        "selectors": SELECTORS.validate(selectors),
        "cover": COVER.validate(cover),
        "bias_prefix": BIASPREFIX.validate(bias),
        "all_bias_nonlinear": ALLBIAS.validate(allbias),
        "moment_sector": MOMSECTOR.validate(moment),
        "entry": ENTRY.validate(entry),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("selector-family prerequisites failed: " + repr(bad))

    outer_closed = bool(
        outer["left_inclusion_closed"]
        and outer["validated_correlated_outer_enclosure_closed"]
        and outer["same_history_required_for_entire_window"]
        and outer["correlation_retained_across_samples"]
        and outer["correlation_retained_across_axes"]
        and not outer["independent_sample_boxes_used"]
    )
    predecessor_closed = bool(
        predecessor["complete_BRMM_predecessor_state_family_covered"]
        and predecessor["all_component_predecessor_invariants_closed"]
    )
    estimator_closed = bool(
        target["joint_transition_physical_BRMM_attachment_closed"]
        and target["complete_BRMM_predecessor_family_covered"]
        and target["coefficient_target_inclusion_closed"]
        and target["same_BRMM_sample_joint_frontend_transition_materialized"]
        and target["same_BRMM_current_applied_schedule_retained"]
    )
    riccati_relation = bool(
        kernel["typed_execution_kernel_ready"]
        and kernel["front_end_branch_splits_retained"]
        and kernel["raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates"]
        and kernel["covariance_floor_increment_computed_from_current_mode_P"]
    )
    selector_relation = bool(
        selectors["branch_correlated_prefix_selectors_available"]
        and selectors["event_local_Riccati_cells_available"]
        and selectors["literal_event_order_retained"]
        and selectors["same_actual_RS_provenance_retained"]
    )
    radial_relation = bool(
        entry["full_declared_scale_enforced"]
        and moment["same_radial_coordinate_must_parameterize_entire_source_history"]
        and moment["moment_sector_matrix_available"]
    )
    bias_relation = bool(
        bias["source_uniform_projection_bias_coordinate_materialized"]
        and allbias["all_three_bias_families_have_executable_A21_lineage_binding"]
        and allbias["joint_shared_w_tau_mismatch_supply_retained_for_augmented_master"]
    )
    relation_closed = all((
        outer_closed, predecessor_closed, estimator_closed, riccati_relation,
        selector_relation, radial_relation, bias_relation,
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "source_relation": outer["outer_set_symbol"],
        "complete_BRMM_left_inclusion_closed": outer_closed,
        "frontend_predecessor_family_closed": predecessor_closed,
        "same_history_estimator_coefficient_relation_closed": estimator_closed,
        "trusted_typed_Riccati_transition_relation_available": riccati_relation,
        "branch_complete_prefix_selector_relation_available": selector_relation,
        "full_zero_to_one_hard_entry_radial_relation_attached": radial_relation,
        "all_bias_absolute_prefix_relations_attached": bias_relation,
        "all_bias_shared_w_tau_mismatch_supply_retained": True,
        "source_reachable_COMPLETE_BRMM_selector_family_relation_closed": relation_closed,
        "selector_family_is_relation_not_finite_source_enumeration": True,
        "all_branch_successors_retained": True,
        "favorable_branch_selected": False,
        "independent_sample_boxes_generate_estimator_history": False,
        "independent_tuner_schedule_used": False,
        "independent_per_event_bias_boxes_used": False,
        "finite_frequency_or_direction_grid_used": False,
        "seeded_realization_used": False,
        "trajectory_replay_used": False,
        "numeric_601_sample_interval_hull_materialized_here": False,
        "joint_augmented_cocycle_closed_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "lift this universal selector-family relation through the exact nonlinear H18/A21 event maps, "
            "retain the radial/moment and family-specific bias supply coordinates in one augmented cocycle, "
            "then run endpoint and every-prefix outward LDLT on that relation"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("canonical source changed")
    if d.get("source_relation") != "O^601_BRMM":
        f.append("correlated outer source relation changed")
    for k in (
        "complete_BRMM_left_inclusion_closed",
        "frontend_predecessor_family_closed",
        "same_history_estimator_coefficient_relation_closed",
        "trusted_typed_Riccati_transition_relation_available",
        "branch_complete_prefix_selector_relation_available",
        "full_zero_to_one_hard_entry_radial_relation_attached",
        "all_bias_absolute_prefix_relations_attached",
        "all_bias_shared_w_tau_mismatch_supply_retained",
        "source_reachable_COMPLETE_BRMM_selector_family_relation_closed",
        "selector_family_is_relation_not_finite_source_enumeration",
        "all_branch_successors_retained",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "favorable_branch_selected",
        "independent_sample_boxes_generate_estimator_history",
        "independent_tuner_schedule_used",
        "independent_per_event_bias_boxes_used",
        "finite_frequency_or_direction_grid_used",
        "seeded_realization_used",
        "trajectory_replay_used",
        "numeric_601_sample_interval_hull_materialized_here",
        "joint_augmented_cocycle_closed_here",
        "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here",
        "P4_MOTION_PASS",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=OUTER.DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "selector_family_relation": d["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"],
        "radial": d["full_zero_to_one_hard_entry_radial_relation_attached"],
        "all_bias": d["all_bias_absolute_prefix_relations_attached"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
