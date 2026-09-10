#!/usr/bin/env python3
"""Bind BIAS0/BIAS1/BIAS2 recursive true-bias prefixes into nonlinear A21.

The existing same-history nonlinear graph consumer requires an absolute
``b_true`` cell at each selector child because the deployed A21 projection is a
map of both estimation error and the absolute physical bias.  Previously its
only executable bias lineage was the homogeneous zero/BIAS1 point helper.

This module adapts the source-uniform recursive outer relations from
``ou3_p4_source_uniform_bias_prefix_lineage`` to that consumer.  Every A21
measurement on a retained selector lineage receives the bias cell generated at
the same prefix by one common-root recurrence.  BIAS0, BIAS1 and BIAS2 are
executed separately; no family is inferred from another and no independent
per-event bias boxes are accepted.

The nonlinear interval evaluation may forget some temporal correlation inside
a finite event Jacobian enclosure, which only enlarges the covered set.  The
shared-w/tau-mismatch relation across predictions is still retained separately
by the joint ISS supply and must be included in the augmented master.  Hence
this bridge closes the absolute-bias/projection attachment, not the endpoint or
every-prefix LDLT and not P4.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping, Sequence

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_bias_family_joint_iss_supply as SUPPLY
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_p4_source_uniform_bias_prefix_lineage as PREFIX

SCHEMA = 1
QUALIFICATION = "OU3_P4_ALL_BIAS_SOURCE_UNIFORM_NONLINEAR_LINEAGE_BINDING_V1"
P3_DELTA = 1.0e-18


@dataclass(frozen=True)
class GraphBiasLineageAdapter:
    family: str
    endpoint_source_cell_id: str
    bias_true_by_source_cell_id: Mapping[str, tuple[Interval, Interval, Interval]]
    prefix_certificate: PREFIX.BiasPrefixLineage
    source_uniform_materialization: bool = True

    def at(self, source_cell_id: str) -> tuple[Interval, Interval, Interval]:
        try:
            return self.bias_true_by_source_cell_id[source_cell_id]
        except KeyError as exc:
            raise RuntimeError(f"{self.family} absolute true-bias prefix missing at {source_cell_id}") from exc

    def validate_homogeneous(self, lineage, constants, projection_limit):
        """Graph hook: validate recursive family lineage, never claim homogeneity."""
        if self.source_uniform_materialization is not True:
            raise RuntimeError("all-family graph adapter must be source-uniform")
        if self.prefix_certificate.family != self.family:
            raise RuntimeError("bias family certificate detached from graph adapter")
        if self.prefix_certificate.endpoint_source_cell_id != self.endpoint_source_cell_id:
            raise RuntimeError("bias endpoint certificate detached from graph lineage")
        ids = [selector.source_cell_id for selector in lineage]
        if not ids or ids[-1] != self.endpoint_source_cell_id:
            raise RuntimeError("bias adapter attached to a different selector endpoint")
        if set(ids) != set(self.bias_true_by_source_cell_id):
            raise RuntimeError("bias adapter must cover every and only selector prefix")
        for selector in lineage:
            if self.at(selector.source_cell_id) != self.prefix_certificate.at(selector.source_cell_id):
                raise RuntimeError("absolute true-bias cell detached from recursive prefix certificate")
        # This adapter intentionally does not require phi_true==phi_hat.  Tau
        # mismatch is a retained supply in the augmented master, not erased.
        supply = SUPPLY.build()["family_supply"][self.family]
        plo, phi = supply["phi_true_interval"]
        cert = self.prefix_certificate.phi_true
        if cert.lo != float(plo) or cert.hi != float(phi):
            raise RuntimeError("bias recurrence factor detached from family joint ISS supply")
        if float(projection_limit) <= 0.0:
            raise RuntimeError("deployed projection radius missing")


def adapter_for_lineage(
    family: str,
    lineage: Sequence[SELECTORS.PrefixSelector],
) -> GraphBiasLineageAdapter:
    cert = PREFIX.attach_to_selector_lineage(family, lineage)
    return GraphBiasLineageAdapter(
        family=family,
        endpoint_source_cell_id=cert.endpoint_source_cell_id,
        bias_true_by_source_cell_id=cert.prefix_boxes,
        prefix_certificate=cert,
    )


def consume_A21_endpoint_all_families(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
    *,
    initial_A_state: Sequence[Interval],
    domain_path: Path = GRAPH.DEFAULT_DOMAIN,
) -> dict[str, GRAPH.NonlinearLineageResult]:
    """Compose the exact A21 lineage separately for every mandatory bias family."""
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint_source_cell_id)
    if not lineage:
        raise ValueError("endpoint nonlinear lineage is empty")
    constants = KERNEL._process_constants(domain_path)
    limit = GRAPH._projection_limit(domain_path)
    out: dict[str, GRAPH.NonlinearLineageResult] = {}
    for family in PREFIX.FAMILIES:
        bias = adapter_for_lineage(family, lineage)
        out[family] = GRAPH._consume_mode_lineage(
            mode="A",
            lineage=lineage,
            initial_state=initial_A_state,
            source_token=f"{GRAPH.CANONICAL_SOURCE}:{endpoint_source_cell_id}:{family}",
            constants=constants,
            bias_lineage=bias,
            projection_limit=limit,
        )
    return out


def build() -> dict:
    prefix = PREFIX.build()
    pf = PREFIX.validate(prefix)
    supply = SUPPLY.build()
    sf = SUPPLY.validate(supply)
    graph = GRAPH.build()
    gf = GRAPH.validate(graph)
    bad = {"prefix": pf, "supply": sf, "graph": gf}
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("all-bias nonlinear binding prerequisites failed: " + repr(bad))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": GRAPH.CANONICAL_SOURCE,
        "required_bias_families": list(PREFIX.FAMILIES),
        "source_uniform_absolute_bias_prefix_relation_consumed": True,
        "same_prefix_bias_cell_feeds_A21_projection": True,
        "all_three_bias_families_have_executable_A21_lineage_binding": True,
        "tau_mismatch_not_erased_by_graph_adapter": True,
        "joint_shared_w_tau_mismatch_supply_retained_for_augmented_master": True,
        "homogeneous_BIAS1_assumption_used_for_all_families": False,
        "independent_per_event_true_bias_boxes_used": False,
        "bias_families_collapsed": False,
        "production_complete_BRMM_selector_family_materialized_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "apply this all-family A21 binding to every source-reachable COMPLETE-BRMM selector lineage, combine its family-specific shared-w/tau-mismatch supply with the nonlinear cocycle, then close endpoint/every-prefix augmented LDLT and first-exit retention"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != GRAPH.CANONICAL_SOURCE:
        f.append("canonical source changed")
    if d.get("required_bias_families") != list(PREFIX.FAMILIES):
        f.append("mandatory bias family set changed")
    for k in (
        "source_uniform_absolute_bias_prefix_relation_consumed",
        "same_prefix_bias_cell_feeds_A21_projection",
        "all_three_bias_families_have_executable_A21_lineage_binding",
        "tau_mismatch_not_erased_by_graph_adapter",
        "joint_shared_w_tau_mismatch_supply_retained_for_augmented_master",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "homogeneous_BIAS1_assumption_used_for_all_families",
        "independent_per_event_true_bias_boxes_used",
        "bias_families_collapsed",
        "production_complete_BRMM_selector_family_materialized_here",
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
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "all_bias_A21_binding": d["all_three_bias_families_have_executable_A21_lineage_binding"],
        "provider_family": d["production_complete_BRMM_selector_family_materialized_here"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
