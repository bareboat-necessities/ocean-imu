#!/usr/bin/env python3
"""Fail-closed numerical gate from V13E radial closure to V14D q<8 closure.

The V14D producer is defined only after the V13E signed/radial prerequisite has
closed.  Its historical CLI assumed the post-prerequisite fields were present
and could therefore raise KeyError when a deliberately coarse V13E cover was
still NOT_ESTABLISHED.  That is a reporting bug, not proof evidence.

This gate evaluates V13E first.  If V13E is not established it emits the exact
radial witness and stops successfully as a proof-status diagnostic.  Only when
V13E passes does it invoke V14D.  Numerical NOT_ESTABLISHED is kept distinct
from validator/infrastructure failure, and neither path promotes P5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_sample1_signed_radial_subcells_v13e as V13E
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D

DEFAULT_DOMAIN = V14D.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 2, parallel_pieces: int = 2) -> dict:
    path = Path(domain_path).resolve()
    kwargs = dict(
        source_pieces=source_pieces,
        source_cell_index=source_cell_index,
        p_pieces=p_pieces,
        tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces,
    )
    radial = V13E.build(path, **kwargs)
    radial_validation = V13E.validate(radial)
    radial_pass = (
        radial.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E") == "PASS"
        and not radial_validation
    )

    if not radial_pass:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V13E_TO_V14D_FAIL_CLOSED_GATE",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "filter_changed": False,
            "source_pieces": int(source_pieces),
            "source_cell_index": int(source_cell_index),
            "p_pieces": int(p_pieces),
            "tangent_pieces": int(tangent_pieces),
            "axial_pieces": int(axial_pieces),
            "residual_x_pieces": int(residual_x_pieces),
            "parallel_pieces": int(parallel_pieces),
            "V13E_status": radial.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E"),
            "V13E_validation_failures": radial_validation,
            "V13E_evaluated_signed_subcells": radial.get("evaluated_signed_subcells"),
            "V13E_above_6rad_subcells": radial.get("above_6rad_subcells"),
            "V13E_unclosed_radial_subcells": radial.get("unclosed_radial_subcells"),
            "V13E_max_radial_upper": radial.get("max_radial_upper"),
            "V13E_minimum_radial_lower_above_6": radial.get("minimum_radial_lower_above_6"),
            "V13E_first_unclosed_radial_subcell": radial.get("first_unclosed_radial_subcell"),
            "V13E_worst_radial_subcell": radial.get("worst_radial_subcell"),
            "V14D_invoked": False,
            "V14D_status": "BLOCKED_BY_V13E",
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_V13E_TO_V14D_GATE": "NOT_ESTABLISHED",
            "next_obligation": "REFINE_SIGNED_RX_PARALLEL_SUBDIVISION_AT_FIRST_V13E_RADIAL_WITNESS",
            "failures": list(radial_validation),
        }

    q8 = V14D.build(path, **kwargs)
    q8_validation = V14D.validate(q8)
    q8_pass = (
        q8.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D") == "PASS"
        and not q8_validation
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V13E_TO_V14D_FAIL_CLOSED_GATE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_pieces": int(source_pieces),
        "source_cell_index": int(source_cell_index),
        "p_pieces": int(p_pieces),
        "tangent_pieces": int(tangent_pieces),
        "axial_pieces": int(axial_pieces),
        "residual_x_pieces": int(residual_x_pieces),
        "parallel_pieces": int(parallel_pieces),
        "V13E_status": "PASS",
        "V13E_validation_failures": [],
        "V13E_evaluated_signed_subcells": radial.get("evaluated_signed_subcells"),
        "V13E_above_6rad_subcells": radial.get("above_6rad_subcells"),
        "V13E_unclosed_radial_subcells": radial.get("unclosed_radial_subcells"),
        "V13E_max_radial_upper": radial.get("max_radial_upper"),
        "V13E_minimum_radial_lower_above_6": radial.get("minimum_radial_lower_above_6"),
        "V13E_first_unclosed_radial_subcell": None,
        "V13E_worst_radial_subcell": radial.get("worst_radial_subcell"),
        "V14D_invoked": True,
        "V14D_status": q8.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D"),
        "V14D_validation_failures": q8_validation,
        "V14D_evaluated_signed_cayley_cells": q8.get("evaluated_signed_cayley_cells"),
        "V14D_product_scalar_antipode_cells": q8.get("product_scalar_antipode_cells"),
        "V14D_unclosed_q8_cells": q8.get("unclosed_q8_cells"),
        "V14D_minimum_abs_product_scalar_lower": q8.get("minimum_abs_product_scalar_lower"),
        "V14D_max_post_sample1_cayley_norm_upper": q8.get("max_post_sample1_cayley_norm_upper"),
        "V14D_first_unclosed_q8_cell": q8.get("first_unclosed_q8_cell"),
        "V14D_worst_q8_cell": q8.get("worst_q8_cell"),
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V13E_TO_V14D_GATE": "PASS" if q8_pass else "NOT_ESTABLISHED",
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if q8_pass else
            "REFINE_CURRENT_CX_OR_CORRECTION_RADIAL_DIRECTION_AT_FIRST_Q8_WITNESS"
        ),
        "failures": list(q8_validation),
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V13E_TO_V14D_FAIL_CLOSED_GATE":
        failures.append("qualification mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("source-generated contract lost")
    for key in ("source_replay_used", "filter_changed", "q8_word_promoted_here",
                "whole_word_promoted_here", "N_H_words_set_here"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")

    radial_pass = d.get("V13E_status") == "PASS" and not d.get("V13E_validation_failures", [])
    if not radial_pass:
        if d.get("V14D_invoked") is not False:
            failures.append("V14D invoked before V13E closure")
        if d.get("V14D_status") != "BLOCKED_BY_V13E":
            failures.append("V14D blocking status mismatch")
        if d.get("P5_SAMPLE1_V13E_TO_V14D_GATE") != "NOT_ESTABLISHED":
            failures.append("gate passed without V13E closure")
        if d.get("V13E_first_unclosed_radial_subcell") is None and not d.get("V13E_validation_failures", []):
            failures.append("V13E numerical nonclosure lacks radial witness")
    else:
        if d.get("V14D_invoked") is not True:
            failures.append("V14D not invoked after V13E closure")
        q8_pass = d.get("V14D_status") == "PASS" and not d.get("V14D_validation_failures", [])
        if (d.get("P5_SAMPLE1_V13E_TO_V14D_GATE") == "PASS") != q8_pass:
            failures.append("gate status does not match V14D closure")
        if not q8_pass and d.get("V14D_first_unclosed_q8_cell") is None and not d.get("V14D_validation_failures", []):
            failures.append("V14D numerical nonclosure lacks q8 witness")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=2)
    ap.add_argument("--parallel-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain,
        source_pieces=args.source_pieces,
        source_cell_index=args.source_cell_index,
        p_pieces=args.p_pieces,
        tangent_pieces=args.tangent_pieces,
        axial_pieces=args.axial_pieces,
        residual_x_pieces=args.residual_x_pieces,
        parallel_pieces=args.parallel_pieces,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_V13E_TO_V14D_GATE"],
        "V13E_status": out["V13E_status"],
        "V13E_unclosed": out.get("V13E_unclosed_radial_subcells"),
        "V13E_first_unclosed": out.get("V13E_first_unclosed_radial_subcell"),
        "V13E_worst": out.get("V13E_worst_radial_subcell"),
        "V14D_invoked": out["V14D_invoked"],
        "V14D_status": out["V14D_status"],
        "V14D_unclosed": out.get("V14D_unclosed_q8_cells"),
        "V14D_first_unclosed": out.get("V14D_first_unclosed_q8_cell"),
        "V14D_worst": out.get("V14D_worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    # Numerical NOT_ESTABLISHED is a proof result, not a CI/infrastructure error.
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
