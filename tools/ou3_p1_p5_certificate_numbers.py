#!/usr/bin/env python3
"""Consolidated P1--P5 certificate numbers and usability verdicts for OU-III.

``docs/ou-iii-certificate-usability-envelope.md`` states what *usable* means for
each stage of the deployed OU-III stability proof.  Until now the numbers behind
that prose lived in eight separate producers, so answering "what are the current
certificate numbers, and are they usable?" meant running them one at a time and
reading fields out of eight JSON files.

This producer recomputes the stages and emits one report: per stage, the
headline physical numbers, the explicit usability checks with their thresholds,
and a per-stage verdict.  It does **not** trust an upstream ``PASS`` field as a
usability statement -- every threshold is re-applied here against the recomputed
number, in the spirit of the final composition gate.

Two stage obligations are open by construction and are reported as such:
the complete-word P4 full-state dissipation and the P5 finite capture into the
inner stochastic localization level.  A stage whose geometry is usable but whose
theorem obligation is open is reported as ``USABLE_GEOMETRY_OPEN_OBLIGATION``,
never as complete.

The slowest inputs -- the rigorous complete-word translation block and the
signed first-accelerometer sector -- may be supplied as pre-computed JSON so
that a report can be produced without re-running them.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_first_accel_sector_budget as BUDGET
import ou3_p4_nonlinear_word_certificate as WORD
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_route_ceiling_certificate as CEILING
import ou3_p4_source_path_reachability as P2
import ou3_p5_outer_sector_capture_certificate as P5
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

# Minimum usability contracts.  These mirror the "Quantitative outer target"
# and "Non-regression rules" sections of the usability envelope note.  CI
# rejects regression below them; a stronger proof may raise them.
MIN_NORMAL_HANDOFF_CAYLEY = 0.20
MIN_TIMEOUT_HANDOFF_CAYLEY = 0.50
MIN_HANDOFF_POSITION_ERROR_M = 20.0
MIN_HANDOFF_VELOCITY_ERROR_MPS = 5.0
MIN_HANDOFF_AW_ERROR_MPS2 = 2.9
MIN_SECTOR_ANGLE_RAD = 0.80
MIN_SECTOR_MONOTONICITY = 0.80
MAX_SECTOR_ETA_RATIO = 0.25
MIN_P5_ENTRANCE_ANGLE_DEG = 45.0
MICROSCOPIC_CAYLEY_THRESHOLD = 1.0e-6


def _check(name: str, value, threshold, ok: bool, note: str = "") -> dict:
    row = {"check": name, "value": value, "threshold": threshold, "pass": bool(ok)}
    if note:
        row["note"] = note
    return row


def _verdict(checks: list[dict], open_obligation: str | None) -> str:
    if not all(c["pass"] for c in checks):
        return "NOT_USABLE"
    return "USABLE_GEOMETRY_OPEN_OBLIGATION" if open_obligation else "USABLE"


def _p1_stage(d: dict) -> dict:
    normal = float(d["P1_normal_gauged_cayley_norm_upper"])
    timeout = float(d["P1_timeout_gauged_cayley_norm_upper"])
    box = d["startup"]["operating_domain"]["physical_handoff_coordinate_bounds"]
    checks = [
        _check("startup certificate validates", d["startup_validation_pass"], True,
               d["startup_validation_pass"] is True),
        _check("normal gauged handoff Cayley norm", normal, MIN_NORMAL_HANDOFF_CAYLEY,
               normal >= MIN_NORMAL_HANDOFF_CAYLEY),
        _check("timeout gauged handoff Cayley norm", timeout, MIN_TIMEOUT_HANDOFF_CAYLEY,
               timeout >= MIN_TIMEOUT_HANDOFF_CAYLEY),
        _check("handoff position error bound (m)",
               box["position_error_norm_upper_m"], MIN_HANDOFF_POSITION_ERROR_M,
               float(box["position_error_norm_upper_m"]) >= MIN_HANDOFF_POSITION_ERROR_M),
        _check("handoff velocity error bound (m/s)",
               box["velocity_error_norm_upper_mps"], MIN_HANDOFF_VELOCITY_ERROR_MPS,
               float(box["velocity_error_norm_upper_mps"]) >= MIN_HANDOFF_VELOCITY_ERROR_MPS),
        _check("handoff latent-acceleration error bound (m/s^2)",
               box["latent_acceleration_error_norm_upper_mps2"], MIN_HANDOFF_AW_ERROR_MPS2,
               float(box["latent_acceleration_error_norm_upper_mps2"]) >= MIN_HANDOFF_AW_ERROR_MPS2),
        _check("handoff family is not microscopic", min(normal, timeout),
               f"> {MICROSCOPIC_CAYLEY_THRESHOLD}",
               min(normal, timeout) > MICROSCOPIC_CAYLEY_THRESHOLD),
    ]
    return {
        "stage": "P1",
        "title": "startup / reset / Live handoff",
        "status": "PASS" if d["startup_validation_pass"] else "NOT_ESTABLISHED",
        "headline": {
            "normal_gauged_handoff_cayley_norm_upper": normal,
            "timeout_gauged_handoff_cayley_norm_upper": timeout,
            "ungauged_timeout_route": "gravity-direction yaw quotient",
            "physical_handoff_coordinate_bounds": box,
            "declared_startup_aw_error_fraction_g":
                d["startup"]["operating_domain"]["latent_acceleration_error_fraction_g"],
        },
        "usability_checks": checks,
        "open_obligation": None,
        "verdict": _verdict(checks, None),
    }


def _p2_stage(d: dict) -> dict:
    p2 = d["P2"]
    states = int(p2["partition"]["states"])
    recurrent = int(p2["recurrent_states"])
    checks = [
        _check("source-path certificate", p2["P2_SOURCE_PATH_CERTIFICATE"], "PASS",
               p2["P2_SOURCE_PATH_CERTIFICATE"] == "PASS"),
        _check("source-only language", p2["source_only"], True, p2["source_only"] is True),
        _check("no replay used", p2["trajectory_replay_used"], False,
               p2["trajectory_replay_used"] is False),
        _check("partition states", states, "> 0", states > 0),
        _check("every state recurrent", recurrent, states, recurrent == states),
        _check("raw tuner sigma subfloor states retained",
               p2["raw_tuner_sigma_subfloor_states_included"], True,
               p2["raw_tuner_sigma_subfloor_states_included"] is True),
        _check("R_S target uses full deployed clamp",
               p2["RS_target_full_deployed_clamp_overapprox"], True,
               p2["RS_target_full_deployed_clamp_overapprox"] is True),
    ]
    return {
        "stage": "P2",
        "title": "complete source-word language / path graph",
        "status": p2["P2_SOURCE_PATH_CERTIFICATE"],
        "headline": {
            "partition_states": states,
            "transition_edges": p2["transition_edges"],
            "strongly_connected_components": p2["strongly_connected_components"],
            "recurrent_states": recurrent,
            "filter_sigma_floor_mps2": p2["filter_sigma_floor_mps2"],
            "raw_tuner_sigma_partition_lower": p2["raw_tuner_sigma_partition_lower"],
        },
        "usability_checks": checks,
        "open_obligation": None,
        "verdict": _verdict(checks, None),
    }


def _p3_stage(d: dict) -> dict:
    word = d["word"]
    modes = {m: word["modes"][m] for m in ("H", "A") if m in word.get("modes", {})}
    deltas = {m: float(v["P3_word_endpoint_delta_lower"]) for m, v in modes.items()}
    prefix = {m: float(v["P3_homogeneous_prefix_information_gain_upper"])
              for m, v in modes.items()}
    checks = [
        _check("H endpoint information margin", deltas.get("H"), "> 0",
               deltas.get("H", 0.0) > 0.0),
        _check("A endpoint information margin", deltas.get("A"), "> 0",
               deltas.get("A", 0.0) > 0.0),
        _check("H prefix information gain bound", prefix.get("H"), "<= 1",
               prefix.get("H", math.inf) <= 1.0),
        _check("A prefix information gain bound", prefix.get("A"), "<= 1",
               prefix.get("A", math.inf) <= 1.0),
    ]
    return {
        "stage": "P3",
        "title": "validated H/A linear information-word certificate",
        "status": "PASS" if all(c["pass"] for c in checks) else "NOT_ESTABLISHED",
        "headline": {
            "endpoint_information_margin_lower": deltas,
            "prefix_information_gain_upper": prefix,
        },
        "interpretation": (
            "the endpoint margin is a relative Riccati/noise comparison constant;"
            " it is not a nonlinear state radius and must not be advertised as one"
        ),
        "usability_checks": checks,
        "open_obligation": None,
        "verdict": _verdict(checks, None),
    }


def _p4_stage(d: dict) -> dict:
    s = d["P4_sector"]
    ceiling = d["route_ceiling"]
    budget = d["P4_first_accel_budget"]
    angle = float(s["design_full_attitude_angle_rad"])
    q = float(s["design_cayley_norm_upper"])
    mono = float(s["exact_vector_strong_monotonicity_factor_lower"])
    eta = float(s["exact_eta_to_rotational_residual_information_ratio_upper"])
    checks = [
        _check("finite-angle sector certificate",
               s["P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE"], "PASS",
               s["P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE"] == "PASS"),
        _check("design sector angle (rad)", angle, MIN_SECTOR_ANGLE_RAD,
               angle >= MIN_SECTOR_ANGLE_RAD),
        _check("Cayley chart stays below one", q, "< 1", q < 1.0),
        _check("exact vector strong monotonicity", mono, MIN_SECTOR_MONOTONICITY,
               mono > MIN_SECTOR_MONOTONICITY),
        _check("exact eta / residual information ratio", eta, MAX_SECTOR_ETA_RATIO,
               eta < MAX_SECTOR_ETA_RATIO),
        _check("normal P1 handoff inside sector",
               s["P1_overlap"]["normal_gauged_inside_sector"], True,
               s["P1_overlap"]["normal_gauged_inside_sector"] is True),
        _check("timeout P1 handoff inside sector",
               s["P1_overlap"]["timeout_gauged_inside_sector"], True,
               s["P1_overlap"]["timeout_gauged_inside_sector"] is True),
    ]
    retired = {m: {
        "certified_attitude_capture_radius_now":
            ceiling["modes"][m]["certified_attitude_capture_radius_now"],
        "route_ceiling_at_shipping_prefix_factor":
            ceiling["modes"][m]["route_ceiling_at_shipping_prefix_factor"],
        "route_can_reach_P1_handoff":
            ceiling["modes"][m]["route_can_reach_P1_handoff"],
    } for m in ceiling["modes"]}
    head = {
        "design_full_attitude_angle_rad": angle,
        "design_full_attitude_angle_deg": s["design_full_attitude_angle_deg"],
        "design_cayley_norm_upper": q,
        "exact_vector_strong_monotonicity_factor_lower": mono,
        "exact_eta_to_rotational_residual_information_ratio_upper": eta,
        "retired_uniform_transport_route": retired,
        "first_accelerometer_sector_budget": {
            "declared_aw_over_lowest_specific_force":
                budget["declared_aw_over_lowest_specific_force"],
            "ladder": [{
                "angle_deg": r["angle_deg"],
                "family": r["ladder_family"],
                "budget_upper_rad": r["sector_invariance_correction_budget_upper_rad"],
                "nuisance_upper_rad": r["nuisance_correction_norm_upper_rad"],
                "nuisance_over_budget_ratio": r["nuisance_over_budget_ratio"],
            } for r in budget["ladder_rows"]],
            "smallest_probe_angle_deg": budget["smallest_probe_angle_deg"],
            "smallest_probe_nuisance_over_budget_ratio":
                budget["smallest_probe_nuisance_over_budget_ratio"],
            "shrinking_the_candidate_angle_alone_can_close_the_budget":
                budget["shrinking_the_candidate_angle_alone_can_close_the_budget"],
        },
    }
    if d.get("first_accel") is not None:
        fa = d["first_accel"]
        head["signed_first_accelerometer_30deg"] = {
            "status": fa.get("P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE"),
            "family_completely_evaluated":
                fa.get("first_accelerometer_family_completely_evaluated"),
            "evaluated_children": fa.get("evaluated_children"),
            "max_correction_norm_upper_rad": fa.get("max_correction_norm_upper_rad"),
            "max_post_update_q_upper":
                fa.get("max_accepted_or_rejected_post_update_q_upper"),
            "minimum_composition_denominator_lower":
                fa.get("minimum_signed_composition_denominator_lower"),
        }
    if d.get("translation") is not None:
        tr = d["translation"]
        head["complete_word_translation_block"] = {
            "status": tr.get("P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS"),
            "modes": {m: {
                "complete_word_translation_margin_lower":
                    tr["modes"][m]["complete_word_translation_margin_lower"],
                "margin_widening_factor_lower":
                    tr["modes"][m]["margin_widening_factor_lower"],
            } for m in tr.get("modes", {})},
        }
    return {
        "stage": "P4",
        "title": "nonlinear finite-angle sector and complete-word dissipation",
        "status": s["P4_OPERATION_MATCHED_FINITE_ANGLE_SECTOR_CERTIFICATE"],
        "headline": head,
        "usability_checks": checks,
        "open_obligation": (
            "complete 18/21-state source-correlated word dissipation is not established;"
            " the first deployed accelerometer update still reports a nuisance correction"
            " above the sector-invariance budget at every candidate angle"
        ),
        "verdict": _verdict(checks, "complete-word dissipation"),
    }


def _p5_stage(d: dict) -> dict:
    p5 = d["P5_outer"]
    e = p5["declared_P5_entrance"]
    checks = [
        _check("outer capture certificate", p5["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"],
               "PASS", p5["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"] == "PASS"),
        _check("declared entrance full-attitude angle (deg)",
               e["gauged_full_attitude_angle_upper_deg"], MIN_P5_ENTRANCE_ANGLE_DEG,
               float(e["gauged_full_attitude_angle_upper_deg"]) >= MIN_P5_ENTRANCE_ANGLE_DEG),
        _check("all source handoff branches enter the outer sector",
               p5["all_source_handoff_branches_enter_outer_sector"], True,
               p5["all_source_handoff_branches_enter_outer_sector"] is True),
        _check("outer words needed from handoff", p5["N_outer_words"], 0,
               p5["N_outer_words"] == 0),
        _check("P1 handoff box not replaced",
               p5["P1_conservative_handoff_box_replaced"], False,
               p5["P1_conservative_handoff_box_replaced"] is False),
        _check("ungauged inclusion uses the upper cosine enclosure",
               p5["upper_cosine_enclosure_used_for_ungauged_inclusion"], True,
               p5["upper_cosine_enclosure_used_for_ungauged_inclusion"] is True),
        _check("no full heading radius assigned to the ungauged branch",
               p5["branches"]["timeout_ungauged"]["full_heading_radius_assigned"], False,
               p5["branches"]["timeout_ungauged"]["full_heading_radius_assigned"] is False),
        _check("legacy microscopic inner seed not used as capture target",
               p5["legacy_microscopic_inner_seed_used_as_outer_capture_target"], False,
               p5["legacy_microscopic_inner_seed_used_as_outer_capture_target"] is False),
    ]
    return {
        "stage": "P5",
        "title": "startup entrance into the outer sector and finite inner capture",
        "status": p5["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"],
        "headline": {
            "declared_entrance_full_attitude_angle_deg":
                e["gauged_full_attitude_angle_upper_deg"],
            "declared_entrance_cayley_norm_upper":
                e["attitude_geometry"]["cayley_norm_upper"],
            "entrance_position_component_abs_error_upper_Hs_factor":
                e["position_component_abs_error_upper_Hs_factor"],
            "entrance_position_norm_upper_Hs_factor":
                e["position_norm_upper_Hs_factor"],
            "outer_sector_angle_deg": p5["outer_sector_angle_deg"],
            "N_outer_words": p5["N_outer_words"],
            "N_inner_words": None,
            "branch_handoff_cayley_norms": {
                b: p5["branches"][b].get("P1_handoff_cayley_norm_upper")
                for b in p5["branches"]
            },
        },
        "usability_checks": checks,
        "open_obligation": (
            "finite capture from the outer sector to the inner stochastic localization"
            " level is downstream of the open P4 complete-word dissipation; N_inner"
            " remains unset"
        ),
        "verdict": _verdict(checks, "finite inner capture"),
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, translation: dict | None = None,
          first_accel: dict | None = None, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()

    startup = P1.build(path)
    startup_failures = P1.validate(startup)
    sector = SECTOR.build(path)
    sector_failures = SECTOR.validate(sector)
    p2 = P2.build(path)
    p2_failures = P2.validate(p2)
    word = WORD.build(path)
    word_failures = WORD.validate(word)
    p5 = P5.build(path)
    p5_failures = P5.validate(p5)
    ceiling = CEILING.build(path)
    ceiling_failures = CEILING.validate(ceiling)
    budget = BUDGET.build(path, source_pieces=source_pieces,
                          alignment_pieces=alignment_pieces,
                          force_magnitude_pieces=force_magnitude_pieces)
    budget_failures = BUDGET.validate(budget)

    overlap = sector["P1_overlap"]
    ctx = {
        "startup": startup,
        "startup_validation_pass": not startup_failures,
        "P1_normal_gauged_cayley_norm_upper": overlap["normal_gauged_cayley_norm_upper"],
        "P1_timeout_gauged_cayley_norm_upper": overlap["timeout_gauged_cayley_norm_upper"],
        "P2": p2,
        "word": word,
        "P4_sector": sector,
        "P4_first_accel_budget": budget,
        "P5_outer": p5,
        "route_ceiling": ceiling,
        "translation": translation,
        "first_accel": first_accel,
    }

    stages = [_p1_stage(ctx), _p2_stage(ctx), _p3_stage(ctx),
              _p4_stage(ctx), _p5_stage(ctx)]

    failures: list[str] = []
    failures += [f"P1 startup: {x}" for x in startup_failures]
    failures += [f"P2 source path: {x}" for x in p2_failures]
    failures += [f"P3/P4 word: {x}" for x in word_failures]
    failures += [f"P4 sector: {x}" for x in sector_failures]
    failures += [f"P4 first-accel budget: {x}" for x in budget_failures]
    failures += [f"P5 outer: {x}" for x in p5_failures]
    failures += [f"route ceiling: {x}" for x in ceiling_failures]

    not_usable = [s["stage"] for s in stages if s["verdict"] == "NOT_USABLE"]
    open_obligations = {s["stage"]: s["open_obligation"]
                        for s in stages if s["open_obligation"]}

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P1_P5_CERTIFICATE_NUMBERS_AND_USABILITY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "upstream_pass_fields_trusted_as_usability": False,
        "minimum_usability_contracts": {
            "normal_handoff_cayley_norm_lower": MIN_NORMAL_HANDOFF_CAYLEY,
            "timeout_handoff_cayley_norm_lower": MIN_TIMEOUT_HANDOFF_CAYLEY,
            "handoff_position_error_lower_m": MIN_HANDOFF_POSITION_ERROR_M,
            "handoff_velocity_error_lower_mps": MIN_HANDOFF_VELOCITY_ERROR_MPS,
            "handoff_latent_acceleration_error_lower_mps2": MIN_HANDOFF_AW_ERROR_MPS2,
            "sector_angle_lower_rad": MIN_SECTOR_ANGLE_RAD,
            "sector_strong_monotonicity_lower": MIN_SECTOR_MONOTONICITY,
            "sector_eta_information_ratio_upper": MAX_SECTOR_ETA_RATIO,
            "P5_entrance_angle_lower_deg": MIN_P5_ENTRANCE_ANGLE_DEG,
            "microscopic_cayley_threshold": MICROSCOPIC_CAYLEY_THRESHOLD,
        },
        "stages": stages,
        "stages_not_usable": not_usable,
        "open_theorem_obligations": open_obligations,
        "all_stage_geometries_usable": not not_usable,
        "P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED": not not_usable and not open_obligations,
        "translation_block_supplied": translation is not None,
        "signed_first_accelerometer_supplied": first_accel is not None,
        "next_obligation": (
            "close the complete 18/21-state P4 word dissipation with an operation-matched"
            " information decrease and a directional block margin, then derive the finite"
            " P5 inner capture word count from the declared 45 deg entrance"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit",):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "upstream_pass_fields_trusted_as_usability"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    stages = d.get("stages", [])
    if [s["stage"] for s in stages] != ["P1", "P2", "P3", "P4", "P5"]:
        f.append("stage list is not P1..P5")
    for s in stages:
        if s["verdict"] not in ("USABLE", "USABLE_GEOMETRY_OPEN_OBLIGATION", "NOT_USABLE"):
            f.append(f"{s['stage']}: unknown verdict")
        for c in s.get("usability_checks", []):
            if not c["pass"]:
                f.append(f"{s['stage']}: usability check failed: {c['check']}")
    if d.get("P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED") is True and d.get(
            "open_theorem_obligations"):
        f.append("complete proof claimed with open obligations")
    return list(dict.fromkeys(f))


def _load(p: Path | None) -> dict | None:
    if p is None:
        return None
    return json.loads(Path(p).read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--translation", type=Path, default=None,
                    help="pre-computed ou3_p4_translation_full_word_rigorous output")
    ap.add_argument("--first-accel", type=Path, default=None,
                    help="pre-computed ou3_p4_30deg_signed_first_accel_sector_v3 output")
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), translation=_load(a.translation),
              first_accel=_load(a.first_accel), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "stages": [{
            "stage": s["stage"],
            "status": s["status"],
            "verdict": s["verdict"],
            "open_obligation": s["open_obligation"],
        } for s in d["stages"]],
        "all_stage_geometries_usable": d["all_stage_geometries_usable"],
        "P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED":
            d["P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
