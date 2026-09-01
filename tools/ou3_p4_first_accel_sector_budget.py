#!/usr/bin/env python3
"""First-accelerometer sector-invariance budget across the P4 candidate ladder.

The P4 candidate ladder in ``docs/ou-iii-p5-sea-scaled-entrance.md`` is ordered
``30 -> 25 -> 20 -> 15`` deg on the assumption that a narrower candidate makes
the complete-word dissipation easier.  That assumption is correct for the
exact Cayley monotonicity and eta ratios, which do improve as the angle drops.
It is not correct for the operation that actually blocks the certificate.

This producer measures, on the same source/alignment/force children the signed
first-accelerometer stage uses, the two quantities that decide whether a single
deployed accelerometer update can keep the state inside the operation-matched
outer sector:

* the **budget** -- the largest correction norm whose worst-case signed Cayley
  composition with the post-prediction candidate norm still lands strictly
  inside the outer sector.  The budget grows as the candidate angle shrinks;
* the **nuisance term** -- the part of the correction driven by the effective
  ``a_w`` input: the declared 0.3 g startup latent-acceleration error, the
  accelerometer bias, and the finite-angle force remainder, multiplied by the
  shared-force-magnitude accelerometer gain.

The nuisance term is dominated by a contribution that does *not* shrink with
the candidate angle, because it is set by the ratio of the declared latent
acceleration error to the lowest admitted specific-force magnitude.  Both the
descending and the ascending ladder are therefore measured here and reported
side by side.

This is a **distance, not a verdict**.  The nuisance term is an outward bound,
so a large value does not prove that any admissible state actually leaves the
sector; the ``a_w`` covariance endpoints used by the gain are enclosure
endpoints, not certified reachable points.  Nothing here promotes or retires a
theorem claim, changes the filter, or shrinks the declared 30 deg / 0.3 g
entrance.  It exists to say which of the three structural changes named by
``ou3_p4_p5_route_ceiling_certificate`` each ladder rung still needs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_30deg_signed_first_accel_sector as V1
import ou3_p4_candidate_aw_capture_budget as AWB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p4_shared_force_gain as SHARED
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_signed_cayley_cell as SIGNED
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DESCENDING_LADDER_DEG = (30.0, 25.0, 20.0, 15.0)
ASCENDING_LADDER_DEG = (35.0, 40.0, 45.0)
# Not candidates: probes that show whether the descending ladder can ever
# bring the nuisance term inside the budget by shrinking the angle alone.
LIMIT_PROBE_DEG = (10.0, 5.0, 1.0)
BUDGET_SEARCH_UPPER_RAD = V1.CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD
BUDGET_BISECTION_STEPS = 60


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _worst_case_qplus(q: float, dnorm: float) -> float:
    """Post-update Cayley upper for a correction of norm ``dnorm`` at worst sign.

    The adverse alignment is ``d`` parallel to ``c``; the signed composition
    helper of the shipping stage is reused unchanged so the budget is measured
    in exactly the arithmetic the sector certificate uses.
    """
    scale = SIGNED.correction_cayley_scale_interval(dnorm)
    adot = up(up(scale.hi * dnorm) * q)
    return float(V1._compose_q_upper(q, dnorm, adot)["post_update_q_upper"])


def _sector_correction_budget(q: float, outer_q: float) -> dict:
    """Smallest correction norm already outside the outer sector (upper budget).

    ``_worst_case_qplus`` is increasing in ``dnorm``, so a bisection on the
    crossing point is well defined.  The reported budget is the *failing* side
    of the bracket, hence an upper bound on any admissible correction norm.
    """
    if q >= outer_q:
        return {
            "sector_invariance_correction_budget_upper_rad": 0.0,
            "candidate_already_outside_outer_sector": True,
            "budget_bracket_rad": [0.0, 0.0],
        }
    lo = 0.0
    hi = BUDGET_SEARCH_UPPER_RAD
    try:
        if _worst_case_qplus(q, down(hi)) < outer_q:
            return {
                "sector_invariance_correction_budget_upper_rad": hi,
                "candidate_already_outside_outer_sector": False,
                "budget_bracket_rad": [hi, hi],
                "budget_exceeds_monotone_cayley_chart": True,
            }
    except RuntimeError:
        pass
    for _ in range(BUDGET_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        try:
            ok = _worst_case_qplus(q, mid) < outer_q
        except RuntimeError:
            ok = False
        if ok:
            lo = mid
        else:
            hi = mid
    return {
        "sector_invariance_correction_budget_upper_rad": up(hi),
        "candidate_already_outside_outer_sector": False,
        "budget_bracket_rad": [down(lo), up(hi)],
    }


def _gain_table(path: Path, domain: dict, *, source_pieces: int,
                alignment_pieces: int, force_magnitude_pieces: int) -> dict:
    """Angle-independent per-cell gain and effective-``a_w`` prefix data."""
    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = VECTOR.build()["configured_measurement_bounds"]
    racc_var = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))[0][0]
    startup = domain["startup"]
    handoff = startup["physical_handoff_coordinate_bounds"]
    aw0 = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    ba = float(handoff["accelerometer_bias_error_norm_upper_mps2"])
    pnorm = AWB._p5_position_norm_upper(domain)
    live = domain["normal_live"]
    force_lower = float(live["specific_force_norm_lower_mps2"])
    force_upper = float(live["specific_force_norm_upper_mps2"])
    force_cells = RG._geom_ranges(force_lower, force_upper, force_magnitude_pieces)
    xcells = SG._linear_cells(alignment_pieces)

    cells = []
    for si, (src, phase) in enumerate(RG._source_phase_children(source_pieces)):
        P0 = FULL._initial_covariance(src, path)
        Fm, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(
            FULL.matrix_mul(FULL.matrix_mul(Fm, P0), FULL.matrix_transpose(Fm)), Q))
        _s, _a, paw_pred = RG._scalar_axis_structure(Pp)
        paw = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)[0] if phase == "due" else paw_pred
        intercept, slope, _d = AWB._s_phase_affine_aw_bound(src, phase, Pp, domain, pnorm)
        aw = up(intercept + up(slope * aw0))
        for xi, x in enumerate(xcells):
            for mi, m in enumerate(force_cells):
                k, _kh, _gd = SHARED.shared_force_structured_gain_bounds(
                    tilt=tilt, yaw=yaw, eps=eps, x=x, m=m,
                    paw=paw, racc_var=racc_var)
                cells.append({
                    "source_phase_cell": si, "pseudo_phase": phase,
                    "alignment_cell": xi, "force_cell": mi,
                    "force_magnitude_mps2": m.as_list(),
                    "alignment_x": x.as_list(),
                    "aw_after_prefix_upper_mps2": aw,
                    "Ktheta_norm_upper": k,
                    "force_upper_mps2": m.hi,
                })
    return {
        "cells": cells,
        "accel_bias_error_norm_upper_mps2": ba,
        "declared_startup_aw_error_norm_upper_mps2": aw0,
        "specific_force_norm_lower_mps2": force_lower,
        "specific_force_norm_upper_mps2": force_upper,
        "prediction_step_s": h,
    }


def _ladder_row(deg: float, table: dict, domain: dict, outer_q: float,
                entrance_q: float, family: str) -> dict:
    geom = SECTOR._validated_design_geometry(math.radians(float(deg)))
    q0 = float(geom["cayley_norm_upper"])
    q = RG._q_after_first_prediction(q0, domain, table["prediction_step_s"])
    budget = _sector_correction_budget(q, outer_q)
    eta_per_vector = VEFF.accel_attitude_eta_per_vector_norm_upper(q)
    ba = table["accel_bias_error_norm_upper_mps2"]

    worst = None
    for c in table["cells"]:
        eta_force = up(eta_per_vector * c["force_upper_mps2"])
        rho = up(c["aw_after_prefix_upper_mps2"] + up(eta_force + ba))
        nuisance = up(c["Ktheta_norm_upper"] * rho)
        if worst is None or nuisance > worst["nuisance_correction_norm_upper_rad"]:
            worst = dict(c)
            worst.pop("force_upper_mps2", None)
            worst["force_attitude_remainder_upper_mps2"] = eta_force
            worst["effective_aw_input_upper_mps2"] = rho
            worst["nuisance_correction_norm_upper_rad"] = nuisance

    b = float(budget["sector_invariance_correction_budget_upper_rad"])
    n = float(worst["nuisance_correction_norm_upper_rad"])
    return {
        "angle_deg": float(deg),
        "ladder_family": family,
        "candidate_q_upper": q0,
        "post_prediction_q_upper": q,
        "operation_matched_outer_q_upper": outer_q,
        "inside_45deg_entrance": q0 < entrance_q,
        "eta_per_vector_norm_upper": eta_per_vector,
        **budget,
        "worst_cell": worst,
        "nuisance_correction_norm_upper_rad": n,
        "nuisance_over_budget_ratio": (math.inf if b <= 0.0 else up(n / b)),
        "nuisance_fits_inside_budget": bool(b > 0.0 and n < b),
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("sector budget domain must not be trajectory fitted")

    sector = SECTOR.build(path)
    entrance = ENTRANCE.build(path)
    failures = [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    outer_q = float(sector["design_cayley_norm_upper"])
    entrance_q = float(
        entrance["P5_entrance"]["attitude_geometry"]["cayley_norm_upper"])

    table = _gain_table(path, domain, source_pieces=source_pieces,
                        alignment_pieces=alignment_pieces,
                        force_magnitude_pieces=force_magnitude_pieces)

    rows = [_ladder_row(d, table, domain, outer_q, entrance_q, "descending")
            for d in DESCENDING_LADDER_DEG]
    rows += [_ladder_row(d, table, domain, outer_q, entrance_q, "ascending")
             for d in ASCENDING_LADDER_DEG]
    rows.sort(key=lambda r: r["angle_deg"])
    probes = [_ladder_row(d, table, domain, outer_q, entrance_q, "limit_probe")
              for d in LIMIT_PROBE_DEG]
    probes.sort(key=lambda r: r["angle_deg"])

    aw0 = table["declared_startup_aw_error_norm_upper_mps2"]
    flo = table["specific_force_norm_lower_mps2"]
    aw_over_force = up(aw0 / down(flo))
    direction_error = (
        math.asin(min(1.0, aw_over_force)) if aw_over_force <= 1.0 else None)

    any_fits = any(r["nuisance_fits_inside_budget"] for r in rows)
    descending_helps = False
    desc = [r for r in rows if r["ladder_family"] == "descending"]
    if len(desc) >= 2:
        by_angle = sorted(desc, key=lambda r: r["angle_deg"])
        descending_helps = (
            by_angle[0]["nuisance_over_budget_ratio"]
            < 0.5 * by_angle[-1]["nuisance_over_budget_ratio"])

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_FIRST_ACCELEROMETER_SECTOR_INVARIANCE_BUDGET",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_entrance_shrunk": False,
        "distance_only_no_verdict_emitted": True,
        "shared_force_magnitude_dependency_preserved": True,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "operation_matched_outer_q_upper": outer_q,
        "declared_P5_entrance_q_upper": entrance_q,
        "declared_startup_aw_error_norm_upper_mps2": aw0,
        "accel_bias_error_norm_upper_mps2": table["accel_bias_error_norm_upper_mps2"],
        "specific_force_norm_lower_mps2": flo,
        "specific_force_norm_upper_mps2": table["specific_force_norm_upper_mps2"],
        "declared_aw_over_lowest_specific_force": aw_over_force,
        "lowest_force_gravity_direction_error_upper_rad": direction_error,
        "audited_cells": len(table["cells"]),
        "ladder_rows": rows,
        "limit_probe_rows": probes,
        "limit_probes_are_not_candidate_angles": True,
        "smallest_probe_angle_deg": probes[0]["angle_deg"],
        "smallest_probe_nuisance_over_budget_ratio": probes[0]["nuisance_over_budget_ratio"],
        "shrinking_the_candidate_angle_alone_can_close_the_budget":
            any(r["nuisance_fits_inside_budget"] for r in probes),
        "any_ladder_rung_fits_nuisance_inside_budget": any_fits,
        "descending_ladder_halves_the_gap": descending_helps,
        "next_obligation": (
            "shrinking the candidate angle widens the budget but the nuisance term stays"
            " above it even at the smallest probe, because that term is set by the declared"
            " latent-acceleration error over the lowest admitted specific force; closing a"
            " rung needs the operation-matched information decrease and a directional block"
            " margin, not a narrower candidate"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit", "distance_only_no_verdict_emitted",
              "shared_force_magnitude_dependency_preserved"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "declared_entrance_shrunk",
              "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
              "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    rows = d.get("ladder_rows", [])
    if [r["angle_deg"] for r in rows] != sorted(r["angle_deg"] for r in rows):
        f.append("ladder rows are not ordered by angle")
    if not rows:
        f.append("no ladder rows emitted")
    for r in rows:
        if not (r["post_prediction_q_upper"] >= r["candidate_q_upper"]):
            f.append(f"{r['angle_deg']}: prediction did not widen the candidate norm")
        if r["sector_invariance_correction_budget_upper_rad"] < 0.0:
            f.append(f"{r['angle_deg']}: negative correction budget")
        if r["nuisance_correction_norm_upper_rad"] <= 0.0:
            f.append(f"{r['angle_deg']}: non-positive nuisance term")
    if d.get("audited_cells", 0) <= 0:
        f.append("no audited cells")
    if d.get("limit_probes_are_not_candidate_angles") is not True:
        f.append("limit probes are not marked as non-candidates")
    if not d.get("limit_probe_rows"):
        f.append("no limit probe rows emitted")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "outer_q": d["operation_matched_outer_q_upper"],
        "aw_over_lowest_force": d["declared_aw_over_lowest_specific_force"],
        "ladder": [{
            "deg": r["angle_deg"],
            "family": r["ladder_family"],
            "q_pred": r["post_prediction_q_upper"],
            "budget_rad": r["sector_invariance_correction_budget_upper_rad"],
            "nuisance_rad": r["nuisance_correction_norm_upper_rad"],
            "ratio": r["nuisance_over_budget_ratio"],
            "fits": r["nuisance_fits_inside_budget"],
        } for r in d["ladder_rows"]],
        "probes": [{
            "deg": r["angle_deg"],
            "budget_rad": r["sector_invariance_correction_budget_upper_rad"],
            "nuisance_rad": r["nuisance_correction_norm_upper_rad"],
            "ratio": r["nuisance_over_budget_ratio"],
            "fits": r["nuisance_fits_inside_budget"],
        } for r in d["limit_probe_rows"]],
        "any_fits": d["any_ladder_rung_fits_nuisance_inside_budget"],
        "angle_shrink_alone_closes": d["shrinking_the_candidate_angle_alone_can_close_the_budget"],
        "descending_halves_gap": d["descending_ladder_halves_the_gap"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
