#!/usr/bin/env python3
"""Maximal chart-retaining P4 basin frontier, as a trade-off surface.

The design goal is NOT the smallest entry set that proves easily.  Every bit of
width certified here is width the Mahony/proxy capture stage does not have to
achieve, so the objective is to make the P4 basin as large as the certificate
allows and to expose the trade-off rather than pick one convenient number.

The retained every-prefix experiment ``entry-block-retention.json`` reports, per
mode, the attitude reach of each declared entry ball taken ALONE and the common
template reach.  Scaling entry ball ``j`` by ``s_j >= 0`` scales exactly that
term of the per-prefix subadditive sum, so at every completed prefix

    total_attitude(s) <= sum_j s_j a_j + t,

with ``a_j`` the single-ball attitude reach at the declared radius and ``t`` the
template reach, both already normalised by the declared attitude radius.  The
chart hypothesis of ``thm:brmm-bounded-bias-motion`` needs the attitude Cayley
excursion to stay inside the declared chart, so the certified basin is the
half-space

    sum_j s_j a_j <= B,    B := chart_cap/R_attitude - t,

a linear trade-off surface in scaled-radius space.  Any point of it is an
admissible P4 basin as far as the chart constraint is concerned, and the basin
is a correlated polytope, never a product box.

Reported frontier points, all maxima rather than convenient choices:

  * the largest UNIFORM inflation of the declared box;
  * the maximum-VOLUME point, which is the barycentric optimum of the simplex;
  * the largest radius of each single coordinate with all others left at their
    declared radius, which for the integral coordinate additionally uses the
    sharper retained joint maximum instead of the subadditive sum;
  * the attitude budget as an explicit function of the other radii, because
    attitude is the coordinate the Mahony/proxy stage has to deliver and the
    only one carrying a chart constraint.

This is a frozen-map binary64 diagnostic on one captured word.  It defines the
CANDIDATE maximal basin that the outward nonlinear certificate must then close;
it is not itself that certificate and promotes nothing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RETENTION_EVIDENCE = REPO / "reports" / "results" / "rao_stability" / "entry-block-retention.json"
QUALIFICATION = "OU3_P4_CHART_RETAINING_BASIN_FRONTIER_V1"


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _mode_frontier(mode_row: dict, radii: dict, entries: list[str]) -> dict:
    a = {j: float(mode_row["single_ball_reach"]["attitude"][j]) for j in entries}
    t = float(mode_row["common_template_reach"]["attitude"])
    r_att = float(radii["attitude"])
    cap = float(mode_row["declared_chart_cayley_norm_upper"])
    budget = down(down(cap / up(r_att)) - t)
    total_a = up(sum(a.values()))
    if budget <= 0.0 or total_a <= 0.0:
        raise RuntimeError("degenerate chart budget or reach row")

    uniform = down(budget / total_a)
    n = len(entries)
    max_volume = {j: down(down(budget / n) / up(a[j])) for j in entries}

    per_axis = {}
    for j in entries:
        rest = up(total_a - a[j])
        head = down(budget - rest)
        per_axis[j] = {
            "others_at_declared_radius": True,
            "max_scale": down(head / up(a[j])) if head > 0.0 else 0.0,
            "feasible_with_others_declared": bool(head > 0.0),
        }
    # The integral coordinate has a sharper retained bound: the joint maximum
    # with its independent ball removed is recorded, so only its own term needs
    # the subadditive relaxation.
    without_s = float(mode_row["without_independent_integral_ball"]["attitude"])
    head_s = down(down(cap / up(r_att)) - without_s)
    per_axis["integral_displacement"]["sharper_max_scale_from_retained_joint_maximum"] = (
        down(head_s / up(a["integral_displacement"])) if head_s > 0.0 else 0.0)
    per_axis["integral_displacement"]["sharper_bound_used"] = True

    declared_total = up(total_a + t)
    return {
        "attitude_entry_radius": r_att,
        "declared_chart_cayley_norm_upper": cap,
        "single_ball_attitude_reach": a,
        "template_attitude_reach": t,
        "half_space_coefficients": a,
        "half_space_budget": budget,
        "declared_box_total_attitude_reach": declared_total,
        "declared_box_retains_chart": bool(up(declared_total * r_att) <= cap),
        "declared_box_overshoot_factor": up(up(declared_total * r_att) / down(cap)),
        "max_uniform_scale": uniform,
        "max_uniform_radii": {j: down(uniform * float(radii[j])) for j in entries},
        "max_volume_scale": max_volume,
        "max_volume_radii": {j: down(max_volume[j] * float(radii[j])) for j in entries},
        "per_axis_max_scale_with_others_declared": per_axis,
        "attitude_budget_functional": (
            "s_attitude <= (B - sum_{j != attitude} s_j a_j) / a_attitude, "
            "with a and B reported above"),
    }


def build(retention_path: Path = RETENTION_EVIDENCE) -> dict:
    payload = json.loads(Path(retention_path).read_text(encoding="utf-8"))
    radii = payload["declared_entry_radii"]
    modes = {}
    for mode, row in payload["modes"].items():
        entries = list(row["entry_order"])
        if set(entries) != set(row["single_ball_reach"]["attitude"]):
            raise RuntimeError("retained reach row does not cover the declared entry order")
        modes[mode] = _mode_frontier(row, radii, entries)

    binding = min(modes, key=lambda m: modes[m]["max_uniform_scale"])
    integral_binding = min(
        modes,
        key=lambda m: modes[m]["per_axis_max_scale_with_others_declared"]
        ["integral_displacement"]["sharper_max_scale_from_retained_joint_maximum"])

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "objective": "maximise the certified P4 basin, never shrink it for proof convenience",
        "basin_is_a_correlated_polytope_not_a_product_box": True,
        "declared_entry_radii": radii,
        "modes": modes,
        "binding_mode_for_uniform_scale": binding,
        "max_uniform_scale_over_modes": modes[binding]["max_uniform_scale"],
        "binding_mode_for_integral_axis": integral_binding,
        "max_integral_radius_with_others_declared_m_s": down(
            modes[integral_binding]["per_axis_max_scale_with_others_declared"]
            ["integral_displacement"]["sharper_max_scale_from_retained_joint_maximum"]
            * float(radii["integral_displacement"])),
        "constraint_modelled_here": "attitude Cayley chart retention only",
        "other_physical_caps_not_applied_here": {
            "accelerometer_bias": ("the deployed radial projection keeps |e_b| <= R + B_true, "
                                   "so a frontier point above that is chart-admissible but not "
                                   "physically attainable"),
            "attitude": "the chart bound is the constraint modelled; no separate cap is added",
            "note": ("a frontier point is a candidate basin for the CHART hypothesis; every "
                     "coordinate still has to satisfy its own deployed physical cap"),
        },
        "frontier_points_are_chart_admissible_not_physically_screened": True,
        "evidence_sha256": hashlib.sha256(Path(retention_path).read_bytes()).hexdigest(),
        "frozen_map_binary64_diagnostic_not_outward_certificate": True,
        "defines_candidate_basin_to_be_certified": True,
        "entry_radii_reduced_for_proof_convenience": False,
        "outward_nonlinear_certificate_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "close the outward nonlinear certificate on the largest frontier point the "
            "deployed Live-entry reachability actually needs, then hand that basin to the "
            "Mahony/proxy finite-time capture proof"),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("source changed")
    if d.get("constraint_modelled_here") != "attitude Cayley chart retention only":
        f.append("frontier scope changed without rebuilding it")
    if not d.get("other_physical_caps_not_applied_here", {}).get("accelerometer_bias"):
        f.append("bias projection cap caveat dropped")
    for k in ("basin_is_a_correlated_polytope_not_a_product_box",
              "frozen_map_binary64_diagnostic_not_outward_certificate",
              "defines_candidate_basin_to_be_certified",
              "frontier_points_are_chart_admissible_not_physically_screened"):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in ("entry_radii_reduced_for_proof_convenience",
              "outward_nonlinear_certificate_closed_here", "P4_promoted_here"):
        if d.get(k) is not False:
            f.append(k + " not false")
    if not d.get("modes"):
        f.append("no mode frontier produced")
    for mode, row in d.get("modes", {}).items():
        budget = float(row.get("half_space_budget", -1.0))
        if not (math.isfinite(budget) and budget > 0.0):
            f.append(mode + " chart budget is not positive")
        coeff = row.get("half_space_coefficients", {})
        if not coeff or any(float(v) < 0.0 for v in coeff.values()):
            f.append(mode + " half-space coefficients invalid")
        uniform = float(row.get("max_uniform_scale", -1.0))
        if not (0.0 < uniform):
            f.append(mode + " uniform scale is not positive")
        # The declared box is the reference the frontier is measured against;
        # if it ever starts retaining the chart the frontier must be rebuilt.
        if row.get("declared_box_retains_chart") is not False and mode == "H18":
            f.append("H18 declared box now retains the chart; rebuild the frontier")
        # A maximum-volume point must saturate the same budget.
        mv = row.get("max_volume_scale", {})
        if set(mv) != set(coeff):
            f.append(mode + " max-volume point does not cover every entry")
        spend = sum(float(mv[j]) * float(coeff[j]) for j in coeff)
        if not (spend <= up(budget * 1.000001)):
            f.append(mode + " max-volume point exceeds its own budget")
        per_axis = row.get("per_axis_max_scale_with_others_declared", {})
        integral = per_axis.get("integral_displacement", {})
        if integral.get("sharper_bound_used") is not True:
            f.append(mode + " integral axis lost its sharper retained bound")
        if float(integral.get("sharper_max_scale_from_retained_joint_maximum", -1.0)) <= 0.0:
            f.append(mode + " integral axis has no admissible radius")
        # The sharper retained bound must never be weaker than the subadditive one.
        if float(integral.get("sharper_max_scale_from_retained_joint_maximum", 0.0)) < \
                float(integral.get("max_scale", 0.0)):
            f.append(mode + " sharper integral bound is weaker than the subadditive one")
    x = float(d.get("max_integral_radius_with_others_declared_m_s", -1.0))
    if not (math.isfinite(x) and x > 0.0):
        f.append("integral frontier radius is not a positive finite bound")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--retention", type=Path, default=RETENTION_EVIDENCE)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.retention)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "max_uniform_scale": d["max_uniform_scale_over_modes"],
        "binding_mode": d["binding_mode_for_uniform_scale"],
        "max_integral_radius_m_s": d["max_integral_radius_with_others_declared_m_s"],
        "H18_max_volume_radii": d["modes"]["H18"]["max_volume_radii"],
        "H18_declared_overshoot": d["modes"]["H18"]["declared_box_overshoot_factor"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
