"""Per-entry-block every-prefix retention of the declared coordinate domain.

`domain-retention.json` already records two facts about the declared initial
set: the closed .4 m/s^2 bias ball plus the common template is retained with a
large margin, and the complete declared product box is not, limited by velocity
in both modes.  It reports the box result as one subadditive total per
coordinate, so the total does not say which declared ball causes it.

This experiment switches each declared ball on ALONE and maximizes over every
completed prefix of the same passively attached capture.  Entry `j` alone
reaches `||Pi_g Phi_ell E_j||_2 R_j / R_g` in coordinate `g`, which is exactly
the term entry `j` contributes to the subadditive total, so the rows sum back
to it and no term is hidden inside another.

Two consequences follow that the total cannot show.

First, the accelerometer-bias ball is not the limiter of ANY coordinate, not
merely of the total: its largest single-ball reach over both modes and all six
motion coordinates is well under one radius.  The separated bias budget was
retired as the chart limiter; this retires the bias ENTRY ball as a limiter of
every individual coordinate too.

Second, the reach is compared against the declared Cayley chart bound rather
than only against the entry radii, because that is what the theorem's level
hypothesis actually constrains.  With the full declared box the H18 attitude
coordinate leaves the chart the frozen map is expanded in, so no chart-valid
L_chart exists and the hypothesis is unsatisfiable at that entry set; with the
independent integral-displacement ball removed it stays inside and a
chart-valid L_chart does exist.

The minimal working radii reported here are the smallest radii this word
retains at every completed prefix.  They are a candidate chart-valid level
L_chart for `thm:brmm-bounded-bias-motion`, whose retention hypothesis is
`Gamma*L + C_p < L_chart` with `L` the entry level: an excursion above an entry
radius is not a failure of that hypothesis, and only the attitude row carries a
chart constraint.  They enlarge the retention target and never reduce the entry
set, so every downstream nonlinear obligation becomes harder, not easier.

Binary64 point arithmetic on one captured word with frozen coefficients.  It is
a diagnostic, not an outward certificate, and promotes nothing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_bounded_bias_cocycle as C
import ou3_p4_motion_gain as G

GROUPS = ("attitude", "gyro_bias", "velocity", "position",
          "integral_displacement", "latent_acceleration")
ENTRIES = GROUPS + ("accelerometer_bias",)
# ou3_p4_closure_domain.json: uniform_nonlinear_enclosure.required_cayley_state_norm_upper
DECLARED_CHART_CAYLEY_NORM_UPPER = 1.0


def block_retention(samples, radii):
    """Worst completed-prefix reach of each declared ball taken alone.

    Returns per coordinate: the single-ball reaches, the common-template reach,
    the subadditive total, and the same total without the independent
    integral-displacement ball.
    """
    n = len(ENTRIES)
    alone = np.zeros((len(GROUPS), n))
    template = np.zeros(len(GROUPS))
    total = np.zeros(len(GROUPS))
    total_without_integral = np.zeros(len(GROUPS))
    integral = ENTRIES.index("integral_displacement")
    transition, response = np.eye(21), np.zeros(21)
    prefixes = 0
    for sample in samples:
        for full, forcing, _ in sample["prefixes"]:
            phi, r = full @ transition, full @ response + forcing
            if not np.isfinite(phi).all() or not np.isfinite(r).all():
                raise ValueError("nonfinite prefix propagation")
            prefixes += 1
            for g in range(len(GROUPS)):
                sl = slice(3*g, 3*g+3)
                reach = [float(np.linalg.norm(phi[sl, 3*j:3*j+3], 2))*radii[j]
                         for j in range(n)]
                t = float(np.linalg.norm(r[sl]))
                alone[g] = np.maximum(alone[g], np.array(reach)/radii[g])
                template[g] = max(template[g], t/radii[g])
                total[g] = max(total[g], (sum(reach)+t)/radii[g])
                total_without_integral[g] = max(
                    total_without_integral[g], (sum(reach)-reach[integral]+t)/radii[g])
        transition, response = full @ transition, full @ response + forcing
    return {"completed_prefix_count": prefixes, "single_ball_reach": alone,
            "template_reach": template, "subadditive_total": total,
            "subadditive_total_without_independent_integral_ball": total_without_integral}


def summarize(result, radii):
    alone = result["single_ball_reach"]
    total = result["subadditive_total"]
    without = result["subadditive_total_without_independent_integral_ball"]
    bias = ENTRIES.index("accelerometer_bias")
    attitude = GROUPS.index("attitude")
    chart_full = float(total[attitude]*radii[attitude])
    chart_without = float(without[attitude]*radii[attitude])
    return {
      "completed_prefix_count": result["completed_prefix_count"],
      "entry_order": list(ENTRIES),
      "single_ball_reach": {g: {e: float(alone[i][j]) for j, e in enumerate(ENTRIES)}
                            for i, g in enumerate(GROUPS)},
      "common_template_reach": {g: float(result["template_reach"][i]) for i, g in enumerate(GROUPS)},
      "declared_box_subadditive_total": {g: float(total[i]) for i, g in enumerate(GROUPS)},
      "without_independent_integral_ball": {g: float(without[i]) for i, g in enumerate(GROUPS)},
      # A minimal working radius is the entry radius times the reach it must absorb.
      "minimal_working_radius_inflation_declared_box": {g: float(max(1.0, total[i])) for i, g in enumerate(GROUPS)},
      "minimal_working_radius_inflation_without_integral_ball": {g: float(max(1.0, without[i])) for i, g in enumerate(GROUPS)},
      "accelerometer_bias_ball_worst_single_reach": float(np.max(alone[:, bias])),
      "accelerometer_bias_ball_limits_no_coordinate": bool(np.max(alone[:, bias]) < 1.0),
      "worst_single_ball_per_coordinate": {g: ENTRIES[int(np.argmax(alone[i]))] for i, g in enumerate(GROUPS)},
      "declared_chart_cayley_norm_upper": DECLARED_CHART_CAYLEY_NORM_UPPER,
      "attitude_cayley_norm_reached_declared_box": chart_full,
      "attitude_cayley_norm_reached_without_integral_ball": chart_without,
      "declared_box_stays_inside_declared_chart": bool(chart_full <= DECLARED_CHART_CAYLEY_NORM_UPPER),
      "without_integral_ball_stays_inside_declared_chart": bool(chart_without <= DECLARED_CHART_CAYLEY_NORM_UPPER),
      "entry_radii_reduced": False,
      "working_domain_is_larger_than_entry_set": True,
      "frozen_point_binary64_diagnostic_not_outward_certificate": True,
      "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--attachment", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    attachment = json.loads(args.attachment.read_text())
    paths = {s: Path(str(args.prefix)+s) for s in
             (".root.json", ".inputs.csv", ".prefixes.jsonl", ".events.jsonl")}
    hashes = {s: hashlib.sha256(p.read_bytes()).hexdigest() for s, p in paths.items()}
    if not attachment["read_only_trace_recovers_baseline_bit_for_bit"]:
        raise ValueError("passive trace parity required")
    if hashes[".events.jsonl"] != attachment["event_trace_sha256"] or any(
            hashes[s] != attachment["baseline_capture_sha256"][s] for s in
            (".root.json", ".inputs.csv", ".prefixes.jsonl")):
        raise ValueError("detached source/word capture")

    root = json.loads(paths[".root.json"].read_text())
    rows = [json.loads(line) for line in paths[".events.jsonl"].read_text().splitlines()]
    points = [json.loads(line) for line in paths[".prefixes.jsonl"].read_text().splitlines()]
    domain_path = Path(__file__).resolve().parents[2] / "tools/stability/ou3_proof_operating_domain.json"
    domain = json.loads(domain_path.read_text())["startup"]["physical_handoff_coordinate_bounds"]
    keys = ("gyro_bias_error_norm_upper_rad_s", "velocity_error_norm_upper_mps",
            "position_error_norm_upper_m", "integral_displacement_error_norm_upper_m_s",
            "latent_acceleration_error_norm_upper_mps2", "accelerometer_bias_error_norm_upper_mps2")
    radii = np.array([2*np.tan(np.pi/12)] + [float(domain[k]) for k in keys])
    if radii.shape != (len(ENTRIES),) or np.any(radii <= 0):
        raise ValueError("positive declared radius required for every entry block")

    report = {"experiment": "DECLARED_ENTRY_BLOCK_EVERY_PREFIX_RETENTION",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "declared_entry_radii": dict(zip(ENTRIES, radii.tolist())),
              "entry_radii_reduced": False, "filter_changed": False,
              "source_or_coefficient_search_used": False,
              "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False, "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached mode: "+mode)
        steps, _, _, defects, counts = G.build_word(root, rows, points, mode)
        if defects.failures:
            raise ValueError("finite factorization failed: "+repr(defects.failures[:1]))
        samples = C.sample_lifts(steps)
        summary = summarize(block_retention(samples, radii), radii)
        summary["sample_count"] = len(samples)
        summary["counts"] = dict(counts)
        report["modes"][mode] = summary
        print("ENTRY_BLOCK_RETENTION", mode,
              "bias_ball_worst=%.6f" % summary["accelerometer_bias_ball_worst_single_reach"],
              "chart_box=%.6f" % summary["attitude_cayley_norm_reached_declared_box"],
              "chart_no_integral=%.6f" % summary["attitude_cayley_norm_reached_without_integral_ball"],
              flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
