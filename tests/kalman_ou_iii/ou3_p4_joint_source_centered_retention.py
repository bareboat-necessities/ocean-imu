"""Source-centered correlated-root retention on the 24D joint P4 point graph.

The preceding joint-storage diagnostic proved that a scalar Young/supply bound
throws away too much correlation.  This route keeps ONE initial correlated
error ellipsoid and ONE attached physical source/bias continuation.  No fresh
bias vector, source port, covariance, or metric is introduced at later samples.

For the frozen attached graph z_j = T_j z_0 + r_j, let zbar_j be the actual
same-history center and let the initial error deviation satisfy

    delta e_0 = L_0 y,  ||y||_2 <= sigma,

where L_0 is the Cholesky factor of the shipping root covariance.  The true
physical-bias root is fixed by BIAS1, so its deviation is zero in this point
experiment.  Every completed-event motion coordinate then obeys exactly

    ||Pi_G z_j|| <= ||Pi_G zbar_j||
                    + sigma ||Pi_G T_j [L_0;0]||_2.

This is a single correlated root-to-prefix map: active projection factors,
motion-to-bias feedback, actual R_S and physical source forcing were already
composed in the same 24D graph.  The same compatible cyclic M_i sequence is
used to report the motion-storage image of that root ellipsoid at every sample
boundary; no isolated per-word spectral metric is selected.

The covariance ellipsoid is still the FILTER'S believed uncertainty, not a
qualified hard true-error set.  Therefore a favorable point margin identifies
a quantitative covariance-consistency/entry-set target but cannot promote P4.
Nonzero physical-bias driver coverage and nonlinear coefficient/projection
outward enclosure also remain open.

For the integral-displacement group this diagnostic additionally reports the
minimum literal-prefix Euclidean headroom remaining after the one-sigma root
image.  That number is only a *candidate scale* for the new origin-invariant
centered-S recurrence D_S.  It is not a D_S certificate: the rigorous source
cover must propagate the same physical p/S_L history through every event and
cannot inject an independent additive S port at each prefix.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_joint_bias_source_storage as J
import ou3_p4_motion_gain as G
import ou3_p4_source_endpoint as SOURCE


def root_covariance(points, mode):
    point = next((p for p in points if p["word"] == mode and p.get("event") == "root"), None)
    if point is None:
        raise ValueError("missing attached word root")
    _, covariance = SOURCE.error_and_covariance(point, 21)
    root = np.linalg.cholesky(covariance)
    lift = np.zeros((24, 21))
    lift[:21, :] = root
    return covariance, lift


def audit_mode(root, rows, points, mode, radii):
    steps21, _, initial21, defects, counts = G.build_word(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite coefficient attachment failed: "+repr(defects.failures[:1]))
    steps24, augmentation = J.augmented_steps(root, rows, steps21, mode, float(root["dt"]))
    samples = J.sample_lifts(steps24)
    initial = J.initial_joint_state(root, initial21)
    states = J.entrance_states(samples, initial)
    metrics, _, metric_summary = J.periodic_metrics(samples, radii)
    covariance, root_lift = root_covariance(points, mode)

    names = ("attitude", "gyro_bias", "velocity", "position",
             "integral_displacement", "latent_acceleration")
    global_transition = np.eye(24)
    prefixes = []
    boundary_storage = []
    # Root boundary before the first sample.
    motion_root = root_lift[:18]
    mroot = np.linalg.cholesky(metrics[0]).T @ motion_root
    boundary_storage.append({"boundary": 0,
                             "one_sigma_storage_level": float(np.linalg.norm(mroot, 2)**2)})

    for i, sample in enumerate(samples):
        for full, forcing, step in sample["prefixes"]:
            transition = full @ global_transition
            center = full @ states[i] + forcing
            groups = {}
            for group, name in enumerate(names):
                sl = slice(3*group, 3*group+3)
                gain = float(np.linalg.norm(transition[sl] @ root_lift, 2))
                center_norm = float(np.linalg.norm(center[sl]))
                radius = float(radii[3*group])
                headroom = radius-center_norm
                one_sigma_total = center_norm+gain
                one_sigma_extra_headroom = radius-one_sigma_total
                critical = (headroom/gain if gain > 0 and headroom > 0
                            else (float("inf") if gain == 0 and headroom > 0 else 0.))
                groups[name] = {"center_norm": center_norm,
                                "one_sigma_deviation_gain": gain,
                                "one_sigma_retention_ratio": one_sigma_total/radius,
                                "one_sigma_extra_euclidean_headroom": one_sigma_extra_headroom,
                                "critical_initial_sigma_level": critical,
                                "declared_radius": radius}
            prefixes.append({"sample": i, "index": step["index"],
                             "stage": step["stage"], "kind": step["kind"],
                             "groups": groups})
        global_transition = sample["full"] @ global_transition
        boundary_motion = global_transition[:18] @ root_lift
        weighted = np.linalg.cholesky(metrics[i+1]).T @ boundary_motion
        boundary_storage.append({"boundary": i+1,
                                 "one_sigma_storage_level": float(np.linalg.norm(weighted, 2)**2)})

    summary = {}
    all_critical = []
    for name in names:
        values = [(p["groups"][name]["one_sigma_retention_ratio"], j, p)
                  for j, p in enumerate(prefixes)]
        ratio, ordinal, worst = max(values, key=lambda x: x[0])
        crit_values = [(p["groups"][name]["critical_initial_sigma_level"], j, p)
                       for j, p in enumerate(prefixes)]
        critical, cordinal, limiting = min(crit_values, key=lambda x: x[0])
        headroom_values = [(p["groups"][name]["one_sigma_extra_euclidean_headroom"], j, p)
                           for j, p in enumerate(prefixes)]
        extra_headroom, hord, hlim = min(headroom_values, key=lambda x: x[0])
        all_critical.append((critical, name, cordinal, limiting))
        summary[name] = {"maximum_one_sigma_retention_ratio": ratio,
                         "worst_prefix_ordinal": ordinal,
                         "worst_index": worst["index"],
                         "worst_stage": worst["stage"],
                         "worst_kind": worst["kind"],
                         "critical_initial_sigma_level": critical,
                         "critical_prefix_ordinal": cordinal,
                         "critical_index": limiting["index"],
                         "critical_stage": limiting["stage"],
                         "critical_kind": limiting["kind"],
                         "minimum_one_sigma_extra_euclidean_headroom": extra_headroom,
                         "headroom_prefix_ordinal": hord,
                         "headroom_index": hlim["index"],
                         "headroom_stage": hlim["stage"],
                         "headroom_kind": hlim["kind"]}
    critical, limiting_name, ordinal, limiting = min(all_critical, key=lambda x: x[0])
    max_storage = max(x["one_sigma_storage_level"] for x in boundary_storage)
    endpoint_storage = boundary_storage[-1]["one_sigma_storage_level"]
    one_sigma_retained = all(v["maximum_one_sigma_retention_ratio"] <= 1. for v in summary.values())
    sdiag = summary["integral_displacement"]
    return {**metric_summary, "counts": dict(counts),
            "joint_augmentation": augmentation,
            "root_covariance_condition": float(np.linalg.cond(covariance)),
            "physical_true_bias_root_fixed_not_independent": True,
            "one_correlated_initial_covariance_root": True,
            "one_sigma_motion_coordinate_retained": one_sigma_retained,
            "critical_initial_sigma_level": critical,
            "limiting_group": limiting_name,
            "limiting_prefix_ordinal": ordinal,
            "limiting_index": limiting["index"],
            "limiting_stage": limiting["stage"],
            "limiting_kind": limiting["kind"],
            "groups": summary,
            "point_centered_S_candidate": {
                "minimum_one_sigma_additive_headroom_m_s": sdiag["minimum_one_sigma_extra_euclidean_headroom"],
                "limiting_prefix_ordinal": sdiag["headroom_prefix_ordinal"],
                "limiting_index": sdiag["headroom_index"],
                "limiting_stage": sdiag["headroom_stage"],
                "limiting_kind": sdiag["headroom_kind"],
                "independent_S_port_used_for_promotion": False,
                "candidate_only_not_D_S_certificate": True,
            },
            "boundary_storage": boundary_storage,
            "maximum_one_sigma_compatible_storage_level": max_storage,
            "endpoint_one_sigma_compatible_storage_level": endpoint_storage,
            "covariance_consistency_is_unproved": True,
            "nonzero_BIAS1_driver_uniformly_enclosed": False,
            "active_projection_uniformly_enclosed": False,
            "source_dependent_coefficients_uniformly_enclosed": False,
            "centered_S_same_history_recurrence_uniformly_enclosed": False,
            "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefix", required=True, type=Path)
    ap.add_argument("--attachment", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
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
            "latent_acceleration_error_norm_upper_mps2")
    radii = np.repeat([2*np.tan(np.pi/12)] + [float(domain[k]) for k in keys], 3)
    report = {"experiment": "JOINT_SOURCE_CENTERED_CORRELATED_ROOT_RETENTION",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "proof_state_dimension": 24,
              "shipping_filter_changed": False,
              "independent_bias_slots_used": False,
              "independent_source_ports_used": False,
              "fresh_metric_per_sample_used": False,
              "covariance_ellipsoid_promoted_to_true_error_set": False,
              "point_centered_S_headroom_can_promote_D_S": False,
              "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False,
              "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached mode: "+mode)
        report["modes"][mode] = audit_mode(root, rows, points, mode, radii)
        m = report["modes"][mode]
        print("JOINT_CORRELATED_ROOT", mode,
              "retained", m["one_sigma_motion_coordinate_retained"],
              "critical_sigma", m["critical_initial_sigma_level"],
              "limiting", m["limiting_group"],
              "S_candidate_headroom", m["point_centered_S_candidate"]["minimum_one_sigma_additive_headroom_m_s"],
              flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
