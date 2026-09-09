"""Nonzero BIAS1 driver on the same-history 24D P4 point graph.

The physical bias is not reduced to a free per-sample error input.  The source
capture declares

    beta(t) = beta0 exp(-t/tau_true) + a sin(omega t + phase)

componentwise.  At every IMU prediction this producer uses the exact declared
affine recurrence

    beta_i = phi_true beta_{i-1} + w_i,
    w_i = beta(t_i) - phi_true beta(t_{i-1}),

and injects the SAME w_i into both the corrected bias error and the proof-only
true-bias coordinate:

    e_b,i = phi_hat e_b,i-1 + (phi_true-phi_hat) beta_i-1 + w_i.

Estimator corrections and radial projection remain in the same graph.  For a
projection multiplier s on the attached history,

    e_b+ = s e_b_pre + (1-s) beta_true.

A local finite 21-error factorization first verifies the driven shipping trace.
Its physical-bias affine terms are then lifted back into the beta_true state,
so they are not double counted as exogenous forcing.  Physical source forcing,
bias recurrence, projection, covariance feedback and all actual R_S therefore
remain one continuation.  No independent source port or fresh bias slot is
introduced.

The compatible cyclic motion metric is reused only as a source-indexed storage
coordinate.  Retention is evaluated directly at every completed event from one
correlated root covariance ellipsoid about the actual driven center.  This is
still a frozen-coefficient point test: that ellipsoid is not a qualified hard
entry set and the source/projection coefficient family is not uniformly
enclosed.  P4/P5 flags therefore stay false.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import ou3_p4_joint_bias_source_storage as J
import ou3_p4_motion_gain as G
import ou3_p4_source_endpoint as SOURCE

I3 = np.eye(3)
I21 = np.eye(21)


def declared_bias(root, t):
    beta0 = np.asarray(root.get("bias_root", []), dtype=float)
    amp = np.asarray(root.get("bias_driver_amplitude", []), dtype=float)
    phase = np.asarray(root.get("bias_driver_phase", []), dtype=float)
    tau = float(root.get("bias_tau_true_s", float("nan")))
    omega = float(root.get("bias_driver_omega_rad_s", float("nan")))
    if root.get("bias_driver") != "DETERMINISTIC_SINUSOIDAL_GM":
        raise ValueError("this diagnostic requires DETERMINISTIC_SINUSOIDAL_GM")
    if beta0.shape != (3,) or amp.shape != (3,) or phase.shape != (3,):
        raise ValueError("complete declared bias driver vectors required")
    if not math.isfinite(tau) or tau <= 0 or not math.isfinite(omega):
        raise ValueError("finite positive tau and finite driver frequency required")
    return beta0*np.exp(-float(t)/tau) + amp*np.sin(omega*float(t)+phase)


def driver_increment(root, t_after, dt):
    tau = float(root["bias_tau_true_s"])
    phi = math.exp(-float(dt)/tau)
    before = declared_bias(root, float(t_after)-float(dt))
    after = declared_bias(root, float(t_after))
    return phi, after-phi*before


def trace_recurrence(root, rows, mode, dt):
    selected = [r for r in rows if r["word"] == mode]
    predictions = [r for r in selected if r["stage"] == "prediction"]
    entrances = [r for r in selected if r["stage"] == "prediction_enter"]
    if len(predictions) != 600 or len(entrances) != 600:
        raise ValueError("600 prediction and entrance events required")
    worst_model = 0.
    worst_recurrence = 0.
    max_driver = 0.
    for enter, pred in zip(entrances, predictions):
        if enter["index"] != pred["index"]:
            raise ValueError("prediction detached from entrance")
        before = np.asarray(enter["true_bias"], dtype=float)
        after = np.asarray(pred["true_bias"], dtype=float)
        t_after = float(pred["source_time"])
        declared_before = declared_bias(root, t_after-dt)
        declared_after = declared_bias(root, t_after)
        phi, w = driver_increment(root, t_after, dt)
        worst_model = max(worst_model,
                          float(np.max(np.abs(before-declared_before))),
                          float(np.max(np.abs(after-declared_after))))
        worst_recurrence = max(worst_recurrence,
                               float(np.max(np.abs(after-(phi*before+w)))))
        max_driver = max(max_driver, float(np.linalg.norm(w)))
    return {"driver": root["bias_driver"],
            "bias_tau_true_s": float(root["bias_tau_true_s"]),
            "sample_count": len(predictions),
            "declared_model_trace_max_abs_defect": worst_model,
            "affine_recurrence_max_abs_defect": worst_recurrence,
            "maximum_driver_increment_norm": max_driver,
            "nonzero_driver_executed": max_driver > 0.}


def _append_affine_port(b, u, covariance, affine):
    """Append one deterministic scalar-one port for finite point parity only."""
    affine = np.asarray(affine, dtype=float)
    if affine.shape != (21,):
        raise ValueError("21-vector affine source required")
    b2 = np.column_stack((np.asarray(b, dtype=float), affine))
    u2 = np.r_[np.asarray(u, dtype=float), 1.]
    n = covariance.shape[0]
    c2 = np.zeros((n+1, n+1))
    c2[:n, :n] = covariance
    return b2, u2, c2


def measurement_factor_nonzero(row, error):
    """Exact finite correction/projection factor at a nonzero true-bias point.

    The returned deterministic affine column is used only to verify the 21-D
    trace.  `projection_bias_affine` is removed again when lifting to 24-D and
    replaced by the beta_true state coupling (1-s)I.
    """
    residual, source, h = G.C.residual_graph(row, error)
    if row["kind"] in ("accelerometer", "magnetometer"):
        h[:, :3] = np.linalg.solve(I3-.5*G.C.skew(error[:3]), h[:, :3])
        if row["kind"] == "accelerometer":
            h[:, 15:18] = G.C.rotation(error[:3])@G.C.mat(row, "R_hat")
    k = G.C.mat(row, "K", 21, 3)
    correction = k@residual
    w, q = G.quat_coefficients(correction[:3]@correction[:3])
    scale_att = 2*q/w
    denominator = 1+scale_att*(error[:3]@correction[:3])/4
    if denominator <= 0:
        raise ValueError("finite correction left the retained Cayley chart")
    a, b = I21-k@h, -k.copy()
    reset_gain = scale_att*(I3+.5*G.C.skew(error[:3]))@k[:3]
    a[:3] = (I21[:3]-reset_gain@h)/denominator
    b[:3] = -reset_gain/denominator

    beta = np.asarray(row["true_bias"], dtype=float)
    error_pre = error[18:21]-correction[18:21]
    estimate_pre = beta-error_pre
    norm = float(np.linalg.norm(estimate_pre))
    radius = float(row["projection_radius"])
    projection = min(1., radius/norm) if norm else 1.
    a[18:21] *= projection
    b[18:21] *= projection
    affine = np.zeros(21)
    affine[18:21] = (1-projection)*beta
    b, source, covariance = _append_affine_port(
        b, source, G.C.mat(row, "R"), affine)
    return a, b, source, covariance, float(denominator), affine, projection


def prediction_factor_nonzero(before, row, dt):
    """Shipping prediction plus the exact physical-bias affine point term."""
    a, b, u, covariance, denominator = G.prediction_factor(before, row, dt)
    phi_hat = math.exp(-dt/float(row["tau_b"])) if row["active"] else 1.
    beta_before = np.asarray(before["true_bias"], dtype=float)
    beta_after = np.asarray(row["true_bias"], dtype=float)
    affine = np.zeros(21)
    affine[18:21] = beta_after-phi_hat*beta_before
    b, u, covariance = _append_affine_port(b, u, covariance, affine)
    return a, b, u, covariance, denominator, affine


def build_word_nonzero(root, rows, points, mode):
    """Build/verify the driven 21-D finite factors without a zero-bias premise."""
    events = [r for r in rows if r["word"] == mode]
    original = [r for r in points if r["word"] == mode]
    if not original or not events or events[0]["stage"] != "prediction_enter":
        raise ValueError("missing connected word root")
    e0, _ = SOURCE.error_and_covariance(original[0], 21)
    running = e0.copy()
    steps, counts, defects = [], Counter(), G.C.ParityChecks()
    previous = original[0]
    entrance = pending = None
    for row in events:
        stage = row["stage"]
        projection_affine = None
        prediction_affine = None
        projection_scale = None
        if stage == "prediction_enter":
            entrance = row
            a, b, u = I21.copy(), np.zeros((21, 0)), np.zeros(0)
            covariance, denominator = np.zeros((0, 0)), 1.
        elif stage == "prediction":
            a, b, u, covariance, denominator, prediction_affine = \
                prediction_factor_nonzero(entrance, row, float(root["dt"]))
        elif stage == "aw_floor":
            a, b, u = I21.copy(), np.zeros((21, 0)), np.zeros(0)
            covariance, denominator = np.zeros((0, 0)), 1.
        elif stage == "measurement":
            pending = row
            continue
        elif stage == "projection":
            if pending is None:
                raise ValueError("projection without pending measurement")
            error, _ = SOURCE.error_and_covariance(pending, 21)
            (a, b, u, covariance, denominator,
             projection_affine, projection_scale) = measurement_factor_nonzero(pending, error)
            if row["kind"] == "S_zero" and pending["R"] != pending["R_S"]:
                raise ValueError("derived gain lost actual R_S")
            counts[row["kind"]] += 1
        else:
            continue

        actual, _ = SOURCE.error_and_covariance(row, 21)
        before, _ = SOURCE.error_and_covariance(previous, 21)
        defects.row = row
        defects("finite_factorization", a@before+b@u, actual,
                1+np.linalg.norm(before, np.inf))
        running = a@running+b@u
        defects("composed_word_factorization", running, actual,
                1+np.linalg.norm(actual, np.inf))
        steps.append({"A": a, "B": b, "u": u, "Uinv": covariance,
                      "bias2": None, "metric": G.metric(row),
                      "bias_cost": stage == "prediction_enter",
                      "index": row["index"], "stage": stage, "kind": row["kind"],
                      "denominator": float(denominator),
                      "prediction_bias_affine": prediction_affine,
                      "projection_bias_affine": projection_affine,
                      "projection_scale": projection_scale})
        counts[stage] += 1
        previous = row
        if stage == "projection":
            pending = None
    if counts["prediction"] != 600 or counts["accelerometer"] != 600:
        raise ValueError("incomplete finite-factor word")
    return steps, G.metric(original[0]), e0, defects, counts


def augmented_steps(root, rows, steps, mode, dt):
    """Lift driven shipping factors into the joint (error,true-bias) graph."""
    # Synthetic unit tests may supply steps without point-factor metadata; in
    # that case recover the attached scale as before.
    fallback_scales = None
    row_by_prediction = {(r["index"], r["stage"]): r for r in rows
                         if r["word"] == mode and r["stage"] == "prediction"}
    lifted = []
    driver_norms = []
    projection_scales = []
    for step in steps:
        a21 = np.asarray(step["A"], dtype=float)
        forcing21 = np.asarray(step["B"]) @ np.asarray(step["u"])
        a = np.eye(24)
        a[:21, :21] = a21
        forcing = np.zeros(24)
        forcing[:21] = forcing21
        if step["stage"] == "prediction":
            point_affine = step.get("prediction_bias_affine")
            if point_affine is not None:
                forcing[:21] -= np.asarray(point_affine, dtype=float)
            phi_hat = float(a21[18, 18])
            if not np.allclose(a21[18:21, 18:21], phi_hat*I3, rtol=0, atol=2e-15):
                raise ValueError("expected isotropic shipping bias prediction")
            row = row_by_prediction[(step["index"], "prediction")]
            phi_true, w = driver_increment(root, float(row["source_time"]), dt)
            a[18:21, 21:24] = (phi_true-phi_hat)*I3
            a[21:24, 21:24] = phi_true*I3
            forcing[18:21] += w
            forcing[21:24] += w
            driver_norms.append(float(np.linalg.norm(w)))
        elif step["stage"] == "projection":
            point_affine = step.get("projection_bias_affine")
            if point_affine is not None:
                forcing[:21] -= np.asarray(point_affine, dtype=float)
            scale = step.get("projection_scale")
            if scale is None:
                if fallback_scales is None:
                    fallback_scales = J.projection_scales(rows, mode)
                scale = fallback_scales[(step["index"], step["kind"])]
            scale = float(scale)
            if not 0 < scale <= 1:
                raise ValueError("invalid radial projection scale")
            a[18:21, 21:24] += (1-scale)*I3
            a[21:24, 21:24] = I3
            projection_scales.append(scale)
        else:
            a[21:24, 21:24] = I3
        lifted.append({**step, "A24": a, "forcing24": forcing})
    return lifted, {"driver_prediction_count": len(driver_norms),
                    "maximum_driver_increment_norm": max(driver_norms, default=0.),
                    "projection_event_count": len(projection_scales),
                    "minimum_projection_scale": min(projection_scales, default=1.),
                    "saturated_projection_count": sum(s < 1. for s in projection_scales),
                    "driver_affine_input_lift_present": True,
                    "projection_beta_coupling_present": bool(projection_scales)}


def first_entrance_beta(rows, mode):
    row = next((r for r in rows if r["word"] == mode and r["stage"] == "prediction_enter"), None)
    if row is None:
        raise ValueError("missing first prediction entrance")
    return np.asarray(row["true_bias"], dtype=float)


def root_covariance(points, mode):
    point = next((p for p in points if p["word"] == mode and p.get("event") == "root"), None)
    if point is None:
        raise ValueError("missing attached word root")
    _, covariance = SOURCE.error_and_covariance(point, 21)
    chol = np.linalg.cholesky(covariance)
    lift = np.zeros((24, 21))
    lift[:21] = chol
    return covariance, lift


def audit_mode(root, rows, points, mode, radii):
    dt = float(root["dt"])
    recurrence = trace_recurrence(root, rows, mode, dt)
    if recurrence["declared_model_trace_max_abs_defect"] > 2e-12:
        raise ValueError("trace is detached from declared physical-bias driver")
    if recurrence["affine_recurrence_max_abs_defect"] > 2e-12:
        raise ValueError("declared physical-bias recurrence does not reconstruct trace")
    if not recurrence["nonzero_driver_executed"]:
        raise ValueError("nonzero driver diagnostic received zero driver")

    steps21, _, initial21, defects, counts = build_word_nonzero(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite coefficient attachment failed: "+repr(defects.failures[:1]))
    steps24, augmentation = augmented_steps(root, rows, steps21, mode, dt)
    samples = J.sample_lifts(steps24)
    initial = np.r_[np.asarray(initial21, dtype=float), first_entrance_beta(rows, mode)]
    states = J.entrance_states(samples, initial)
    metrics, _, metric_summary = J.periodic_metrics(samples, radii)

    entrances = [r for r in rows if r["word"] == mode and r["stage"] == "prediction_enter"]
    if len(entrances) != len(samples):
        raise ValueError("joint lift/trace sample mismatch")
    beta_defect = max(float(np.max(np.abs(states[i][21:24]-np.asarray(row["true_bias"], dtype=float))))
                      for i, row in enumerate(entrances))
    if beta_defect > 2e-10:
        raise ValueError("24D recurrence lost true-bias continuation")

    covariance, root_lift = root_covariance(points, mode)
    names = ("attitude", "gyro_bias", "velocity", "position",
             "integral_displacement", "latent_acceleration")
    global_transition = np.eye(24)
    prefixes = []
    boundary_storage = []
    weighted_root = np.linalg.cholesky(metrics[0]).T @ root_lift[:18]
    boundary_storage.append({"boundary": 0,
                             "one_sigma_storage_level": float(np.linalg.norm(weighted_root, 2)**2)})

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
                critical = (headroom/gain if gain > 0 and headroom > 0
                            else (float("inf") if gain == 0 and headroom > 0 else 0.))
                groups[name] = {"driven_center_norm": center_norm,
                                "one_sigma_correlated_root_gain": gain,
                                "one_sigma_retention_ratio": (center_norm+gain)/radius,
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
    critical_all = []
    for name in names:
        worst = max(((p["groups"][name]["one_sigma_retention_ratio"], j, p)
                     for j, p in enumerate(prefixes)), key=lambda x: x[0])
        limiting = min(((p["groups"][name]["critical_initial_sigma_level"], j, p)
                        for j, p in enumerate(prefixes)), key=lambda x: x[0])
        ratio, ordinal, wp = worst
        critical, cordinal, cp = limiting
        critical_all.append((critical, name, cordinal, cp))
        summary[name] = {"maximum_one_sigma_retention_ratio": ratio,
                         "worst_prefix_ordinal": ordinal,
                         "worst_index": wp["index"], "worst_stage": wp["stage"],
                         "worst_kind": wp["kind"],
                         "critical_initial_sigma_level": critical,
                         "critical_prefix_ordinal": cordinal,
                         "critical_index": cp["index"], "critical_stage": cp["stage"],
                         "critical_kind": cp["kind"]}
    critical, limiting_name, limiting_ordinal, limiting_prefix = min(critical_all, key=lambda x: x[0])
    retained = all(v["maximum_one_sigma_retention_ratio"] <= 1 for v in summary.values())
    return {**metric_summary, "counts": dict(counts),
            "physical_bias_recurrence": recurrence,
            "joint_augmentation": augmentation,
            "joint_physical_bias_lift_defect": beta_defect,
            "finite_factorization_normalized_defects": dict(defects.defects),
            "root_covariance_condition": float(np.linalg.cond(covariance)),
            "one_correlated_initial_covariance_root": True,
            "physical_bias_driver_is_same_history_not_independent_port": True,
            "active_projection_is_same_24D_continuation": True,
            "completed_prefix_count": len(prefixes),
            "one_sigma_motion_coordinate_retained": retained,
            "critical_initial_sigma_level": critical,
            "limiting_group": limiting_name,
            "limiting_prefix_ordinal": limiting_ordinal,
            "limiting_index": limiting_prefix["index"],
            "limiting_stage": limiting_prefix["stage"],
            "limiting_kind": limiting_prefix["kind"],
            "groups": summary,
            "maximum_one_sigma_compatible_storage_level": max(x["one_sigma_storage_level"] for x in boundary_storage),
            "endpoint_one_sigma_compatible_storage_level": boundary_storage[-1]["one_sigma_storage_level"],
            "covariance_consistency_is_unproved": True,
            "nonzero_BIAS1_driver_point_covered": True,
            "nonzero_BIAS1_driver_uniformly_enclosed": False,
            "active_projection_point_graph_retained": True,
            "active_projection_uniformly_enclosed": False,
            "source_dependent_coefficients_uniformly_enclosed": False,
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
    report = {"experiment": "JOINT_NONZERO_BIAS_DRIVER_SOURCE_CENTERED_RETENTION_24D",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "shipping_filter_changed": False,
              "proof_state_dimension": 24,
              "independent_bias_slots_used": False,
              "independent_source_ports_used": False,
              "fresh_metric_per_sample_used": False,
              "nonzero_driver_required": True,
              "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False,
              "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached mode: "+mode)
        report["modes"][mode] = audit_mode(root, rows, points, mode, radii)
        m = report["modes"][mode]
        print("JOINT_NONZERO_DRIVER", mode,
              "retained", m["one_sigma_motion_coordinate_retained"],
              "critical_sigma", m["critical_initial_sigma_level"],
              "limiting", m["limiting_group"],
              "projection_min", m["joint_augmentation"]["minimum_projection_scale"],
              flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
