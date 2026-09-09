"""Joint same-history physical-bias/source recurrence with compatible storage.

This is the next P4-motion feasibility screen after the bounded-bias cocycle.
It does NOT select a fresh Lyapunov matrix for each isolated word and it does
NOT relax the accelerometer-bias error to an independent vector every sample.

The proof graph is augmented from the shipping 21-error coordinates e to

    z = (e_H[18], e_b[3], beta_true[3]) in R^24.

The shipping estimator is unchanged.  The extra beta_true coordinate is a
proof/source coordinate only.  It makes the BIAS1 recurrence and the active
radial estimate projection live in the SAME graph:

 prediction:
   beta+ = phi_true beta + w,
   e_b+  = phi_hat e_b + (phi_true-phi_hat) beta + w,

 correction/projection with scale s:
   e_b+  = s e_b_pre + (1-s) beta.

Thus a nonzero true-bias root is not silently treated as a free bias-error
history.  The current attached source declares beta_root=0 and driver=ZERO;
that literal history is verified exactly, but it cannot establish nonzero
BIAS1 driver coverage.  Future nonzero captures must declare a true-bias time
constant/driver model; this producer refuses to infer those parameters from a
trace.

For the attached coefficient/source history, full 24-state subevents are
composed within each sample before selecting motion rows.  A cyclic metric
sequence M_i is then constructed with

    M_i = A_i^T M_{i+1} A_i + Q,   M_N=M_0,

so consecutive samples use compatible storage by construction.  The actual
joint internal state (e_b,beta_true) and physical source forcing stay together
in every boundary and completed-event prefix budget.  No independent source
ports, no independent bias slots, and no source-fitted replacement observer
are introduced.

This remains a frozen-coefficient/source point test.  Uniform enclosure of the
source-dependent metric sequence, nonzero BIAS1 driver, active projection
sectors and first-exit retention is still required before any P4/P5 promotion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import ou3_p4_motion_gain as G
import ou3_p4_source_endpoint as SOURCE

I3 = np.eye(3)


def sym(a):
    return (a + a.T) / 2


def solve_stein(transition, rhs):
    """Solve M-T^T M T=rhs for one Schur transition, point arithmetic only."""
    n = transition.shape[0]
    system = np.eye(n*n) - np.kron(transition.T, transition.T)
    result = np.linalg.solve(system, np.asarray(rhs).reshape(-1)).reshape(n, n)
    return sym(result)


def projection_scales(rows, mode):
    """Recover the actual radial-projection multiplier at each correction.

    The estimate before projection is beta_true-e_b+d_b, because
    e_b=beta_true-b_hat and the state injection adds d_b to b_hat.
    """
    pending = None
    scales = {}
    for row in (r for r in rows if r["word"] == mode):
        if row["stage"] == "measurement":
            if pending is not None:
                raise ValueError("overlapping measurement before projection")
            pending = row
        elif row["stage"] == "projection":
            if pending is None or pending["index"] != row["index"] or pending["kind"] != row["kind"]:
                raise ValueError("projection detached from measurement")
            error, _ = SOURCE.error_and_covariance(pending, 21)
            correction = G.C.mat(pending, "K", 21, 3) @ np.asarray(pending["r"], dtype=float)
            beta = np.asarray(pending["true_bias"], dtype=float)
            estimate = beta - error[18:21] + correction[18:21]
            length = float(np.linalg.norm(estimate))
            radius = float(pending["projection_radius"])
            scale = min(1., radius/length) if length else 1.
            scales[(row["index"], row["kind"])] = scale
            pending = None
    if pending is not None:
        raise ValueError("unterminated measurement")
    return scales


def declared_true_phi(root, row, dt):
    """Return the declared BIAS1 physical decay for this source point.

    A zero root with ZERO driver is identically zero for any positive tau.  In
    that one degenerate case we may use the filter phi as a harmless algebraic
    representative while recording that true tau is unidentified.  A nonzero
    root is never allowed without an explicit physical tau declaration.
    """
    root_bias = np.asarray(root.get("bias_root", []), dtype=float)
    driver = root.get("bias_driver")
    tau = root.get("bias_tau_true_s")
    if tau is not None:
        tau = float(tau)
        if not math.isfinite(tau) or tau <= 0:
            raise ValueError("positive finite bias_tau_true_s required")
        return math.exp(-dt/tau), False
    if driver == "ZERO" and root_bias.shape == (3,) and np.all(root_bias == 0):
        return math.exp(-dt/float(row["tau_b"])), True
    raise ValueError("nonzero/undetermined physical bias needs declared bias_tau_true_s")


def augmented_steps(root, rows, steps, mode, dt):
    """Lift attached 21-state factors to the joint (error,true-bias) graph."""
    scales = projection_scales(rows, mode)
    result = []
    unidentified_tau = False
    for step in steps:
        a21 = np.asarray(step["A"], dtype=float)
        forcing21 = np.asarray(step["B"]) @ np.asarray(step["u"])
        a = np.eye(24)
        a[:21, :21] = a21
        forcing = np.zeros(24)
        forcing[:21] = forcing21
        if step["stage"] == "prediction":
            phi_hat = float(a21[18, 18])
            if not np.allclose(a21[18:21, 18:21], phi_hat*I3, rtol=0, atol=2e-15):
                raise ValueError("expected isotropic shipping bias prediction")
            row = next(r for r in rows if r["word"] == mode and r["index"] == step["index"]
                       and r["stage"] == "prediction")
            phi_true, unknown = declared_true_phi(root, row, dt)
            unidentified_tau = unidentified_tau or unknown
            a[18:21, 21:24] = (phi_true-phi_hat)*I3
            a[21:24, 21:24] = phi_true*I3
            if root.get("bias_driver") != "ZERO":
                raise ValueError("nonzero BIAS1 driver requires an explicit driver-input lift")
        elif step["stage"] == "projection":
            scale = scales[(step["index"], step["kind"])]
            a[18:21, 21:24] += (1-scale)*I3
            # Physical true bias is unchanged by estimator correction.
            a[21:24, 21:24] = I3
        else:
            a[21:24, 21:24] = I3
        result.append({**step, "A24": a, "forcing24": forcing})
    return result, {"true_tau_unidentified_only_because_zero_root": unidentified_tau,
                    "nonzero_driver_lift_present": False}


def sample_lifts(steps):
    """Compose full 24-state subevents before selecting performance rows."""
    samples = []
    current = None
    for step in steps:
        a = np.asarray(step["A24"], dtype=float)
        forcing = np.asarray(step["forcing24"], dtype=float)
        if a.shape != (24, 24) or forcing.shape != (24,):
            raise ValueError("joint 24-state factor required")
        if step["bias_cost"]:
            if current is not None:
                if step["index"] != current["index"]+1:
                    raise ValueError("nonconsecutive sample entrance")
                samples.append(current)
            current = {"index": step["index"], "full": np.eye(24),
                       "forcing": np.zeros(24), "prefixes": []}
        if current is None or step["index"] != current["index"]:
            raise ValueError("factor detached from sample")
        current["full"] = a @ current["full"]
        current["forcing"] = a @ current["forcing"] + forcing
        current["prefixes"].append((current["full"].copy(), current["forcing"].copy(), step))
    if current is None:
        raise ValueError("empty joint word")
    samples.append(current)
    return samples


def true_bias_recurrence(root, rows, mode):
    """Verify the literal physical root/driver continuation serialized by trace."""
    events = [r for r in rows if r["word"] == mode and r["stage"] == "prediction_enter"]
    if len(events) != 600:
        raise ValueError("complete sample-entrance physical-bias history required")
    root_bias = np.asarray(root.get("bias_root", []), dtype=float)
    driver = root.get("bias_driver")
    actual = np.asarray([r["true_bias"] for r in events], dtype=float)
    if root_bias.shape != (3,) or actual.shape != (600, 3):
        raise ValueError("invalid physical-bias coordinates")
    if driver == "ZERO" and np.all(root_bias == 0):
        defect = float(np.max(np.abs(actual)))
        passed = defect == 0.
    else:
        # Nonzero roots need a time origin/tau; nonzero drivers need their
        # serialized driver law.  Never infer either from observed extrema.
        defect = None
        passed = False
    return {"root": root_bias.tolist(), "driver": driver,
            "sample_entrance_count": len(actual),
            "recurrence_reconstruction_max_abs": defect,
            "recurrence_reconstruction_pass": passed,
            "nonzero_driver_coverage": False,
            "zero_driver_is_not_uniform_BIAS1_coverage": driver == "ZERO"}


def initial_joint_state(root, initial21):
    beta = np.asarray(root.get("bias_root", []), dtype=float)
    if beta.shape != (3,):
        raise ValueError("missing physical-bias root")
    # Current attached root is at zero and remains zero.  A future nonzero root
    # must carry a source time origin before this point helper can be promoted.
    if np.any(beta != 0):
        raise ValueError("nonzero root needs explicit root-to-word transport")
    return np.r_[np.asarray(initial21, dtype=float), beta]


def entrance_states(samples, initial):
    states = [np.asarray(initial, dtype=float).copy()]
    for sample in samples:
        states.append(sample["full"] @ states[-1] + sample["forcing"])
    return states


def periodic_metrics(samples, radii):
    """Construct one cyclic source-indexed metric sequence over all samples."""
    q = np.diag(1/np.asarray(radii, dtype=float)**2)
    prefix = np.eye(18)
    rhs = np.zeros((18, 18))
    transitions = []
    for sample in samples:
        a = sample["full"][:18, :18]
        transitions.append(a)
        rhs = sym(rhs + prefix.T @ q @ prefix)
        prefix = a @ prefix
    word = prefix
    spectral = float(np.max(np.abs(np.linalg.eigvals(word))))
    if not np.isfinite(spectral) or spectral >= 1:
        raise ValueError("joint sample motion word is not Schur")
    m0 = solve_stein(word, rhs)
    np.linalg.cholesky(m0)
    metrics = [None]*(len(samples)+1)
    metrics[-1] = m0
    for i in reversed(range(len(samples))):
        metrics[i] = sym(transitions[i].T @ metrics[i+1] @ transitions[i] + q)
        np.linalg.cholesky(metrics[i])
    closure = float(np.linalg.norm(metrics[0]-m0, 2)/np.linalg.norm(m0, 2))
    if closure > 2e-8:
        raise ValueError("cyclic compatible metric failed numerical closure")
    return metrics, q, {"word_spectral_radius": spectral,
                        "cyclic_metric_closure_relative": closure,
                        "metric_condition_max": max(float(np.linalg.cond(m)) for m in metrics)}


def generalized_floor(q, metric):
    root = np.linalg.cholesky(metric)
    congruence = np.linalg.solve(root, q) @ np.linalg.solve(root.T, np.eye(len(metric)))
    return max(0., float(np.linalg.eigvalsh(sym(congruence))[0]))


def one_step_budget(a, g, d, internal, mnext, q, theta=.5):
    """Dissipative compatible-storage bound with one correlated joint input."""
    if not 0 < theta < 1:
        raise ValueError("theta must be strictly between zero and one")
    mi = sym(a.T @ mnext @ a + q)
    mu = generalized_floor(q, mi)
    rate = 1-(1-theta)*mu
    qinv = np.linalg.solve(q, np.eye(len(q)))
    h = sym(mnext + (mnext @ a @ qinv @ a.T @ mnext)/theta)
    v = g @ internal + d
    supply = float(v @ h @ v)
    return rate, supply, v, mu


def prefix_bound(sample, metric, level, internal, radii):
    """Explicit every-completed-event physical-coordinate bound."""
    root = np.linalg.cholesky(metric)
    invroot_t = np.linalg.solve(root.T, np.eye(18))
    names = ("attitude", "gyro_bias", "velocity", "position",
             "integral_displacement", "latent_acceleration")
    result = []
    for full, forcing, step in sample["prefixes"]:
        a, g = full[:18, :18], full[:18, 18:24]
        v = g @ internal + forcing[:18]
        groups = {}
        for group, name in enumerate(names):
            sl = slice(3*group, 3*group+3)
            homogeneous = np.sqrt(max(level, 0.))*float(np.linalg.norm(a[sl] @ invroot_t, 2))
            correlated = float(np.linalg.norm(v[sl]))
            bound = homogeneous + correlated
            radius = float(radii[3*group])
            groups[name] = {"homogeneous_from_boundary_level": homogeneous,
                            "joint_bias_source_excursion": correlated,
                            "bound": bound, "declared_radius": radius,
                            "retention_ratio": bound/radius}
        result.append({"stage": step["stage"], "kind": step["kind"],
                       "index": step["index"], "groups": groups})
    return result


def joint_lift_defect(samples, states21, rows, mode):
    """Compare augmented error/beta propagation with attached 21-error history."""
    worst_error = 0.
    worst_beta = 0.
    entrances = [r for r in rows if r["word"] == mode and r["stage"] == "prediction_enter"]
    if len(entrances) != len(samples):
        raise ValueError("joint lift/trace sample mismatch")
    for i, row in enumerate(entrances):
        beta = np.asarray(row["true_bias"], dtype=float)
        worst_beta = max(worst_beta, float(np.max(np.abs(states21[i][21:24]-beta))))
    for i, sample in enumerate(samples):
        expected = sample["full"] @ states21[i] + sample["forcing"]
        worst_error = max(worst_error, float(np.max(np.abs(expected[:21]-states21[i+1][:21]))))
    return worst_error, worst_beta


def audit_mode(root, rows, points, mode, radii):
    steps21, _, initial21, defects, counts = G.build_word(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite coefficient attachment failed: "+repr(defects.failures[:1]))
    steps24, augmentation = augmented_steps(root, rows, steps21, mode, float(root["dt"]))
    samples = sample_lifts(steps24)
    initial = initial_joint_state(root, initial21)
    states = entrance_states(samples, initial)
    metrics, q, metric_summary = periodic_metrics(samples, radii)
    recurrence = true_bias_recurrence(root, rows, mode)
    error_lift, beta_lift = joint_lift_defect(samples, states, rows, mode)
    if max(error_lift, beta_lift) > 1e-10:
        raise ValueError("augmented joint graph lost the attached history")

    levels = [float(states[0][:18] @ metrics[0] @ states[0][:18])]
    exact = [levels[0]]
    boundary = []
    prefixes = []
    cumulative_supply = 0.
    product_rate = 1.
    for i, sample in enumerate(samples):
        a = sample["full"][:18, :18]
        g = sample["full"][:18, 18:24]
        d = sample["forcing"][:18]
        internal = states[i][18:24]
        rate, supply, v, mu = one_step_budget(a, g, d, internal, metrics[i+1], q)
        cumulative_supply += supply
        product_rate *= rate
        next_level = rate*levels[-1] + supply
        levels.append(next_level)
        actual_next = states[i+1][:18]
        exact.append(float(actual_next @ metrics[i+1] @ actual_next))
        boundary.append({"sample": i, "index": sample["index"],
                         "rate": rate, "dissipation_floor_mu": mu,
                         "joint_bias_source_supply": supply,
                         "joint_input_norm": float(np.linalg.norm(v)),
                         "bias_error_norm": float(np.linalg.norm(internal[:3])),
                         "physical_bias_norm": float(np.linalg.norm(internal[3:])),
                         "level_bound": next_level, "actual_storage": exact[-1]})
        prefixes.extend(prefix_bound(sample, metrics[i], levels[-2], internal, radii))

    max_prefix = {}
    for name in ("attitude", "gyro_bias", "velocity", "position",
                 "integral_displacement", "latent_acceleration"):
        candidates = [(p["groups"][name]["retention_ratio"], j, p)
                      for j, p in enumerate(prefixes)]
        ratio, j, p = max(candidates, key=lambda x: x[0])
        max_prefix[name] = {"retention_ratio": ratio, "prefix_ordinal": j,
                            "index": p["index"], "stage": p["stage"], "kind": p["kind"],
                            **p["groups"][name]}
    worst_storage_slack = max(e-b for e, b in zip(exact, levels))
    if worst_storage_slack > 5e-8*max(1., max(exact)):
        raise ValueError("compatible-storage supply failed to bound attached history")
    retained = all(item["retention_ratio"] <= 1. for item in max_prefix.values())
    return {**metric_summary, "counts": dict(counts),
            "physical_bias_recurrence": recurrence,
            "joint_augmentation": augmentation,
            "joint_error_lift_defect": error_lift,
            "joint_physical_bias_lift_defect": beta_lift,
            "boundary_count": len(boundary), "completed_prefix_count": len(prefixes),
            "compatible_storage_product_rate": product_rate,
            "total_joint_bias_source_supply": cumulative_supply,
            "initial_storage": levels[0], "endpoint_level_bound": levels[-1],
            "endpoint_actual_storage": exact[-1],
            "max_actual_minus_bound": worst_storage_slack,
            "max_prefix_coordinate_budgets": max_prefix,
            "point_motion_coordinate_retained": retained,
            "boundary_budget": boundary,
            "frozen_coefficients_only": True,
            "source_dependent_metric_uniformly_enclosed": False,
            "nonzero_BIAS1_driver_uniformly_enclosed": False,
            "active_projection_uniformly_enclosed": False,
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
    report = {"experiment": "JOINT_BIAS_SOURCE_COMPATIBLE_STORAGE_24D",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "same_history_source_bias_projection_required": True,
              "proof_state_dimension": 24,
              "shipping_error_state_dimension": 21,
              "shipping_filter_changed": False,
              "independent_bias_energy_slots_used": False,
              "independent_source_ports_used": False,
              "fresh_metric_per_sample_used": False,
              "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False,
              "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached mode: "+mode)
        report["modes"][mode] = audit_mode(root, rows, points, mode, radii)
        m = report["modes"][mode]
        print("JOINT_BIAS_SOURCE_STORAGE", mode,
              "rate", m["compatible_storage_product_rate"],
              "retained", m["point_motion_coordinate_retained"],
              "max_prefix", max(v["retention_ratio"] for v in m["max_prefix_coordinate_budgets"].values()),
              flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
