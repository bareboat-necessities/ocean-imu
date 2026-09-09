"""Joint same-history bias/source recurrence with compatible motion storage.

This is the next P4-motion feasibility screen after the bounded-bias cocycle.
It deliberately does NOT pick a fresh Lyapunov matrix for each isolated word.
For one SHA-attached coefficient/source history it constructs a cyclic metric
sequence M_i satisfying

    M_i = A_i^T M_{i+1} A_i + Q,       M_N = M_0,

where A_i is the motion block of the FULL composed 21-state shipping sample.
Thus storage before/after consecutive samples is compatible by construction.
The actual corrected bias at each sample entrance and the exact correlated
physical source continuation are then retained together in

    x_{i+1} = A_i x_i + G_i b_i + d_i.

No independent per-sample bias ball or source port is introduced.  Every
completed-event prefix receives an explicit physical-coordinate bound from
the same boundary level and the same actual bias/source continuation.

The current source endpoint has the literal BIAS1 root/driver ZERO.  This file
checks that recurrence rather than silently replacing it with arbitrary bias
slots.  A zero-driver point does not certify nonzero BIAS1 coverage; all P4/P5
promotion flags remain false until source-uniform driver/projection enclosure
and retained levels are proved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_bounded_bias_cocycle as C
import ou3_p4_motion_gain as G
import ou3_p4_source_endpoint as SOURCE


def sym(a):
    return (a + a.T) / 2


def solve_stein(transition, rhs):
    """Solve M-T^T M T=rhs for one Schur transition, point arithmetic only."""
    n = transition.shape[0]
    system = np.eye(n*n) - np.kron(transition.T, transition.T)
    result = np.linalg.solve(system, np.asarray(rhs).reshape(-1)).reshape(n, n)
    return sym(result)


def periodic_metrics(samples, radii):
    """Construct one cyclic metric sequence over the complete sample history.

    Q is fixed in physical coordinates.  The full-period Stein equation gives
    M_0.  Backward recursion gives each M_i and returns to the SAME M_0, so
    adjacent samples never get independently selected endpoint metrics.
    """
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
        raise ValueError("sample motion word is not Schur; cyclic metric unavailable")
    m0 = solve_stein(word, rhs)
    np.linalg.cholesky(m0)
    metrics = [None] * (len(samples)+1)
    metrics[-1] = m0
    for i in reversed(range(len(samples))):
        metrics[i] = sym(transitions[i].T @ metrics[i+1] @ transitions[i] + q)
        np.linalg.cholesky(metrics[i])
    closure = float(np.linalg.norm(metrics[0]-m0, 2)/np.linalg.norm(m0, 2))
    if closure > 2e-8:
        raise ValueError("cyclic metric failed numerical closure")
    return metrics, q, {"word_spectral_radius": spectral,
                        "cyclic_metric_closure_relative": closure,
                        "metric_condition_max": max(float(np.linalg.cond(m)) for m in metrics)}


def true_bias_recurrence(root, rows, mode):
    """Verify the literal physical bias root/driver represented by this trace.

    This source capture currently declares a zero root and ZERO driver.  That
    is a valid BIAS1 history but not a substitute for nonzero driver coverage.
    The check is intentionally exact at the serialized binary64 coordinates.
    """
    events = [row for row in rows if row["word"] == mode and row["stage"] == "prediction_enter"]
    if len(events) != 600:
        raise ValueError("complete sample-entrance bias history required")
    root_bias = np.asarray(root.get("bias_root", []), dtype=float)
    driver = root.get("bias_driver")
    if root_bias.shape != (3,) or not np.isfinite(root_bias).all():
        raise ValueError("missing finite BIAS1 physical root")
    actual = np.asarray([row["true_bias"] for row in events], dtype=float)
    if actual.shape != (600, 3) or not np.isfinite(actual).all():
        raise ValueError("missing true-bias continuation")
    if driver == "ZERO":
        expected = np.repeat(root_bias[None, :], len(actual), axis=0)
        defect = float(np.max(np.abs(actual-expected)))
        passed = defect == 0.
    else:
        # A future nonzero source must declare enough driver parameters for an
        # independent recurrence reconstruction.  Do not infer them from the
        # observed trajectory or silently downgrade to a point fit.
        defect = None
        passed = False
    return {"root": root_bias.tolist(), "driver": driver,
            "sample_entrance_count": len(actual),
            "recurrence_reconstruction_max_abs": defect,
            "recurrence_reconstruction_pass": passed,
            "nonzero_driver_coverage": bool(driver != "ZERO" and passed),
            "zero_driver_is_not_uniform_BIAS1_coverage": driver == "ZERO"}


def entrance_states(samples, initial):
    """Propagate the exact full affine coefficient history at sample boundaries."""
    states = [np.asarray(initial, dtype=float).copy()]
    for sample in samples:
        states.append(sample["full"] @ states[-1] + sample["forcing"])
    return states


def generalized_floor(q, metric):
    root = np.linalg.cholesky(metric)
    scaled = np.linalg.solve(root, q) @ np.linalg.solve(root.T, np.eye(len(metric)))
    return max(0., float(np.linalg.eigvalsh(sym(scaled))[0]))


def one_step_budget(a, g, d, bias, mnext, q, theta=.5):
    """Compatible-storage dissipative bound with correlated actual input.

    From M_i=A^T M_next A+Q and Young's inequality,

      W+ <= [1-(1-theta) mu_i] W
            + v^T[M_next + theta^-1 M_next A Q^-1 A^T M_next]v,

    v=G*b+d and Q >= mu_i M_i.  Only v from the actual joint continuation is
    evaluated here; no independent bias/source energy relaxation is used.
    """
    if not 0 < theta < 1:
        raise ValueError("theta must be strictly between zero and one")
    mi = sym(a.T @ mnext @ a + q)
    mu = generalized_floor(q, mi)
    rate = 1-(1-theta)*mu
    qi = np.linalg.solve(q, np.eye(len(q)))
    h = sym(mnext + (mnext @ a @ qi @ a.T @ mnext)/theta)
    v = g @ bias + d
    supply = float(v @ h @ v)
    return rate, supply, v, mu


def prefix_bound(sample, metric, level, bias, radii):
    """Bound every completed event from one compatible sample-boundary level.

    For W_i<=L and x_l=A_l x_i+v_l,
      ||Pi_G x_l|| <= sqrt(L)||Pi_G A_l M_i^-1/2||_2 + ||Pi_G v_l||.
    This keeps the actual same-history bias/source vector together.  It is a
    point coefficient bound, not a nonlinear source-uniform first-exit proof.
    """
    root = np.linalg.cholesky(metric)
    invroot_t = np.linalg.solve(root.T, np.eye(18))
    names = ("attitude", "gyro_bias", "velocity", "position",
             "integral_displacement", "latent_acceleration")
    report = []
    for full, forcing, step in sample["prefixes"]:
        a, g = full[:18, :18], full[:18, 18:]
        v = g @ bias + forcing[:18]
        groups = {}
        for group, name in enumerate(names):
            sl = slice(3*group, 3*group+3)
            homogeneous = np.sqrt(max(level, 0.))*float(np.linalg.norm(a[sl] @ invroot_t, 2))
            correlated = float(np.linalg.norm(v[sl]))
            bound = homogeneous + correlated
            groups[name] = {"homogeneous_from_boundary_level": homogeneous,
                            "joint_bias_source_excursion": correlated,
                            "bound": bound,
                            "declared_radius": float(radii[3*group]),
                            "retention_ratio": bound/float(radii[3*group])}
        report.append({"stage": step["stage"], "kind": step["kind"],
                       "index": step["index"], "groups": groups})
    return report


def audit_mode(root, rows, points, mode, radii):
    steps, _, initial, defects, counts = G.build_word(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite coefficient attachment failed: "+repr(defects.failures[:1]))
    samples = C.sample_lifts(steps)
    states = entrance_states(samples, initial)
    metrics, q, metric_summary = periodic_metrics(samples, radii)
    recurrence = true_bias_recurrence(root, rows, mode)

    # Verify that the bias used by the reduced motion identity is the ACTUAL
    # corrected bias error from the full same-history factorization.
    lift = C.reconstruction_defect(samples, initial)
    if lift > 1e-10:
        raise ValueError("joint motion/bias lift lost the full-state history")

    levels = [float(states[0][:18] @ metrics[0] @ states[0][:18])]
    exact = [levels[0]]
    boundary = []
    prefixes = []
    cumulative_supply = 0.
    product_rate = 1.
    for i, sample in enumerate(samples):
        a = sample["full"][:18, :18]
        g = sample["full"][:18, 18:]
        d = sample["forcing"][:18]
        bias = states[i][18:]
        rate, supply, v, mu = one_step_budget(a, g, d, bias, metrics[i+1], q)
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
                         "bias_error_norm": float(np.linalg.norm(bias)),
                         "level_bound": next_level,
                         "actual_storage": exact[-1]})
        prefixes.extend(prefix_bound(sample, metrics[i], levels[-2], bias, radii))

    max_prefix = {}
    for name in ("attitude", "gyro_bias", "velocity", "position",
                 "integral_displacement", "latent_acceleration"):
        candidates = [(p["groups"][name]["retention_ratio"], j, p)
                      for j, p in enumerate(prefixes)]
        ratio, j, p = max(candidates, key=lambda x: x[0])
        max_prefix[name] = {"retention_ratio": ratio, "prefix_ordinal": j,
                            "index": p["index"], "stage": p["stage"], "kind": p["kind"],
                            **p["groups"][name]}

    # Periodic compatible storage should upper-bound the exact point history.
    worst_storage_slack = max(e-b for e, b in zip(exact, levels))
    if worst_storage_slack > 5e-8*max(1., max(exact)):
        raise ValueError("compatible-storage supply failed to bound attached history")
    motion_retained = all(item["retention_ratio"] <= 1. for item in max_prefix.values())
    return {**metric_summary, "counts": dict(counts),
            "physical_bias_recurrence": recurrence,
            "actual_bias_sequence_lift_defect": lift,
            "boundary_count": len(boundary), "completed_prefix_count": len(prefixes),
            "compatible_storage_product_rate": product_rate,
            "total_joint_bias_source_supply": cumulative_supply,
            "initial_storage": levels[0], "endpoint_level_bound": levels[-1],
            "endpoint_actual_storage": exact[-1],
            "max_actual_minus_bound": worst_storage_slack,
            "max_prefix_coordinate_budgets": max_prefix,
            "point_motion_coordinate_retained": motion_retained,
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
    report = {"experiment": "JOINT_BIAS_SOURCE_COMPATIBLE_STORAGE",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "same_history_source_bias_projection_required": True,
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
