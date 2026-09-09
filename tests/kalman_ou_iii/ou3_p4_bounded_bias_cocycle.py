"""Bounded internal bias reduction on an attached finite-coefficient word.

All full 21-state subevents are composed BEFORE selecting the motion block
at each sample boundary. This preserves the sample-entrance D_b convention
and within-sample corrected-bias feedback. No shipping state is reset.
The selected products are variation-of-constants operators, not estimators.

A Schur matrix for one frozen word does not establish uniform P4. The report
separates that point result, quantitative supplies, and unresolved nonlinear
coverage/retention. Decimal checks concern the binary64 coefficients only.
"""
from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

import numpy as np


def sym(a):
    return (a + a.T) / 2


def project_ball(value, radius):
    """Reference exact-real projection formula; not shipping float validation."""
    value = np.asarray(value, dtype=float)
    if value.ndim != 1 or not np.isfinite(value).all():
        raise ValueError("projection needs a finite vector")
    if not np.isfinite(radius) or radius < 0:
        raise ValueError("projection radius must be finite and nonnegative")
    scale = float(np.max(np.abs(value), initial=0.))
    if scale == 0 or radius == 0:
        return np.zeros_like(value)
    direction = value / scale
    length = float(np.linalg.norm(direction))
    if scale <= radius / length:
        return value.copy()
    return direction * (radius / length)


def sample_lifts(steps):
    """Compose each full sample, retaining every completed-event prefix.

    The returned blocks are products of full 21-state factors, not products
    of their 18-state principal blocks. A sample ends just before the next
    prediction_enter, after its last asynchronous correction if present.
    """
    samples = []
    current = None
    for step in steps:
        a = np.asarray(step["A"], dtype=float)
        forcing = np.asarray(step["B"]) @ np.asarray(step["u"])
        if a.shape != (21, 21) or forcing.shape != (21,):
            raise ValueError("full 21-state factors required")
        if not np.isfinite(a).all() or not np.isfinite(forcing).all():
            raise ValueError("nonfinite coefficient")
        if step["bias_cost"]:
            if not np.array_equal(a, np.eye(21)) or np.any(forcing):
                raise ValueError("sample entrance must preserve the mean")
            if current is not None:
                if step["index"] != current["index"] + 1:
                    raise ValueError("nonconsecutive sample entrance")
                samples.append(current)
            current = {"index": step["index"], "prefixes": [], "steps": [],
                       "full": np.eye(21), "forcing": np.zeros(21)}
        if current is None or step["index"] != current["index"]:
            raise ValueError("factor detached from its sample entrance")
        current["full"] = a @ current["full"]
        current["forcing"] = a @ current["forcing"] + forcing
        current["steps"].append(step)
        current["prefixes"].append((current["full"].copy(),
                                     current["forcing"].copy(), step))
    if current is None:
        raise ValueError("empty word")
    samples.append(current)
    return samples


def motion_scan(samples, dt):
    """Motion propagator, bias-energy Gramian and ONE source-template column.

    Bias is charged once with weight dt at the sample entrance. The Gramian
    relaxes the internal bias sequence only for a sufficient norm bound; it
    never asserts independence from the same-history motion/source graph.
    """
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("positive finite sample duration required")
    transition, gramian, response = np.eye(18), np.zeros((18, 18)), np.zeros(18)
    result = []
    for number, sample in enumerate(samples):
        for full, forcing, step in sample["prefixes"]:
            a, b = full[:18, :18], full[:18, 18:]
            item = {"transition": a @ transition,
                    "gramian": sym(a @ gramian @ a.T + b @ b.T / dt),
                    "response": a @ response + forcing[:18],
                    "sample": number, "index": sample["index"],
                    "stage": step["stage"], "kind": step["kind"]}
            result.append(item)
        transition, gramian, response = (item[k] for k in
                                         ("transition", "gramian", "response"))
    return result


def reconstruction_defect(samples, initial):
    """Check the lift with the ACTUAL full-factor corrected-bias sequence."""
    state = np.asarray(initial, dtype=float).copy()
    motion = state[:18].copy()
    worst = 0.
    for sample in samples:
        bias = state[18:].copy()
        expected = sample["full"] @ state + sample["forcing"]
        motion = (sample["full"][:18, :18] @ motion
                  + sample["full"][:18, 18:] @ bias + sample["forcing"][:18])
        worst = max(worst, float(np.max(np.abs(motion - expected[:18]))))
        state = expected
    return worst


def stein_metric(transition, radii):
    """One constructive point metric, not a grid or a source-uniform fit.

    For q=(1+spectral_radius(T))/2, solve M-(T/q)^T M(T/q)=D^-2.
    Exact-real Schur stability proves existence by the convergent series.
    This binary64 candidate and its residual are not an outward certificate.
    """
    transition = np.asarray(transition, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(transition)
    if transition.shape != (n, n) or radii.shape != (n,):
        raise ValueError("incompatible metric dimensions")
    if not np.isfinite(transition).all() or not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError("finite matrix and positive physical radii required")
    spectral = float(np.max(np.abs(np.linalg.eigvals(transition))))
    if spectral >= 1:
        raise ValueError("point motion cocycle is not Schur")
    q = (1 + spectral) / 2
    scaled = transition / q
    metric = sym(np.linalg.solve(np.eye(n*n) - np.kron(scaled.T, scaled.T),
                                 np.diag(1/radii**2).reshape(-1)).reshape(n, n))
    root = np.linalg.cholesky(metric).T
    weighted = root @ transition @ np.linalg.solve(root, np.eye(n))
    ratio = float(np.linalg.norm(weighted, 2)**2)
    residual = metric - scaled.T @ metric @ scaled - np.diag(1/radii**2)
    relative = float(np.linalg.norm(residual, 2) / np.linalg.norm(metric, 2))
    if not np.isfinite(ratio) or ratio >= 1 or relative > 1e-9:
        raise ValueError("point Stein solve failed its numerical checks")
    return metric, {"spectral_radius": spectral, "stein_q": q,
                    "complete_word_ratio": ratio, "distance_to_one": 1-ratio,
                    "relative_stein_residual": relative,
                    "metric_condition": float(np.linalg.cond(metric))}


def weighted_direction(transition, metric):
    root = np.linalg.cholesky(metric).T
    _, singular, right = np.linalg.svd(root @ transition @ np.linalg.solve(root, np.eye(len(root))))
    return float(singular[0]**2), np.linalg.solve(root, right[0])


def decimal_direction(samples, initial, metric, stop=None):
    """80-digit coefficient-direction check with signed operation costs.

    At sample boundaries E_H^T injects the performance coordinate with zero
    bias into the homogeneous operator. This is the analytical cocycle, NOT
    a modification, replay, or reset of the shipping filter.
    """
    with localcontext() as ctx:
        ctx.prec = 80
        def dec(value):
            array = np.asarray(value)
            return np.array([Decimal.from_float(float(x)) for x in array.flat],
                            dtype=object).reshape(array.shape)
        m = dec(metric)
        motion = dec(initial)
        start = motion @ m @ motion
        before = start
        costs = Counter()
        prefix, count = Decimal(1), 0
        for sample in samples:
            state = np.r_[motion, [Decimal(0)]*3]
            for step in sample["steps"]:
                state = dec(step["A"]) @ state
                energy = state[:18] @ m @ state[:18]
                costs[step["stage"] + ":" + step["kind"]] += (energy-before)/start
                before = energy
                prefix = max(prefix, energy/start)
                if stop is not None and count == stop:
                    return {"ratio": str(energy/start), "direction_prefix_max": str(prefix),
                            "signed_operation_costs": {k: str(v) for k, v in costs.items()}}
                count += 1
            motion = state[:18]
        return {"ratio": str(before/start), "direction_prefix_max": str(prefix),
                "signed_operation_costs": {k: str(v) for k, v in costs.items()}}


def supply_constants(endpoint, metric, ratio):
    """Construct separate bias/template gains by an endpoint Schur complement.

    X = M^-1 - T (rho M)^-1 T^T > 0, rho=(1+ratio)/2.
    G/gamma_b + ff^T/gamma_s < X proves the dense quadratic supply bound.
    The factor 2.01 gives strict room without any gain grid or replay fitting.
    """
    rho = (1+ratio)/2
    inverse = np.linalg.solve(metric, np.eye(18))
    t = endpoint["transition"]
    x = sym(inverse - t @ inverse @ t.T/rho)
    root = np.linalg.cholesky(x)
    whiten = np.linalg.solve(root, np.eye(18))
    gram = sym(whiten @ endpoint["gramian"] @ whiten.T)
    forcing = whiten @ endpoint["response"]
    gb = 2.01*max(0., float(np.linalg.eigvalsh(gram)[-1]))
    gs = 2.01*float(forcing @ forcing)
    normalized = np.zeros((18, 18))
    if gb:
        normalized += gram/gb
    if gs:
        normalized += np.outer(forcing, forcing)/gs
    return {"rho": rho, "gamma_bias": gb, "gamma_template": gs,
            "normalized_supply_ratio": float(np.linalg.eigvalsh(sym(normalized))[-1]),
            "template_amplitude_is_not_a_BRMM_uniform_budget": True}


def audit_steps(steps, initial, dt, radii):
    samples = sample_lifts(steps)
    scan = motion_scan(samples, dt)
    full = np.eye(21)
    for sample in samples:
        full = sample["full"] @ full
    transition = scan[-1]["transition"]
    metric, summary = stein_metric(transition, radii)
    ratios = [weighted_direction(item["transition"], metric)[0] for item in scan]
    worst = int(np.argmax(ratios))
    endpoint_ratio, endpoint_direction = weighted_direction(transition, metric)
    _, prefix_direction = weighted_direction(scan[worst]["transition"], metric)
    summary.update({"sample_count": len(samples), "completed_prefix_count": len(scan),
                    "full_word_motion_principal_spectral_radius": float(np.max(np.abs(np.linalg.eigvals(full[:18, :18])))),
                    "full_word_principal_vs_sample_cocycle_norm": float(np.linalg.norm(full[:18, :18]-transition, 2)),
                    "actual_bias_sequence_lift_defect": reconstruction_defect(samples, initial),
                    "worst_prefix_ratio": max(ratios),
                    "worst_prefix": {k: scan[worst][k] for k in ("sample", "index", "stage", "kind")},
                    "endpoint_direction": endpoint_direction.tolist(),
                    "prefix_direction": prefix_direction.tolist(),
                    "endpoint_high_precision": decimal_direction(samples, endpoint_direction, metric),
                    "prefix_high_precision": decimal_direction(samples, prefix_direction, metric, stop=worst),
                    "endpoint_supply": supply_constants(scan[-1], metric, endpoint_ratio),
                    "point_complete_word_feasible": True,
                    "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False})
    # Direct coordinate gains avoid charging all storage to the attitude row.
    names = ("attitude", "gyro_bias", "velocity", "position", "integral_displacement", "latent_acceleration")
    coordinates = {}
    for group, name in enumerate(names):
        sl = slice(3*group, 3*group+3)
        bias_gains = [float(np.sqrt(max(0., np.linalg.eigvalsh(item["gramian"][sl, sl])[-1]))) for item in scan]
        forcing = [float(np.linalg.norm(item["response"][sl])) for item in scan]
        coordinates[name] = {"max_bias_energy_gain": max(bias_gains),
                             "bias_gain_prefix": int(np.argmax(bias_gains)),
                             "max_template_excursion": max(forcing),
                             "bounded_bias_plus_template_outer_ratio": max(
                                 .4*np.sqrt(len(samples)*dt)*b + f for b, f in zip(bias_gains, forcing))/float(radii[3*group]),
                             "zero_true_bias_point_budget_only": True}
    summary["coordinate_prefix_supplies"] = coordinates
    return summary


def main():
    import ou3_p4_motion_gain as G

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
            "latent_acceleration_error_norm_upper_mps2")
    radii = np.repeat([2*np.tan(np.pi/12)] + [float(domain[k]) for k in keys], 3)
    report = {"experiment": "BOUNDED_INTERNAL_BIAS_SAMPLE_COCYCLE",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "domain_sha256": hashlib.sha256(domain_path.read_bytes()).hexdigest(),
              "frozen_coefficients_only": True, "same_history_factorization_retained": True,
              "bias_cost_convention": "dt times squared bias error at sample entrances only",
              "no_shipping_bias_reset_or_filter_change": True,
              "BIAS2_uniform_sector_used": False,
              "nonlinear_projection_invariance_falsified": False,
              "projection_ball_statement": "exact-real projection; float rounding is a separate enclosure",
              "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False, "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached mode: "+mode)
        steps, _, initial, defects, counts = G.build_word(root, rows, points, mode)
        if defects.failures:
            raise ValueError("finite factorization failed: "+repr(defects.failures[:1]))
        result = audit_steps(steps, initial, root["dt"], radii)
        result["counts"] = dict(counts)
        result["factorization_normalized_defects"] = dict(defects.defects)
        report["modes"][mode] = result
        print("BOUNDED_BIAS_COCYCLE", mode, result["spectral_radius"],
              result["complete_word_ratio"], result["worst_prefix_ratio"], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
