"""Three complete-word storage diagnostics on one attached physical history.

All finite 21-state factors, actual covariances, corrected residuals, resets
and projections come from the shipping trace. Coefficients are frozen only
for this diagnostic. No metric/source fit, domain reduction or promotion.

Route 1 uses x_j = T_j x_0 + r_j, r_0=0, and separately budgets initial
motion, initial bias and the entire common forcing template. Route 2 tests
fixed physical storages; changing coordinates alone is an invariant control.
Route 3 transports actual accepted vector information through T and retains
the exact signed remainder in M_0 - T_N^T M_N T_N = J + remainder.
Neither a positive Gramian nor a contracting centered endpoint bounds r_j
uniformly over physical sources or the nonlinear coefficient dependencies.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_connected_motion as C
import ou3_p4_motion_gain as G


def whitener(metric):
    return np.linalg.solve(np.linalg.cholesky(G.sym(metric)).T,
                           np.eye(len(metric)))


def largest_ratio(transition, before, after):
    """Generalized eigenvector in physical root coordinates, before energy 1."""
    root = whitener(before)
    normalized = np.linalg.cholesky(G.sym(after)).T @ transition @ root
    values, vectors = np.linalg.eigh(G.sym(normalized.T @ normalized))
    return float(values[-1]), root @ vectors[:, -1]


def physical_scales(gravity, seconds):
    """Fixed units: 1 rad, 1/T rad/s, g*T, g*T^2, g*T^3, g, g."""
    if gravity <= 0 or seconds <= 0:
        raise ValueError("physical units require positive gravity and time")
    return np.repeat([1., 1/seconds, gravity*seconds, gravity*seconds**2,
                      gravity*seconds**3, gravity, gravity], 3)


def finite_h(row):
    error, _ = C.SOURCE.error_and_covariance(row, 21)
    _, _, h = C.residual_graph(row, error)
    if row["kind"] in ("accelerometer", "magnetometer"):
        h[:, :3] = np.linalg.solve(np.eye(3)-.5*C.skew(error[:3]), h[:, :3])
        if row["kind"] == "accelerometer":
            h[:, 15:18] = C.rotation(error[:3]) @ C.mat(row, "R_hat")
    return h


def compose(steps, events):
    """All source ports share one amplitude; none are dropped or decorrelated."""
    transitions, responses, metrics, measurements = [], [], [], []
    transition, response = np.eye(21), np.zeros(21)
    completed = [r for r in events if r["stage"] in
                 ("prediction_enter", "prediction", "aw_floor", "projection")]
    measured = {(r["index"], r["kind"]): r for r in events
                if r["stage"] == "measurement"}
    if len(completed) != len(steps):
        raise ValueError("incomplete event/step alignment")
    for step, event in zip(steps, completed, strict=True):
        if any(step[k] != event[k] for k in ("index", "stage", "kind")):
            raise ValueError("reordered source events")
        measurement = None
        if step["stage"] == "projection":
            row = measured[(step["index"], step["kind"])]
            measurement = {"H": finite_h(row), "R": C.mat(row, "R"),
                           "kind": row["kind"]}
        transition = step["A"] @ transition
        response = step["A"] @ response + step["B"] @ step["u"]
        _, covariance = C.SOURCE.error_and_covariance(event, 21)
        metrics.append(G.sym(np.linalg.solve(covariance, np.eye(21))))
        transitions.append(transition.copy())
        responses.append(response.copy())
        measurements.append(measurement)
    return transitions, responses, metrics, measurements


def decimal_array(value):
    value = np.asarray(value)
    return np.array([Decimal.from_float(float(x)) for x in value.flat],
                    dtype=object).reshape(value.shape)


def witness(steps, m0, metrics, direction, stop, measurements=None):
    """80-digit propagation of the supplied maximizing direction, not intervals."""
    n = len(m0)
    with localcontext() as ctx:
        ctx.prec = 80
        x = decimal_array(np.r_[direction, np.zeros(21-n)])
        initial = x[:n] @ decimal_array(m0) @ x[:n]
        previous = initial
        changes, information = defaultdict(Decimal), defaultdict(Decimal)
        for i in range(stop+1):
            step = steps[i]
            measurement = measurements[i] if measurements is not None else None
            if measurement is not None:
                y = decimal_array(measurement["H"]) @ x
                rinv = decimal_array(np.linalg.solve(measurement["R"], np.eye(3)))
                information[measurement["kind"]] += y @ rinv @ y
            x = decimal_array(step["A"]) @ x
            energy = x[:n] @ decimal_array(metrics[i]) @ x[:n]
            changes[step["stage"]+":"+step["kind"]] += (energy-previous)/initial
            previous = energy
        ratio = previous/initial
        info = {k: str(v/initial) for k, v in information.items()}
        return {"root_direction_full21": np.r_[direction, np.zeros(21-n)].tolist(),
                "terminal_direction_full21": [str(v) for v in x],
                "ratio_80_digit": str(ratio), "distance_to_one_80_digit": str(1-ratio),
                "operation_energy_changes_80_digit": {k: str(v) for k, v in changes.items()},
                "transported_information_80_digit": info,
                "signed_remainder_80_digit": str(1-ratio-sum(information.values())/initial),
                "nonlinear_admissibility_of_direction_established": False}


def evaluate(steps, transitions, responses, m0, metrics, measurements=None):
    n = len(m0)
    ratios, forced, bias_gains = [], [], []
    worst, worst_direction = -1., None
    for i, (transition, response, metric) in enumerate(zip(transitions, responses, metrics, strict=True)):
        ratio, direction = largest_ratio(transition[:n, :n], m0, metric)
        ratios.append(ratio)
        if ratio > worst:
            worst, worst_direction, worst_index = ratio, direction.copy(), i
        forced.append(float(response[:n] @ metric @ response[:n]))
        bias = transition[:n, 18:]
        bias_gains.append(float(np.linalg.eigvalsh(G.sym(bias.T @ metric @ bias))[-1]))
    endpoint = witness(steps, m0, metrics, direction, len(steps)-1, measurements)
    prefix = witness(steps, m0, metrics, worst_direction, worst_index, measurements)
    for expected, checked in ((ratios[-1], endpoint), (worst, prefix)):
        if abs(float(checked["ratio_80_digit"])-expected) > 1e-7*max(1., expected):
            raise ValueError("80-digit complete-word direction check disagrees")
    report = {"dimension": n, "endpoint_ratio": ratios[-1],
              "endpoint_distance_to_one": 1-ratios[-1], "prefix_ratio_max": worst,
              "limiting_prefix": {k: steps[worst_index][k] for k in ("index", "stage", "kind")},
              "endpoint_witness": endpoint, "prefix_witness": prefix,
              "particular_response_energy_endpoint": forced[-1],
              "particular_response_energy_max": max(forced),
              "particular_response_endpoint_full21": responses[-1].tolist(),
              "P4_PASS": False}
    if n == 18:
        report["separate_budget_inequality"] = (
            "sqrt(W_j) <= sqrt(rho_j)*sqrt(W0_motion) + "
            "sqrt(bias_gain_j)*|b0/(1 m/s^2)| + sqrt(response_energy_j)*|alpha|; "
            "alpha scales the entire same-history forcing template")
        report["initial_bias_gain_squared_endpoint"] = bias_gains[-1]
        report["initial_bias_gain_squared_max"] = max(bias_gains)
        report["zero_bias_restriction_only_in_homogeneous_motion_ratio"] = True
        report["bias_columns_retained_in_separate_budget"] = True
        radius = 2*np.tan(np.deg2rad(30.)/2)
        chart = min(radius**2/np.linalg.eigvalsh(np.linalg.solve(m, np.eye(n))[:3, :3])[-1]
                    for m in metrics)
        # Zero true bias is part of this attached capture, hence the closed
        # .4 estimate ball supplies this diagnostic initial-error budget.
        # No such substitution is asserted for other true-bias histories.
        budget = max((.4*np.sqrt(b)+np.sqrt(f))**2
                     for b, f in zip(bias_gains, forced, strict=True))
        report["point_30_degree_chart_storage_min"] = float(chart)
        report["point_zero_true_bias_ball_plus_template_prefix_bound"] = float(budget)
        report["sufficient_prefix_bound_over_point_chart"] = float(budget/chart)
    return report


def transported_information(steps, transitions, metrics, measurements, m0, centered):
    """Whole-word information minus its signed finite-map storage remainder."""
    n = len(m0)
    channels = {kind: np.zeros((n, n)) for kind in
                ("accelerometer", "magnetometer", "S_zero")}
    before = np.eye(21)
    for transition, measurement in zip(transitions, measurements, strict=True):
        if measurement is not None:
            h = measurement["H"] @ before[:, :n]
            channels[measurement["kind"]] += h.T @ np.linalg.solve(measurement["R"], h)
        before = transition
    root = whitener(m0)
    normalized = {k: G.sym(root.T @ value @ root) for k, value in channels.items()}
    vector = normalized["accelerometer"] + normalized["magnetometer"]
    total = vector + normalized["S_zero"]
    terminal = transitions[-1][:n, :n]
    loss = G.sym(root.T @ (m0-terminal.T @ metrics[-1] @ terminal) @ root)
    remainder = G.sym(loss-total)
    values, directions = np.linalg.eigh(vector)
    eta6 = G.sym(channels["accelerometer"][:6, :6] + channels["magnetometer"][:6, :6])
    eta_values, eta_directions = np.linalg.eigh(eta6)
    all_values = np.linalg.eigvalsh(total)
    minimum_loss = float(np.linalg.eigvalsh(loss)[0])
    if abs(minimum_loss-centered["endpoint_distance_to_one"]) > 1e-7:
        raise ValueError("transported ledger detached from controlling endpoint")
    # A scalarization control demonstrates why separately bounding these two
    # matrices can destroy the positive signed complete-word margin.
    scalar_lower = float(all_values[0] + np.linalg.eigvalsh(remainder)[0])
    return {"transport": "every full corrected finite coefficient, before each executed measurement",
            "vector_information_normalized_eigenvalues": values.tolist(),
            "vector_information_min_direction_full21": np.r_[root @ directions[:, 0], np.zeros(21-n)].tolist(),
            "eta6_SI_information_eigenvalues": eta_values.tolist(),
            "eta6_SI_min_direction": eta_directions[:, 0].tolist(),
            "vector_plus_actual_S_information_normalized_eigenvalues": all_values.tolist(),
            "signed_complete_word_loss_eigenvalues": np.linalg.eigvalsh(loss).tolist(),
            "signed_remainder_eigenvalues": np.linalg.eigvalsh(remainder).tolist(),
            "separately_scalarized_loss_lower": scalar_lower,
            "matrix_identity_max_defect": float(np.max(np.abs(loss-total-remainder))),
            "endpoint_ratio": centered["endpoint_ratio"],
            "endpoint_distance_to_one": minimum_loss,
            "endpoint_witness": centered["endpoint_witness"],
            "prefix_ratio_max": centered["prefix_ratio_max"],
            "prefix_witness": centered["prefix_witness"],
            "complete_signed_matrix_required": True,
            "point_Gramian_is_not_uniform_PE_or_P3_admission": True,
            "P3_THRESHOLD_UNCHANGED": 1e-18, "P4_PASS": False}


def audit_mode(root, rows, points, mode):
    steps, _, _, defects, counts = G.build_word(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite factorization failed: "+repr(defects.failures[:1]))
    events = [r for r in rows if r["word"] == mode]
    point = next(r for r in points if r["word"] == mode)
    _, covariance = C.SOURCE.error_and_covariance(point, 21)
    full_m0 = G.sym(np.linalg.solve(covariance, np.eye(21)))
    transitions, responses, full_metrics, measurements = compose(steps, events)
    # Motion storage in BOTH modes; full active-error storage is an additional
    # diagnostic and cannot replace the bounded-bias BRMM motion theorem.
    m0 = full_m0[:18, :18]
    metrics = [m[:18, :18] for m in full_metrics]
    centered = evaluate(steps, transitions, responses, m0, metrics, measurements)
    result = {"counts": dict(counts), "factorization_defects": dict(defects.defects),
              "source_time_start": point["source_time"],
              "source_time_end": events[-1]["source_time"],
              "source_centered_motion": centered,
              "transported_vector_information": transported_information(
                  steps, transitions, metrics, measurements, m0, centered),
              "physical_storages": {}}
    if mode == "A21":
        result["source_centered_full_active_error"] = evaluate(
            steps, transitions, responses, full_m0, full_metrics, measurements)
    scales = {"SI": np.ones(21), "gravity_word_units": physical_scales(root["gravity"], 3.)}
    for name, units in scales.items():
        physical = np.diag(1/units[:18]**2)
        item = evaluate(steps, transitions, responses, physical, [physical]*len(steps))
        item["coordinate_units_full21"] = units.tolist()
        item["metric_selected_without_replay_optimization"] = True
        result["physical_storages"][name] = item
    # Congruent physical-coordinate changes of the original metric must leave
    # the complete-word ratio invariant; they cannot manufacture contraction.
    units = np.diag(scales["gravity_word_units"][:18])
    transformed = np.linalg.solve(units, transitions[-1][:18, :18] @ units)
    ratio, _ = largest_ratio(transformed, units @ m0 @ units, units @ metrics[-1] @ units)
    defect = abs(ratio-centered["endpoint_ratio"])
    if defect > 1e-7:
        raise ValueError("coordinate invariance control failed")
    result["physical_coordinate_invariance_defect"] = defect
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--attachment", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    attachment = json.loads(args.attachment.read_text())
    paths = {suffix: Path(str(args.prefix)+suffix) for suffix in
             (".root.json", ".inputs.csv", ".prefixes.jsonl", ".events.jsonl")}
    hashes = {suffix: hashlib.sha256(path.read_bytes()).hexdigest() for suffix, path in paths.items()}
    if not attachment["read_only_trace_recovers_baseline_bit_for_bit"]:
        raise ValueError("passive source trace parity is required")
    if hashes[".events.jsonl"] != attachment["event_trace_sha256"]:
        raise ValueError("detached event trace")
    if any(hashes[k] != attachment["baseline_capture_sha256"][k] for k in
           (".root.json", ".inputs.csv", ".prefixes.jsonl")):
        raise ValueError("detached source/root/prefix")
    rows = [json.loads(line) for line in paths[".events.jsonl"].read_text().splitlines()]
    points = [json.loads(line) for line in paths[".prefixes.jsonl"].read_text().splitlines()]
    root = json.loads(paths[".root.json"].read_text())
    report = {"experiment": "THREE_CONNECTED_COMPLETE_WORD_STORAGE_ROUTES",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "frozen_coefficients_only": True, "same_history_forcing_preserved": True,
              "source_uniform_response_bound_established": False,
              "nonlinear_coefficient_cover_established": False,
              "physical_source_admission_pass": False, "P4_PASS": False, "P5_MAY_START": False,
              "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached word: "+mode)
        report["modes"][mode] = audit_mode(root, rows, points, mode)
        print("COMPLETE_WORD_ROUTES", mode, "completed", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
