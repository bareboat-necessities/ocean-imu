#!/usr/bin/env python3
"""Non-promoting prescribed-storage test on unreseeded native quiet histories.

The source is p=v=a=0, identity physical attitude, zero BIAS0 and constant
horizontal magnetic field. The native probe rejects nonzero nominal states or
innovations. Consequently the zero-residual tangent uses the actual prediction
factors and I-KH, including all 21 covariance coordinates and active-bias outputs.
It is NOT a finite-error word or a universal source certificate. No metric is
fitted: every endpoint uses M(P)=diag(P21^-1,I3). Uniform coercivity is unproved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import mpmath as mp
import numpy as np

from tools.stability.ou3_alt_contraction import magnetic_service_formulation as FORM
from tools.stability.ou3_alt_contraction import proof_plan as PLAN
from tools.stability.ou3_alt_contraction import shipping_finite_identity as NATIVE

QUALIFICATION = "OU3_ALT_NATIVE_CARRIED_COVARIANCE_STORAGE_DIAGNOSTIC_V1"
SOURCE = NATIVE.ROOT / "tests/ou3_alt_contraction/carried_storage_probe.cpp"
WINDOW = 600


def covariance_metric(P):
    """One fixed storage law, never an eigenvalue floor or wordwise fit."""
    P = np.asarray(P, dtype=float)
    if P.shape != (21, 21) or not np.isfinite(P).all():
        raise ValueError("full finite 21-state covariance required")
    if not np.array_equal(P, P.T):
        raise ValueError("native covariance is not symmetric")
    np.linalg.cholesky(P)
    M = np.eye(24)
    inverse = np.linalg.solve(P, np.eye(21))
    M[:21, :21] = (inverse + inverse.T) / 2
    return M


def lift(A, kind):
    J = np.eye(24)
    J[:21, :21] = A
    if kind == 2:
        # BIAS0: beta_true is constant; e_ba = beta_true - ba_hat.
        J[18:21, 21:24] = np.eye(3) - A[18:21, 18:21]
    return J


def high_precision_ratio(A, P0, P1, digits=60):
    """Re-evaluate storage and terminal eigensolve at arbitrary precision.

    A is independently accumulated in extended precision from exact binary32
    prediction/gain operands. This check does not enclose its roundoff.
    """
    with mp.workdps(digits):
        def matrix(x):
            return mp.matrix([[mp.mpf(str(v)) for v in row] for row in x])
        # Only initial motion columns are unsupplied; ALL 21 output rows stay.
        B = matrix(A[:21, :18])
        M0 = matrix(P0) ** -1
        M1 = matrix(P1) ** -1
        den = M0[:18, :18]
        num = B.T * M1 * B
        L = mp.cholesky(den)
        Li = L ** -1
        W = Li * num * Li.T
        values = mp.eigsy((W + W.T) / 2, eigvals_only=True)
        return mp.nstr(values[values.rows - 1, 0], digits)


def analyze(path, metadata, window_samples=WINDOW, include_map=False):
    if not isinstance(window_samples, int) or window_samples <= 0:
        raise ValueError("positive integer window size required")
    lines = iter(Path(path).read_text().splitlines())

    def parse(line):
        words = line.split()
        return int(words[0]), int(words[1]), np.array(list(map(float, words[2:])))

    step, kind, values = parse(next(lines))
    if (step, kind) != (0, 0) or len(values) != 442:
        raise ValueError("native initial covariance missing")
    P = values[1:].reshape(21, 21)
    mode = bool(values[0])
    P0 = P.copy()
    M0 = covariance_metric(P)
    Phi = np.eye(24)
    Phi_high = np.eye(24, dtype=np.longdouble)
    events = []
    operations = []
    output = []
    counts = {}
    last_boundary = 0
    changes = 0
    pending = None
    for line in lines:
        step, kind, values = parse(line)
        if step != last_boundary + 1:
            raise ValueError("native history is not consecutive")
        if kind in (2, 12, 22, 32):
            size = 441 if kind == 2 else 513
            if len(values) != size or pending is not None:
                raise ValueError("incomplete or detached event")
            A = values[:441].reshape(21, 21)
            if kind == 22:
                H = values[441:504].reshape(3, 21)
                R = values[504:].reshape(3, 3)
                # Actual H/R and ALL intervening updates, not [1,t].
                B = np.linalg.solve(np.linalg.cholesky(R), H @ Phi[:21, [2, 5]])
                events.append(FORM.MagneticEvent(step * metadata["dt_s"], True, True, B))
            J = lift(A, kind)
            Phi = J @ Phi
            Phi_high = J.astype(np.longdouble) @ Phi_high
            pending = (step, kind, J)
            counts[str(kind)] = counts.get(str(kind), 0) + 1
        elif kind in (100, 1002, 1003, 1004, 1011, 1021, 1031, 1061):
            if len(values) != 442 or values[0] not in (0, 1):
                raise ValueError("invalid covariance endpoint")
            next_mode = bool(values[0])
            next_P = values[1:].reshape(21, 21)
            next_M = covariance_metric(next_P)
            if kind in (1002, 1011, 1021, 1031):
                expected = {1002: 2, 1011: 12, 1021: 22, 1031: 32}[kind]
                if pending is None or pending[:2] != (step, expected):
                    raise ValueError("covariance does not complete its literal event")
                J = pending[2]
                pending = None
            else:
                if pending is not None:
                    raise ValueError("uncompleted event before covariance-only operation")
                J = np.eye(24)
            if next_mode != mode:
                if kind != 1061 or mode:
                    raise ValueError("unrepresented mode transition")
                changes += 1
            operations.append((step, kind, J, next_M))
            P, mode = next_P, next_mode
            if kind == 100:
                last_boundary = step
                if step % window_samples == 0:
                    ratio = FORM.projected_storage_ratio(Phi, M0, next_M,
                                                         FORM.joint24_motion_injection())
                    hp = high_precision_ratio(Phi_high, P0, P)
                    x = np.array(ratio["maximizing_direction"])
                    initial_energy = float(x @ M0 @ x)
                    energy = initial_energy
                    consumption = {}
                    for _, op, Jop, M_after in operations:
                        x = Jop @ x
                        after = float(x @ M_after @ x)
                        consumption[str(op)] = consumption.get(str(op), 0.0) + (after - energy) / initial_energy
                        energy = after
                    info = sum((e.transported_whitened_heading_rows.T @
                                e.transported_whitened_heading_rows for e in events),
                               start=np.zeros((2, 2)))
                    start = (step - window_samples) * metadata["dt_s"]
                    end = step * metadata["dt_s"]
                    gaps = np.diff([start, *(e.time_s for e in events), end])
                    output.append({
                        "start_sample": step - window_samples, "end_sample": step,
                        "mode_at_end": "A" if mode else "H", "release_edges": changes,
                        "storage_ratio": ratio, "rho_60_digit_terminal_check": hp,
                        "double_vs_extended_terminal_difference": abs(float(hp) - ratio["rho_point"]),
                        "identity_storage_ratio": FORM.projected_storage_ratio(
                            Phi, np.eye(24), np.eye(24), FORM.joint24_motion_injection())["rho_point"],
                        "accepted_magnetic_events": len(events),
                        "minimum_transported_heading_response": min(
                            (float(np.linalg.norm(e.transported_whitened_heading_rows[:, 0]))
                             for e in events), default=0.0),
                        "transported_heading_bias_gramian": info.tolist(),
                        "transported_heading_bias_min_eigenvalue": float(np.linalg.eigvalsh(info)[0]),
                        "max_observed_magnetic_gap_s": float(max(gaps)),
                        "event_counts": counts,
                        "maximizing_direction_energy_change_by_operation": consumption,
                        "energy_accounting_residual": abs(sum(consumption.values()) - (ratio["rho_point"] - 1)),
                        "covariance_before_eigen_extrema": np.linalg.eigvalsh(P0)[[0, -1]].tolist(),
                        "covariance_after_eigen_extrema": np.linalg.eigvalsh(P)[[0, -1]].tolist(),
                    })
                    if include_map:
                        output[-1]["joint24_tangent_map"] = Phi.tolist()
                        output[-1]["covariance_before"] = P0.tolist()
                        output[-1]["covariance_after"] = P.tolist()
                    # Reset the measurement origin ONLY. Native covariance,
                    # frontend, tuner, clock and source history keep running.
                    P0, M0 = P.copy(), next_M
                    Phi = np.eye(24)
                    Phi_high = np.eye(24, dtype=np.longdouble)
                    events, operations, counts, changes = [], [], {}, 0
        else:
            raise ValueError("unrecognized native operation")
    if pending is not None or last_boundary != metadata["samples"] or last_boundary % window_samples:
        raise ValueError("incomplete native superword")
    return output


def run(work, samples=1800):
    if samples < 1200 or samples % WINDOW:
        raise ValueError("at least two complete 600-sample windows required")
    PLAN.require_theorem_task(
        obligation="prescribed compatible storage on carried informative quiet words",
        evidence_kind="analytic_stationary_source_native_tangent", complete_physical_word=False,
        requested_phase="feasibility_diagnostic")
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    manifest = NATIVE.make_overlay(work / "include")
    binaries = {}
    for instrumented in (True, False):
        binary = work / ("probe" if instrumented else "plain")
        command = [os.environ.get("CXX", "g++"), "-std=c++20", "-O1", "-ffp-contract=off",
                   "-fno-fast-math", "-DEIGEN_DONT_VECTORIZE"]
        if instrumented:
            command.append("-I" + str(work / "include"))
        command += ["-I" + str(NATIVE.eigen_include()), "-I" + str(NATIVE.ROOT / "src"),
                    str(SOURCE), "-o", str(binary)]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=240)
        binaries[instrumented] = binary
    histories = {}
    for mode in ("H", "A", "HA"):
        meta = None
        for instrumented in (True, False):
            label = mode + ("-observed" if instrumented else "-plain")
            proc = subprocess.run([str(binaries[instrumented]), mode, str(samples),
                                   str(work / (label + ".txt")), str(work / (label + ".bin"))],
                                  check=True, capture_output=True, text=True, timeout=120)
            current = json.loads(proc.stdout)
            if meta is not None and current != meta:
                raise ValueError("instrumentation changed startup/runtime metadata")
            meta = current
        if not NATIVE.same_files(work / (mode + "-observed.bin"), work / (mode + "-plain.bin")):
            raise ValueError("instrumentation changed shipping states")
        histories[mode] = {"native": meta, "windows": analyze(work / (mode + "-observed.txt"), meta)}
    worst = max((row["storage_ratio"]["rho_point"], mode, row["end_sample"])
                for mode, history in histories.items() for row in history["windows"])
    report = {
        "qualification": QUALIFICATION, "phase": "feasibility_diagnostic",
        "metric_law": "diag(inverse(full_carried_P21), I3)",
        "source": "stationary p=v=a=S=0, zero BIAS0, identity attitude, B=(32,0,0) uT",
        "source_manifest": manifest, "probe_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "histories": histories, "worst_rho_point": worst[0], "worst_mode": worst[1],
        "worst_window_end_sample": worst[2],
        "passive_instrumentation_bit_identity": True, "same_native_history_across_windows": True,
        "nominal_state_and_innovations_exactly_zero": True, "metric_fitted": False,
        "tangent_only": True, "finite_error_map_certified": False,
        "uniform_metric_coercivity_certified": False, "source_uniform_service_certified": False,
        "source_uniform_rho_certified": False, "interval_enclosure_authorized": False,
        "storage_search_allowed": False, "ALT_STARTUP_PASS": False,
        "ALT_LIVE_PASS": False, "ALT_END_TO_END_PASS": False,
    }
    PLAN.assert_non_promoting_report(report)
    failures = validate(report)
    if failures:
        raise ValueError("native diagnostic validation failed: " + repr(failures))
    return report


def validate(report):
    failures = []
    if report.get("qualification") != QUALIFICATION:
        failures.append("qualification mismatch")
    for key in ("finite_error_map_certified", "uniform_metric_coercivity_certified",
                "source_uniform_service_certified", "source_uniform_rho_certified",
                "interval_enclosure_authorized", "storage_search_allowed", "ALT_STARTUP_PASS",
                "ALT_LIVE_PASS", "ALT_END_TO_END_PASS", "metric_fitted"):
        if report.get(key) is not False:
            failures.append(key + " not false")
    for key in ("passive_instrumentation_bit_identity", "same_native_history_across_windows",
                "nominal_state_and_innovations_exactly_zero", "tangent_only"):
        if report.get(key) is not True:
            failures.append(key + " not true")
    histories = report.get("histories", {})
    if set(histories) != {"H", "A", "HA"}:
        failures.append("H/A/release family missing")
    for mode, history in histories.items():
        windows = history.get("windows", [])
        if len(windows) < 2:
            failures.append(mode + " consecutive windows missing")
        for i, row in enumerate(windows):
            if (row["start_sample"], row["end_sample"]) != (i * WINDOW, (i + 1) * WINDOW):
                failures.append(mode + " windows are not consecutive")
            expected_mode = "H" if mode == "H" or (mode == "HA" and i == 0) else "A"
            if row["mode_at_end"] != expected_mode:
                failures.append(mode + " native mode mismatch")
            if any(row["event_counts"].get(kind) != WINDOW for kind in ("2", "12")):
                failures.append(mode + " prediction/accelerometer word incomplete")
            if row["accepted_magnetic_events"] != WINDOW // 8:
                failures.append(mode + " accepted native service incomplete")
            if row["transported_heading_bias_min_eigenvalue"] <= 0:
                failures.append(mode + " heading/bias pair uninformative")
            if row["release_edges"] != int(mode == "HA" and i == 1):
                failures.append(mode + " release accounting mismatch")
            if row["energy_accounting_residual"] > 1e-8:
                failures.append(mode + " energy accounting failed")
            if row["double_vs_extended_terminal_difference"] > 1e-8:
                failures.append(mode + " numerical precision disagreement")
            # No requirement that rho pass: a failed mathematical attempt is
            # evidence to analyze, not grounds for fabricating a new metric.
            if not np.isfinite(row["storage_ratio"]["rho_point"]):
                failures.append(mode + " nonfinite ratio")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=1800)
    args = parser.parse_args()
    report = run(args.work, args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: report[k] for k in ("worst_rho_point", "worst_mode",
                                          "source_uniform_rho_certified")}))


if __name__ == "__main__":
    main()
