"""Connected finite error/source/forcing attachment before the motion master.

Uses a read-only trace of ONE shipping observer, never a reset-deleted or
reseeded replay. Verifies every prediction, Joseph injection, finite Cayley
reset and bias projection. S=0 contributes -S_true even with sensor noise off.
This numerical attachment test does not invent graph sectors, fit gains to a
trajectory, or turn a forced energy ratio into a contraction certificate.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_source_endpoint as SOURCE

ROUNDING_PARITY_TOL = 128*np.finfo(np.float32).eps


class ParityChecks:
    """Explicit per-word/current-event context, never a loop closure."""
    def __init__(self):
        self.defects = defaultdict(float)
        self.failures = []
        self.row = None

    def __call__(self, name, got, expected, scale=None):
        got, expected = np.asarray(got), np.asarray(expected)
        normalizer = max(float(np.max(np.abs(got))), float(np.max(np.abs(expected))),
                         1e-30, 0. if scale is None else float(scale))
        defect = float(np.max(np.abs(got-expected)))/normalizer
        self.defects[name] = max(self.defects[name], defect)
        if not np.isfinite(defect) or defect > ROUNDING_PARITY_TOL:
            self.failures.append({"check": name, "index": self.row["index"], "stage": self.row["stage"],
                                  "normalized_defect": defect, "normalizer": normalizer})


def mat(row, key, n=3, m=None):
    return np.asarray(row[key], dtype=float).reshape(n, n if m is None else m)


def skew(x):
    a, b, c = x
    return np.array([[0., -c, b], [c, 0., -a], [-b, a, 0.]])


def rotation(c):
    s = skew(c)
    return np.eye(3)+(s+.5*s@s)/(1+np.dot(c, c)/4)


def cayley(r):
    s = 2*(r-np.eye(3)) @ np.linalg.inv(r+np.eye(3))
    return np.array([s[2, 1]-s[1, 2], s[0, 2]-s[2, 0], s[1, 0]-s[0, 1]])/2


def deployed_rotation(d):
    # Literal normalized polynomial/trigonometric quat_from_delta_theta
    # branches; the trace owns the actual float/FMA arithmetic.
    angle = np.linalg.norm(d)
    if angle < 1e-2:
        t2, t4 = angle**2, angle**4
        w, k = 1-t2/8+t4/384, .5-t2/48+t4/3840
    else:
        w, k = np.cos(angle/2), np.sin(angle/2)/angle
    return rotation(2*k*np.asarray(d)/w)


def residual_graph(row, error):
    """Full physical finite residual = homogeneous graph + explicit source.

    True source outputs and the actual frontend-delivered measurement remain
    separate. No 'measurement noise = 0' assumption removes these defects.
    """
    rh, rt = mat(row, "R_hat"), mat(row, "R_true")
    estimate = np.asarray(row["x_hat"])
    true = np.asarray(row["linear_true"])
    measured = np.asarray(row["measured"])
    er = rotation(error[:3])
    h = np.zeros((3, 21))
    if row["kind"] == "accelerometer":
        gravity = np.array([0., 0., row["gravity"]])
        f = rh @ (estimate[15:18]-gravity)
        h[:, :3], h[:, 15:18], h[:, 18:21] = -skew(f), rh, np.eye(3)
        homogeneous = (er-np.eye(3))@f + er@rh@error[15:18] + error[18:21]
        source = measured-rt@(true[9:12]-gravity)-row["true_bias"]-row["temperature_mean"]
    elif row["kind"] == "magnetometer":
        reference = np.asarray(row["mag_reference"])
        f = rh@reference
        h[:, :3] = -skew(f)
        homogeneous = (er-np.eye(3))@f
        source = measured-rt@reference
    elif row["kind"] == "S_zero":
        h[:, 12:15] = np.eye(3)
        homogeneous = error[12:15]
        source = -true[6:9]
    else:
        raise ValueError("unknown connected measurement")
    return homogeneous+source, source, h


def correct(error, gain, residual, true_bias, radius):
    correction = gain@residual
    out = error-correction
    out[:3] = cayley(rotation(error[:3]) @ deployed_rotation(correction[:3]).T)
    unprojected = np.asarray(true_bias)-out[18:21]
    length = np.linalg.norm(unprojected)
    projected = unprojected*min(1., radius/length) if length else unprojected
    out[18:21] = np.asarray(true_bias)-projected
    return out


def motion_energy(row):
    """18-error storage using the principal block of the FULL information.

    Solving the full 21-state covariance is only storage evaluation: the gain
    and filter are never replaced by an 18-state recursion.
    """
    error, p = SOURCE.error_and_covariance(row, 21)
    error[18:21] = 0
    return float(error @ np.linalg.solve(p, error))


def exact_motion_energy(row):
    # Zero only the bias PERFORMANCE coordinate, not any state in the trace.
    copy = dict(row, true_bias=list(row["x_hat"][18:21]))
    return SOURCE.high_precision_energy(copy, 21)


def audit(root, rows, checkpoints):
    dt = root["dt"]
    result = {}
    for mode in ("H18", "A21"):
        events = [r for r in rows if r["word"] == mode]
        points = [r for r in checkpoints if r["word"] == mode]
        if not points or points[0]["event"] != "root":
            raise ValueError("missing baseline word root")
        if not events or events[0]["stage"] != "prediction_enter":
            raise ValueError("missing prediction entrance")
        check = ParityChecks()
        defects, failures, counts = check.defects, check.failures, Counter()
        forcing = defaultdict(float)
        energy_by_stage = defaultdict(float)
        previous_committed = points[0]
        previous_energy = motion_energy(previous_committed)
        energies = [previous_energy]
        initial_bias_energy = 0.
        pending = None
        predicted = None
        group = None
        last_index = None

        for row in events:
            check.row = row
            stage = row["stage"]
            counts[stage] += 1
            try:
                e, p = SOURCE.error_and_covariance(row, 21)
            except (ValueError, np.linalg.LinAlgError) as exc:
                raise ValueError(f"{mode} sample {row['index']} {stage}: {exc}") from exc
            linear, rt = SOURCE.stokes(root, row["source_time"])
            heel = row["wind_heel"]
            unheel = np.array([[1., 0, 0], [0, np.cos(heel), np.sin(heel)],
                              [0, -np.sin(heel), np.cos(heel)]])
            check("source_linear", row["linear_true"], linear, 1)
            check("source_rotation", mat(row, "R_true"), unheel@rt, 1)
            true_bias = np.asarray(row["true_bias"], dtype=float)
            if row["active"] != (mode == "A21") or true_bias.shape != (3,) or not np.isfinite(true_bias).all():
                raise ValueError("mode or finite BIAS1 physical-bias history detached")

            if stage == "prediction_enter":
                if pending is not None or (last_index is not None and row["index"] != last_index+1):
                    raise ValueError("incomplete/reordered previous sample")
                if group is not None and not group.get("accelerometer_done"):
                    raise ValueError("missing completed accelerometer event")
                # Tuner commits may change P between samples, but not x/q.
                check("inter_sample_nominal", row["x_hat"], previous_committed["x_hat"], 1)
                check("inter_sample_attitude", mat(row, "R_hat"), mat(previous_committed, "R_hat"), 1)
                group = {"entrance": row, "accelerometer_done": False}
                last_index = row["index"]
                initial_bias_energy += dt*float(e[18:21]@e[18:21])
            elif stage == "prediction":
                if group is None or predicted is not None:
                    raise ValueError("prediction order mismatch")
                before = group["entrance"]
                eb, pb = SOURCE.error_and_covariance(before, 21)
                fll, faa = mat(row, "F_LL", 12), mat(row, "F_AA", 6)
                phi = np.exp(-dt/row["tau_b"]) if row["active"] else 1.
                f = np.zeros((21, 21))
                f[:6, :6], f[6:18, 6:18], f[18:, 18:] = faa, fll, phi*np.eye(3)
                q = np.zeros((21, 21))
                q[:6, :6], q[6:18, 6:18] = mat(row, "Q_AA", 6), mat(row, "Q_LL", 12)
                if row["active"]:
                    q[18:, 18:] = -.5*row["tau_b"]*np.expm1(-2*dt/row["tau_b"])*mat(row, "Q_b")
                check("prediction_covariance", p, f@pb@f.T+q,
                      np.linalg.norm(f, np.inf)**2*np.linalg.norm(pb, np.inf)+np.linalg.norm(q, np.inf))
                u = np.asarray(row["linear_true"])-fll@before["linear_true"]
                expected = eb.copy()
                expected[6:18] = fll@eb[6:18]+u
                # e_b = beta_true-b_hat.  The shipping estimate predicts
                # b_hat+ = phi*b_hat, while the physical BIAS1 history keeps
                # its own same-source beta transition.  Retain the resulting
                # exact forcing instead of assuming beta_true==0.
                beta_before = np.asarray(before["true_bias"], dtype=float)
                beta_after = np.asarray(row["true_bias"], dtype=float)
                expected[18:] = phi*eb[18:] + beta_after - phi*beta_before
                unbiased = np.asarray(row["gyro_measured"])
                source_step_defect = mat(row, "R_true") @ (
                    deployed_rotation(-dt*unbiased)@mat(before, "R_true")).T
                rn = deployed_rotation(-dt*np.asarray(row["omega_hat"]))
                rs = deployed_rotation(-dt*(np.asarray(row["omega_hat"])-eb[3:6]))
                expected[:3] = cayley(source_step_defect@rs@rotation(eb[:3])@rn.T)
                check("prediction_finite_error", e, expected, 1+np.linalg.norm(eb, np.inf))
                check("prediction_nominal_rotation", mat(row, "R_hat"), rn@mat(before, "R_hat"), 1)
                forcing["latent_increment_squared_SI_coordinate_norm"] += float(u@u)
                angle = cayley(source_step_defect)
                forcing["gyro_discretization_Cayley_squared_rad"] += float(angle@angle)
                predicted = row
            elif stage == "aw_floor":
                if predicted is None:
                    raise ValueError("floor detached from prediction")
                check("aw_floor_preserves_mean", row["x_hat"], predicted["x_hat"], 1)
                check("aw_floor_preserves_attitude", mat(row, "R_hat"), mat(predicted, "R_hat"), 1)
                predicted = None
            elif stage == "measurement":
                if group is None or pending is not None or predicted is not None:
                    raise ValueError("measurement detached from completed prediction/event")
                for key, n in (("x_hat", None), ("R_hat", 3), ("P", 21)):
                    current = row[key] if n is None else mat(row, key, n)
                    previous = previous_committed[key] if n is None else mat(previous_committed, key, n)
                    check("same_history_measurement_"+key, current, previous)
                kind = row["kind"]
                if kind == "S_zero" and group["accelerometer_done"]:
                    raise ValueError("S event moved after accelerometer")
                if kind == "magnetometer" and not group["accelerometer_done"]:
                    raise ValueError("magnetic event moved before accelerometer")
                residual, source, h = residual_graph(row, e)
                k, pct = mat(row, "K", 21, 3), mat(row, "PCt", 21, 3)
                s, rnoise = mat(row, "innovation_cov"), mat(row, "R")
                scale = 1+np.linalg.norm(row["measured"])+np.linalg.norm(row["linear_true"])
                check("finite_"+kind+"_residual", row["r"], residual, scale)
                gain_h = h.copy()
                if not row["active"]:
                    gain_h[:, 18:21] = 0
                expected_pct = p@gain_h.T
                if not row["active"]:
                    expected_pct[18:21] = 0
                check("actual_"+kind+"_PCt", pct, expected_pct,
                      np.linalg.norm(p, np.inf)*np.linalg.norm(h.T, np.inf))
                check("actual_"+kind+"_innovation_cov", s, h@p@h.T+rnoise,
                      np.linalg.norm(h, np.inf)**2*np.linalg.norm(p, np.inf)+np.linalg.norm(rnoise, np.inf))
                check("actual_"+kind+"_gain_equation", k@s, pct,
                      np.linalg.norm(k, np.inf)*np.linalg.norm(s, np.inf)+np.linalg.norm(pct, np.inf))
                if kind == "S_zero":
                    if row["R"] != row["R_S"]:
                        raise ValueError("S event detached from actual applied R_S")
                forcing[kind+"_source_energy_actual_R_inverse"] += float(source@np.linalg.solve(rnoise, source))
                forcing[kind+"_source_max_norm"] = max(forcing[kind+"_source_max_norm"], float(np.linalg.norm(source)))
                pending = {"row": row, "error": e, "P": p, "K": k, "PCt": pct,
                           "S": s, "expected_residual": residual, "phase": "measurement"}
                counts[kind] += 1
            elif stage == "injected":
                if pending is None or pending["phase"] != "measurement" or row["kind"] != pending["row"]["kind"]:
                    raise ValueError("Joseph injection detached from measurement")
                before, k = pending["row"], pending["K"]
                delta = k@before["r"]
                check("state_injection", row["x_hat"], np.asarray(before["x_hat"])+delta,
                      1+np.linalg.norm(before["x_hat"], np.inf)+np.linalg.norm(delta, np.inf))
                pj = pending["P"]-k@pending["PCt"].T-pending["PCt"]@k.T+k@pending["S"]@k.T
                check("Joseph_full_covariance", p, pj, np.linalg.norm(pending["P"], np.inf))
                check("pre_reset_attitude_unchanged", mat(row, "R_hat"), mat(before, "R_hat"), 1)
                pending.update(injected=row, phase="injected")
            elif stage == "reset":
                if pending is None or pending["phase"] != "injected":
                    raise ValueError("finite reset detached from Joseph injection")
                injected = pending["injected"]
                delta = np.asarray(injected["x_hat"][:3])
                reset = np.eye(21)
                reset[:3, :3] += .5*skew(delta)
                check("shipping_covariance_reset", p, reset@mat(injected, "P", 21)@reset.T,
                      np.linalg.norm(mat(injected, "P", 21), np.inf)*np.linalg.norm(reset, np.inf)**2)
                check("finite_quaternion_injection", mat(row, "R_hat"),
                      deployed_rotation(delta)@mat(injected, "R_hat"), 1)
                check("reset_clears_error_register", row["x_hat"][:3], np.zeros(3), 1)
                pending.update(reset=row, phase="reset")
            elif stage == "projection":
                if pending is None or pending["phase"] != "reset":
                    raise ValueError("projection detached from finite reset")
                before = pending["row"]
                expected = correct(pending["error"], pending["K"], pending["expected_residual"],
                                   row["true_bias"], row["projection_radius"])
                check("complete_nonlinear_correction", e, expected, 1+np.linalg.norm(pending["error"], np.inf))
                check("projection_preserves_covariance", p, mat(pending["reset"], "P", 21))
                check("projection_preserves_other_estimates", row["x_hat"][:18], pending["reset"]["x_hat"][:18], 1)
                bound = np.linalg.norm(row["x_hat"][18:21])
                if bound > row["projection_radius"]*(1+ROUNDING_PARITY_TOL):
                    raise ValueError("committed estimate exceeds projected float ball")
                if before["kind"] == "accelerometer":
                    if group["accelerometer_done"]:
                        raise ValueError("duplicate accelerometer")
                    group["accelerometer_done"] = True
                pending = None
            else:
                raise ValueError("unrecognized shipping event")

            # Includes prediction, floor, each completed measurement, and tuner
            # entry metric changes. Auxiliary injection/reset energies are not
            # mislabelled as complete-event prefix values.
            if stage in ("prediction_enter", "prediction", "aw_floor", "projection"):
                energy = motion_energy(row)
                energy_by_stage[stage+":"+row["kind"]] += energy-previous_energy
                previous_energy = energy
                energies.append(energy)
                previous_committed = row

        if pending is not None or predicted is not None or not group["accelerometer_done"]:
            raise ValueError("incomplete terminal event")
        if counts["prediction"] != 600 or counts["accelerometer"] != 600 or counts["aw_floor"] != 600:
            raise ValueError("not one complete 600-sample shipping word")
        if counts["S_zero"] != sum(r["s_due"] for r in points) or counts["magnetometer"] != sum(r["event"] == "mag" for r in points):
            raise ValueError("trace dropped or inserted a shipping S/magnetic event")
        for key in ("x_hat", "R_hat", "P"):
            if previous_committed[key] != points[-1][key]:
                raise ValueError("terminal state detached from baseline capture: "+key)
        result[mode] = {
            "decision": "CONNECTED_POINT_ATTACHMENT_PASS" if not failures else "CONNECTED_POINT_ATTACHMENT_FAIL",
            "counts": dict(counts), "maximum_normalized_parity_defects": dict(defects),
            "first_failed_attachment": failures[0] if failures else None,
            "failed_attachment_checks": len(failures),
            "motion_storage": "principal 18x18 block of full 21x21 P inverse; actual gains unchanged",
            "W0_80_digit": exact_motion_energy(points[0]),
            "WN_80_digit": exact_motion_energy(points[-1]),
            "actual_forced_motion_ratio": energies[-1]/energies[0],
            "actual_completed_event_prefix_ratio_max": max(energies)/energies[0],
            "signed_motion_energy_by_event": dict(energy_by_stage),
            "signed_energy_telescope_residual": sum(energy_by_stage.values())-(energies[-1]-energies[0]),
            "bias_sample_entrance_energy": initial_bias_energy,
            "source_forcing_channels": dict(forcing),
            "forcing_channels_have_distinct_declared_units_not_summed": True,
            "source_uniform_motion_gains": None, "nonlinear_graph_sector_witness": None,
            "coefficients_or_filter_architecture_falsified": False,
        }
    return {"experiment": "COMPLETE_SOURCE_CONNECTED_SUBEVENT_ERROR_FORCING_ATTACHMENT",
            "numerical_point_parity_tolerance": float(ROUNDING_PARITY_TOL),
            "float_parity_is_an_outward_certificate": False,
            "measurement_S_true_forcing_retained": True,
            "bias_prediction_correction_projection_coupling_retained": True,
            "nonzero_physical_bias_history_supported": True,
            "no_gain_fit_metric_search_or_source_window_reselection": True,
            "modes": result, "P4_MOTION_PASS": False, "P5_MOTION_MAY_START": False,
            "remaining_gate": "source/error joint graph sectors and channel-gain augmented master; then source-uniform prefix/retention cover"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--baseline-prefix", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = json.loads(Path(str(args.prefix)+".root.json").read_text())
    rows = [json.loads(line) for line in Path(str(args.prefix)+".events.jsonl").read_text().splitlines()]
    points = [json.loads(line) for line in Path(str(args.prefix)+".prefixes.jsonl").read_text().splitlines()]
    baseline_paths = {suffix: Path(str(args.baseline_prefix)+suffix) for suffix in
                      (".root.json", ".inputs.csv", ".prefixes.jsonl")}
    traced_paths = {suffix: Path(str(args.prefix)+suffix) for suffix in baseline_paths}
    baseline_hashes = {suffix: hashlib.sha256(path.read_bytes()).hexdigest()
                       for suffix, path in baseline_paths.items()}
    traced_hashes = {suffix: hashlib.sha256(path.read_bytes()).hexdigest()
                     for suffix, path in traced_paths.items()}
    parity = baseline_hashes == traced_hashes
    if not parity:
        raise ValueError("read-only trace changed the shipping source/root/endpoint capture")
    report = audit(root, rows, points)
    report.update({"root_sha256": traced_hashes[".root.json"],
                   "input_trace_sha256": traced_hashes[".inputs.csv"],
                   "prefix_trace_sha256": traced_hashes[".prefixes.jsonl"],
                   "event_trace_sha256": hashlib.sha256(Path(str(args.prefix)+".events.jsonl").read_bytes()).hexdigest(),
                   "baseline_capture_sha256": baseline_hashes,
                   "read_only_trace_recovers_baseline_bit_for_bit": parity,
                   "trace_overlay_manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                   "P4_MOTION_PASS": False, "P5_MOTION_MAY_START": False})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    for mode, value in report["modes"].items():
        print("CONNECTED_MOTION", mode, value["decision"],
              value["actual_forced_motion_ratio"],
              value["actual_completed_event_prefix_ratio_max"], flush=True)


if __name__ == "__main__":
    main()
