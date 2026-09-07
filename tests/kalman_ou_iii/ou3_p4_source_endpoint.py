"""Audit an actual source-connected finite endpoint without hiding forcing.

The float shipping recursion owns P, resets, nominal state and tuner history.
High precision here verifies endpoint storage evaluation, NOT a high-precision
reimplementation of the observer and NOT the unforced nonlinear P4 master.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).with_name("fixtures")


def stokes(root, t):
    atoms = np.asarray(root["atoms"], dtype=float)
    omega, k, a, phase, dx, dy = atoms.T
    coeff = np.stack((a, k*a*a/2, 3*k*k*a**3/8))
    n = np.arange(1, 4)[:, None]
    frequency = n*omega
    theta = n*(phase - omega*t)
    directions = np.stack((dx, dy))
    disp = np.r_[-directions @ (coeff*np.cos(theta)).sum(axis=0),
                 (coeff*np.sin(theta)).sum()]
    velocity = np.r_[-directions @ (frequency*coeff*np.sin(theta)).sum(axis=0),
                     -(frequency*coeff*np.cos(theta)).sum()]
    acc = np.r_[directions @ (frequency**2*coeff*np.cos(theta)).sum(axis=0),
                -(frequency**2*coeff*np.sin(theta)).sum()]
    integral = np.r_[directions @ (coeff*np.sin(theta)/frequency).sum(axis=0),
                     (coeff*np.cos(theta)/frequency).sum()]
    displaced_theta = n*(k*(dx*disp[0]+dy*disp[1]) - omega*t + phase)
    sx, sy = directions @ (n*k*coeff*np.cos(displaced_theta)).sum(axis=0)
    pitch = math.atan(-sx)
    roll = math.atan2(sy, math.sqrt(1+sx*sx))
    cp, sp, cr, sr = math.cos(pitch), math.sin(pitch), math.cos(roll), math.sin(roll)
    rotation = np.array([[cp, 0, -sp], [sp*sr, cr, cp*sr], [sp*cr, -sr, cp*cr]])
    basis = np.array([[0., 1, 0], [1, 0, 0], [0, 0, -1]])
    return (np.concatenate([basis @ x for x in (velocity, disp, integral, acc)]),
            basis @ rotation @ basis.T)


def error_and_covariance(row, dimension):
    rt = np.asarray(row["R_true"]).reshape(3, 3)
    rh = np.asarray(row["R_hat"]).reshape(3, 3)
    rotation = rt @ rh.T
    # Exact rational Cayley chart formula; no small-angle surrogate.
    chart = 2*(rotation - np.eye(3)) @ np.linalg.inv(rotation + np.eye(3))
    cayley = np.array([chart[2, 1]-chart[1, 2], chart[0, 2]-chart[2, 0],
                       chart[1, 0]-chart[0, 1]])/2
    estimated = np.asarray(row["x_hat"])
    error = -estimated.copy()
    error[:3] = cayley
    error[6:18] += row["linear_true"]
    error[18:21] += row["true_bias"]
    p = np.asarray(row["P"]).reshape(21, 21)[:dimension, :dimension]
    p = (p+p.T)/2
    np.linalg.cholesky(p)
    return error[:dimension], p


def rational_solve(matrix, rhs):
    """Exact elimination for endpoint diagnostics, with explicit pivoting."""
    a = [list(map(Fraction, values))+[Fraction(value)]
         for values, value in zip(matrix, rhs)]
    n = len(a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError("singular exact endpoint matrix")
        a[j], a[pivot] = a[pivot], a[j]
        for i in range(j+1, n):
            scale = a[i][j]/a[j][j]
            for k in range(j+1, n+1):
                a[i][k] -= scale*a[j][k]
            a[i][j] = Fraction(0)
    x = [Fraction(0)]*n
    for j in reversed(range(n)):
        x[j] = (a[j][n]-sum(a[j][k]*x[k] for k in range(j+1, n)))/a[j][j]
    return x


def decimal_text(value, digits=80):
    value = Fraction(value)
    with localcontext() as ctx:
        ctx.prec = digits
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def high_precision_energy(row, dimension, digits=80):
    # Recompute the chart and inverse in high precision from the captured
    # floating-point coordinates, not merely the final double energy ratio.
    rt = np.asarray([Fraction(v) for v in row["R_true"]], dtype=object).reshape(3, 3)
    rh = np.asarray([Fraction(v) for v in row["R_hat"]], dtype=object).reshape(3, 3)
    rotation = rt@rh.T
    identity = np.eye(3, dtype=int)
    chart = np.asarray([rational_solve((rotation+identity).T, r)
                        for r in 2*(rotation-identity)], dtype=object)
    x = [Fraction(-v) for v in row["x_hat"][:dimension]]
    for i, (j, k) in enumerate(((2, 1), (0, 2), (1, 0))):
        x[i] = (chart[j, k]-chart[k, j])/2
    for i, value in enumerate(row["linear_true"], 6):
        x[i] += Fraction(value)
    if dimension == 21:
        for i, value in enumerate(row["true_bias"], 18):
            x[i] += Fraction(value)
    p = [[(Fraction(row["P"][21*i+j])+Fraction(row["P"][21*j+i]))/2
          for j in range(dimension)] for i in range(dimension)]
    value = sum(u*v for u, v in zip(x, rational_solve(p, x)))
    if value <= 0:
        raise ValueError("nonpositive endpoint storage")
    return decimal_text(value, digits)


def decision(*, source_matches, common_premises, forcing_zero):
    if not source_matches:
        return "SOURCE_ATTACHMENT_REJECTED"
    if not common_premises:
        return "FIXED_WORD_OUTSIDE_CHECKED_NORMAL_LIVE_PREMISES"
    if not forcing_zero:
        return "ACTUAL_FORCED_ENDPOINT_NOT_HOMOGENEOUS_P4_TEST"
    return "SUBEVENT_GRAPH_AND_HOMOGENEOUS_ATTACHMENT_STILL_REQUIRED"


def audit(root, rows):
    normal = json.loads((ROOT / "tools/stability/ou3_proof_operating_domain.json").read_text())["normal_live"]
    atoms = np.asarray(root["atoms"])
    if atoms.shape != (128, 6) or not np.isfinite(atoms).all():
        raise ValueError("missing complete finite source root")
    omega, k, a, _phase, dx, dy = atoms.T
    if (not np.allclose(k, omega**2/root["gravity"], rtol=2e-14, atol=0)
            or np.any(k*a > 0.2) or np.any(a < 0)
            or not np.allclose(dx*dx+dy*dy, 1, rtol=2e-14, atol=2e-14)):
        raise ValueError("invalid Stokes dependency graph")
    modes = {}
    for mode, dimension in (("H18", 18), ("A21", 21)):
        word = [r for r in rows if r["word"] == mode]
        if not word or word[0]["event"] != "root":
            raise ValueError("missing fixed word entrance")
        imu = [r for r in word if r["event"] == "imu"]
        if len(imu) != 600 or [r["index"] for r in imu] != list(range(imu[0]["index"], imu[0]["index"]+600)):
            raise ValueError("incomplete 600-sample word")
        energies, source_defects, directions = [], [], []
        for row in word:
            linear, rotation = stokes(root, row["source_time"])
            heel = row["wind_heel"]
            unheel = np.array([[1, 0, 0], [0, math.cos(heel), math.sin(heel)],
                              [0, -math.sin(heel), math.cos(heel)]])
            # Float unheel in the shipping implementation is measured below;
            # this numerical tolerance is a point parity check, not outward.
            source_defects.append(max(float(np.max(np.abs(linear-row["linear_true"]))),
                                      float(np.max(np.abs(unheel@rotation-np.asarray(row["R_true"]).reshape(3, 3))))))
            e, p = error_and_covariance(row, dimension)
            energies.append(float(e @ np.linalg.solve(p, e)))
            directions.append(e)
        v0, vn = high_precision_energy(word[0], dimension), high_precision_energy(word[-1], dimension)
        rho = Fraction(vn)/Fraction(v0)
        rho_text, distance = decimal_text(rho, 65), decimal_text(1-rho, 65)
        force = np.asarray([r["prediction_forcing"] for r in imu])
        norms = {
            "physical_acceleration": max(np.linalg.norm(r["physical_acceleration"]) for r in word),
            "physical_body_rate_deg_s": max(np.linalg.norm(r["physical_gyro"]) for r in word)*180/math.pi,
            "estimated_bias": max(np.linalg.norm(r["x_hat"][18:21]) for r in word),
            "attitude_error_deg": max(2*math.atan(np.linalg.norm(e[:3])/2) for e in directions)*180/math.pi,
        }
        checks = {
            "all_Live": all(r["live"] for r in word),
            "same_requested_mode": all(r["active"] == (mode == "A21") for r in word),
            "no_magnetic_lock_or_refinement_change": len({(r["mag_lock"], r["mag_refined"]) for r in word}) == 1,
            "all_accelerometers_accepted": all(r["acc_accepted"] for r in imu),
            "physical_acceleration_cap": norms["physical_acceleration"] <= normal["non_gravitational_cog_acceleration_norm_upper_mps2"],
            "physical_body_rate_cap": norms["physical_body_rate_deg_s"] <= normal["body_rate_norm_upper_deg_s"],
            "A21_bias_interior": mode == "H18" or norms["estimated_bias"] <= normal["active_accelerometer_bias_state_norm_upper_mps2"],
            "finite_Cayley_45_degree_chart": norms["attitude_error_deg"] <= 45,
            "BIAS1_zero_true_root_and_driver": all(r["true_bias"] == [0, 0, 0] and r["tau_b"] == 5000 for r in word),
        }
        checks = {name: bool(value) for name, value in checks.items()}
        source_matches = max(source_defects) < 2e-7
        forcing_zero = bool(np.max(np.abs(force)) <= 1e-13)
        margins = {event: float(sum(energies[i]-energies[i-1]
                    for i in range(1, len(word)) if word[i]["event"] == event))
                   for event in ("imu", "mag")}
        modes[mode] = {
            "decision": decision(source_matches=source_matches,
                                 common_premises=all(checks.values()), forcing_zero=forcing_zero),
            "fixed_first_sample_index": imu[0]["index"],
            "fixed_last_sample_index": imu[-1]["index"],
            "source_time_root": word[0]["source_time"],
            "source_time_endpoint": word[-1]["source_time"],
            "source_coordinate_max_parity_defect": max(source_defects),
            "source_root_point_parity_pass": source_matches,
            "checked_common_premises": checks,
            "full_Normal_Live_membership_certified": False,
            "maxima": {key: float(value) for key, value in norms.items()},
            "events": {"imu": len(imu), "mag_calls": len(word)-len(imu)-1,
                       "due_S_calls": sum(r["s_due"] for r in imu)},
            "V0_80_digit": v0, "VN_80_digit": vn,
            "actual_forced_endpoint_ratio": rho_text,
            "one_minus_actual_forced_ratio": distance,
            "actual_sample_prefix_ratio_max": max(energies)/energies[0],
            "initial_error_not_optimized": directions[0].tolist(),
            "margin_by_shipping_call": margins,
            "margin_telescope_residual": sum(margins.values())-(energies[-1]-energies[0]),
            "truth_minus_shipping_FLL_prediction_max_abs_by_coordinate": np.max(np.abs(force), axis=0).tolist(),
            "homogeneous_prediction_forcing_zero": forcing_zero,
            "BIAS2_uniform_sector_available": False,
            "homogeneous_P4_endpoint_ratio": None,
            "canonical_P4_falsified": False,
        }
    return {
        "experiment": "SINGLE_SHIPPING_OBSERVER_SOURCE_CONNECTED_ACTUAL_FINITE_ENDPOINT",
        "branch": root["branch"], "new_source_capture_not_original_retained_payload": True,
        "original_retained_witness_modified_or_reselected": False,
        "initial_covariance_from_actual_shipping_history": True,
        "nominal_and_covariance_reset_feedback_owned_by_shipping": True,
        "physical_noise_driver": "ZERO_BUT_WAVE_MODEL_AND_GYRO_DISCRETIZATION_FORCING_REMAIN",
        "hardware_qualification_required_for_conditional_run": False,
        "subevent_nonlinear_graph_captured": False,
        "modes": modes, "P4_promoted": False, "P5_may_start": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--generator-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    provenance = json.loads((FIXTURES / "ou3_physical_generator_provenance.json").read_text())
    commit = subprocess.check_output(["git", "-C", str(args.generator_root), "rev-parse", "HEAD"], text=True).strip()
    if commit != provenance["generator_commit"]:
        raise ValueError("generator commit mismatch")
    for name, digest in provenance["generator_file_sha256"].items():
        if hashlib.sha256((args.generator_root/name).read_bytes()).hexdigest() != digest:
            raise ValueError("generator source digest mismatch: " + name)
    paths = {suffix: Path(str(args.prefix)+suffix) for suffix in (".root.json", ".inputs.csv", ".prefixes.jsonl")}
    root = json.loads(paths[".root.json"].read_text())
    rows = [json.loads(line) for line in paths[".prefixes.jsonl"].read_text().splitlines()]
    report = audit(root, rows)
    report["capture_sha256"] = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}
    report["generator_commit"] = commit
    report["observer_repository_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    for mode, result in report["modes"].items():
        print("SOURCE_CONNECTED_ENDPOINT", mode, json.dumps(result, allow_nan=False))


if __name__ == "__main__":
    main()
