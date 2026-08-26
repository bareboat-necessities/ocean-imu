#!/usr/bin/env python3
"""Evaluate the source-correlated P5 sample-1 S/accelerometer prefix.

This is the first stage after the certified sample-1 entry that advances an
actual second measurement packet.  It deliberately keeps three quantities
separate:

* P: the shipping 18x18 covariance enclosure;
* e: the physical non-attitude error enclosure;
* xhat: the estimator's non-attitude mean enclosure.

That separation matters because the S=0 pseudo measurement uses the shipping
residual ``-xhat_S``, not the physical S error.  At sample zero xhat is exactly
zero before the accelerometer; an accepted first accelerometer creates the
linear mean K*r, while the safe-LDLT identity branch leaves it zero.  The next
prediction propagates both branches before sample 1.

Pseudo timing is also source-correlated.  If sample zero was due, the elapsed
clock resets to zero and the same minimum-cadence source is due again at sample
1.  If sample zero was not due, elapsed=h and sample 1 is due only for source
periods no larger than 2h plus the exact float tolerance; otherwise it remains
not due.  Ambiguous source intervals retain both branches.

The sample-1 accelerometer diagnostic uses an exact finite-rotation residual
norm bound

    ||r_a|| <= ||R_e-I|| ||f_pred|| + ||e_aw|| + ||e_ba||,
    ||R_e-I||_2 = 2q/sqrt(4+q^2),

with ||f_pred||<=g+||xhat_aw||.  The shipping H/P/R/K/Joseph/reset calculation
is still evaluated with the full matrix backend.  Numerical nonclosure is a
valid output with a concrete first source witness; this stage never promotes a
complete word and never changes the deployed 6-rad correction primitive.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_post_reset as POST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_entry_v3 as ENTRY3
import ou3_source_reachable_matrix_p3 as P3CELL
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = ENTRY.DEFAULT_DOMAIN
SCHEMA = 1
N = ENTRY.N
FULL = ENTRY.FULL
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0
FLOAT_EPS = 2.0 ** -23


def _zero_state():
    return [FULL.I(0.0) for _ in range(N)]


def _state_hull(a, b):
    return [hull(x, y) for x, y in zip(a, b)]


def _predict_estimator_state(x, F):
    out = list(x)
    # Local attitude error is injected/cleared after every correction.  Gyro
    # bias is constant in the deterministic mean model over one prediction.
    for i in FULL.TH:
        out[i] = FULL.I(0.0)
    for i in FULL.BG:
        out[i] = x[i]
    lin = list(FULL.V)+list(FULL.P)+list(FULL.SS)+list(FULL.AW)
    for i in lin:
        y = FULL.I(0.0)
        for j in lin:
            y = y + F[i][j]*x[j]
        out[i] = y
    return out


def _phase1_options(src: dict, phase0: str) -> list[str]:
    """Source-complete sample-1 pseudo phase conditioned on sample-0 phase."""
    h = float(src["dt_s"])
    tol = 16.0 * FLOAT_EPS
    p = src["pseudo_period_s"]
    if phase0 == "due":
        # Conditioning on sample-0 due implies p<=h+tol.  With the elapsed
        # clock reset to zero, exactly the same condition is tested again.
        return ["due"]
    if phase0 != "not_due":
        raise ValueError("invalid sample-0 pseudo phase")
    # sample-0 not-due conditions p>h+tol and leaves elapsed=h.  At sample 1,
    # total=2h; ambiguous source intervals retain both possibilities.
    threshold2 = 2.0*h + tol
    out = []
    if p.lo <= threshold2:
        out.append("due")
    if p.hi > threshold2:
        out.append("not_due")
    return out or ["not_due"]


def _apply_measurement_state(Pm, e, xhat, H, R, r, *, allow_identity: bool, theta_cap: float | None = None):
    cell = FULL._measurement_cell(Pm, H, R, r)
    dx = list(cell["dx"])
    if theta_cap is not None:
        cap = Interval(-FULL.up(theta_cap), FULL.up(theta_cap))
        for i in range(3):
            dx[i] = FULL._intersect(dx[i], cap)
        # Recompute reset covariance with the certified theta component cap;
        # Joseph itself remains exactly the shipping enclosure.
        PHt, S = FULL._innovation(Pm, H, R)
        Sinv, backend = FULL._spd_inverse_enclosure(S, R)
        K = matrix_mul(PHt, Sinv)
        Pj = FULL._shipping_joseph(Pm, K, S, PHt)
        cell["P_accepted"] = FULL._reset_covariance(Pj, dx[0:3])
        cell["inverse_backend"] = backend
        cell["K"] = K
        cell["S"] = S
        cell["dx"] = dx

    xacc = list(xhat)
    eacc = list(e)
    for i in range(3, N):
        xacc[i] = xhat[i] + dx[i]
        eacc[i] = e[i] - dx[i]
    for i in range(3):
        xacc[i] = FULL.I(0.0)

    if allow_identity:
        Pout = FULL._psd_tighten(FULL._mat_hull(Pm, cell["P_accepted"]))
        eout = FULL._vec_hull(e, eacc)
        xout = _state_hull(xhat, xacc)
    else:
        Pout, eout, xout = cell["P_accepted"], eacc, xacc
    return Pout, eout, xout, cell


def _post_correction_q_upper(q: float, d: float) -> float:
    theta = POST._theta_from_q_upper(q)
    total = FULL.up(theta + d)
    if total >= math.pi:
        return math.inf
    return POST._q_from_theta_upper(total)


def _rotation_difference_upper(q: float) -> float:
    den = FULL.down(math.sqrt(FULL.down(4.0 + FULL.down(q*q))))
    if not den > 0.0:
        raise RuntimeError("finite rotation difference denominator lost positivity")
    return min(2.0, FULL.up(FULL.up(2.0*q)/den))


def _H_acc_bound(force_norm_upper: float):
    H = FULL._zero(3, N)
    b = FULL._box(force_norm_upper)
    H[0][1] = b; H[0][2] = b
    H[1][0] = b; H[1][2] = b
    H[2][0] = b; H[2][1] = b
    # Source orientation remains an actual rotation, but after the first
    # accepted/reset branch the full covariance is no longer axis-isotropic.
    # This diagnostic therefore uses the conservative entrywise orientation
    # enclosure and reports if that becomes the next dependency obstruction.
    for i in range(3):
        for j in range(3):
            H[i][15+j] = Interval(-1.0, 1.0)
    return H


def _norm_group(x, idxs) -> float:
    return FULL._norm_upper([x[i] for i in idxs])


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("sample-1 prefix domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("sample-1 prefix requires lever arm disabled")

    ENTRY.FULL3._install_backend()
    first = FIRST.build(domain_path, source_pieces=source_pieces)
    post = POST.build(domain_path, source_pieces=source_pieces)
    vector = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"post: {x}" for x in POST.validate(post)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = RG._source_phase_children(source_pieces)
    first_rows = first.get("source_cells", [])
    if len(src_phases) != len(first_rows):
        failures.append("sample-0 source cell ordering/count mismatch")

    h = float(FULL._source_cell()["dt_s"])
    gravity = float(domain["startup"]["gravity_mps2"])
    H0 = ENTRY._canonical_first_H(gravity)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    e0 = FULL._initial_error(domain)
    x0 = _zero_state()
    d0cap = float(first["max_first_accelerometer_correction_norm_upper_rad"])
    q0post = float(post["post_accel_cayley_norm_upper"])
    q1entry = RG._q_after_first_prediction(q0post, domain, h)
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])

    inverse_counts = {"FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": 0, "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE": 0}
    rows = []
    first_failure = None
    max_s_corr = 0.0
    max_acc_corr = 0.0
    max_q_after_s = q1entry
    max_residual = 0.0
    phase1_counts = {"due": 0, "not_due": 0}

    for si, (src, phase0) in enumerate(src_phases):
        try:
            F, Q, _ = FULL._transition_and_Q(src, domain)
            P0 = FULL._initial_covariance(src, domain_path)
            Pp = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))
            Pp = ENTRY._canonicalize_first_attitude_covariance(Pp, domain_path, h)
            ep = FULL._predict_error(e0, F)

            s0_backend = "NOT_DUE_IDENTITY"
            if phase0 == "due":
                Ppre, s0_backend = ENTRY._zero_residual_S_covariance(Pp, src)
                inverse_counts[s0_backend] += 1
            else:
                Ppre = Pp

            rho0 = float(first_rows[si]["combined_useful_residual_norm_upper_mps2"])
            Pa, ea, acell = ENTRY3._first_accel_covariance_and_state_tight(Ppre, ep, H0, Racc, rho0, d0cap)
            inverse_counts[acell["inverse_backend"]] += 1

            # Estimator mean after sample-0 accelerometer: accepted branch is
            # x+=dx; safe-LDLT identity branch is x unchanged.  Attitude is
            # immediately injected and cleared in both cases.
            xacc = _zero_state()
            for i in range(3, N):
                xacc[i] = acell["dx"][i]
            xpost0 = _state_hull(x0, xacc)

            P1 = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pa), matrix_transpose(F)), Q))
            e1 = FULL._predict_error(ea, F)
            x1 = _predict_estimator_state(xpost0, F)

            for phase1 in _phase1_options(src, phase0):
                phase1_counts[phase1] += 1
                P2, e2, x2 = P1, e1, x1
                q_after_s = q1entry
                s1_backend = "NOT_DUE_IDENTITY"
                ds = 0.0
                if phase1 == "due":
                    rS = [-x1[12+i] for i in range(3)]
                    P2, e2, x2, Scell = _apply_measurement_state(
                        P1, e1, x1, FULL._H_S(), FULL._R_S(src), rS,
                        allow_identity=True,
                    )
                    s1_backend = Scell["inverse_backend"]
                    inverse_counts[s1_backend] += 1
                    ds = FULL._norm_upper([-Scell["dx"][i] for i in range(3)])
                    max_s_corr = max(max_s_corr, ds)
                    q_after_s = max(q1entry, _post_correction_q_upper(q1entry, ds))
                    max_q_after_s = max(max_q_after_s, q_after_s)
                    if ds > DEPLOYED_CORRECTION_LIMIT_RAD:
                        raise RuntimeError(f"sample-1 S correction exceeds deployed 6-rad range: {ds}")
                    if not math.isfinite(q_after_s) or q_after_s >= 8.0:
                        raise RuntimeError(f"sample-1 S prefix leaves q<8 chart: q={q_after_s}")

                xaw = _norm_group(x2, FULL.AW)
                eaw = _norm_group(e2, FULL.AW)
                fpred = FULL.up(gravity + xaw)
                rot = _rotation_difference_upper(q_after_s)
                rho1 = FULL.up(FULL.up(rot*fpred) + FULL.up(eaw + ba))
                max_residual = max(max_residual, rho1)
                H1 = _H_acc_bound(fpred)
                r1 = FULL._vec_box(rho1)
                Acell1 = FULL._measurement_cell(P2, H1, Racc, r1)
                inverse_counts[Acell1["inverse_backend"]] += 1
                da = FULL._norm_upper([-Acell1["dx"][i] for i in range(3)])
                max_acc_corr = max(max_acc_corr, da)

                row = {
                    "source_phase_cell": si,
                    "sample0_pseudo_phase": phase0,
                    "sample1_pseudo_phase": phase1,
                    "tau_s": src["tau_s"].as_list(),
                    "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                    "R_S_filter_std": src["R_S_filter_std"].as_list(),
                    "sample0_S_inverse_backend": s0_backend,
                    "sample1_S_inverse_backend": s1_backend,
                    "sample1_S_correction_norm_upper_rad": ds,
                    "sample1_q_after_S_upper": q_after_s,
                    "sample1_estimator_aw_norm_upper_mps2": xaw,
                    "sample1_physical_aw_error_norm_upper_mps2": eaw,
                    "sample1_predicted_force_norm_upper_mps2": fpred,
                    "sample1_acc_residual_norm_upper_mps2": rho1,
                    "sample1_acc_inverse_backend": Acell1["inverse_backend"],
                    "sample1_acc_correction_norm_upper_rad": da,
                }
                rows.append(row)
                if Acell1["inverse_backend"] != "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN":
                    raise RuntimeError("sample-1 accelerometer innovation requires spectral fallback")
                if da > DEPLOYED_CORRECTION_LIMIT_RAD:
                    raise RuntimeError(f"sample-1 accelerometer correction exceeds deployed 6-rad range: {da}")
                q_after_a = _post_correction_q_upper(q_after_s, da)
                if not math.isfinite(q_after_a) or q_after_a >= 8.0:
                    raise RuntimeError(f"sample-1 accelerometer prefix leaves q<8 chart: q={q_after_a}")
        except Exception as exc:
            first_failure = {
                "source_phase_cell": si,
                "sample0_pseudo_phase": phase0,
                "reason": f"{type(exc).__name__}: {exc}",
                "rows_completed_before_failure": len(rows),
            }
            break

    closed = first_failure is None and bool(rows) and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SOURCE_CORRELATED_S_ACCEL_PREFIX_DIAGNOSTIC",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "estimator_mean_tracked_separately_from_physical_error": True,
        "sample1_S_residual_uses_estimator_mean_not_physical_error": True,
        "sample1_pseudo_phase_conditioned_on_sample0_timer_branch": True,
        "periodic_update_due_float_tolerance_retained": True,
        "sample1_accel_exact_rotation_difference_norm_used": True,
        "sample1_accel_generic_post_reset_orientation_box_used": True,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "sample1_entry_cayley_norm_upper": q1entry,
        "evaluated_sample1_paths": len(rows),
        "sample1_phase_counts": phase1_counts,
        "inverse_backend_counts": inverse_counts,
        "max_sample1_S_correction_norm_upper_rad": max_s_corr,
        "max_sample1_q_after_S_upper": max_q_after_s,
        "max_sample1_acc_residual_norm_upper_mps2": max_residual,
        "max_sample1_acc_correction_norm_upper_rad": max_acc_corr,
        "source_paths": rows,
        "first_failure": first_failure,
        "sample1_S_accel_prefix_closed": closed,
        "sample1_magnetometer_evaluated_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "PROPAGATE_SAMPLE1_ACCEPTED_RESET_CHILDREN_TO_MAGNETOMETER_AND_SAMPLE2"
            if closed else
            "REFINE_FIRST_FAILURE_WITH_POST_RESET_ROTATION_GAUGE_AND_DIRECTIONAL_SAMPLE1_ACCELERATOR_GAIN"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "estimator_mean_tracked_separately_from_physical_error",
        "sample1_S_residual_uses_estimator_mean_not_physical_error",
        "sample1_pseudo_phase_conditioned_on_sample0_timer_branch",
        "periodic_update_due_float_tolerance_retained",
        "sample1_accel_exact_rotation_difference_norm_used",
        "sample1_accel_generic_post_reset_orientation_box_used",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "sample1_magnetometer_evaluated_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction range changed")
    if int(d.get("evaluated_sample1_paths", 0)) < 0:
        failures.append("invalid sample-1 path count")
    status = d.get("P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE")
    if status == "PASS":
        if d.get("first_failure") is not None or d.get("sample1_S_accel_prefix_closed") is not True:
            failures.append("PASS retains a sample-1 failure")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_failure") is None:
            failures.append("sample-1 nonclosure lacks a concrete source witness")
    else:
        failures.append("invalid sample-1 prefix status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"],
        "paths": out["evaluated_sample1_paths"],
        "phase_counts": out["sample1_phase_counts"],
        "inverse_backends": out["inverse_backend_counts"],
        "max_S_d": out["max_sample1_S_correction_norm_upper_rad"],
        "max_q_after_S": out["max_sample1_q_after_S_upper"],
        "max_acc_residual": out["max_sample1_acc_residual_norm_upper_mps2"],
        "max_acc_d": out["max_sample1_acc_correction_norm_upper_rad"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
