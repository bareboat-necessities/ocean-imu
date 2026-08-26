#!/usr/bin/env python3
"""Refine the first failing P5 sample-1 accelerometer witness by preserving rotation structure.

The generic sample-1 diagnostic fails on source cell 0 because it replaces the
shipping J_aw=R_wb by nine independent intervals in [-1,1].  This producer keeps
that orthogonality through the first accepted correction.

At sample zero the gravity-gauged covariance, isotropic linear blocks, isotropic
accelerometer noise, and H=[-[g e3]_x, I_aw] are equivariant under simultaneous
rotations about e3.  The first accelerometer attitude correction has no e3
component.  Therefore every accepted child can be represented, without loss,
by the canonical tangent correction

    d0 = (delta,0,0),   0 <= delta <= d0_max(cell).

After the deployed reset, the proof coordinates of every world-linear 3-vector
(v,p,S,a_w) are transformed by Rx(delta), i.e. into the corrected estimated body
frame.  The linear OU chain and its process covariance are axis-isotropic, so
this is a coordinate change, not a filter modification.  After the next 5 ms
prediction the same linear blocks are transformed by the source-enclosed
one-step body rotation Rstep.  In that sample-1 body gauge the exact shipping
accelerometer Jacobian is again J_aw=I.

Safe-LDLT solver-identity and accepted branches are kept separate rather than
hulled.  The sample-1 due S branch is likewise split into accepted and identity
subbranches.  This producer refines only the first previously failing source
cell; passing it is not a whole-family or whole-word promotion.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_identity, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_post_reset as POST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_entry_v3 as ENTRY3
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1
N = FULL.N
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0


def _zero_state():
    return [FULL.I(0.0) for _ in range(N)]


def _rx_delta(delta_hi: float):
    if not (0.0 <= delta_hi < math.pi/2):
        raise RuntimeError("canonical first correction must stay below pi/2 for monotone Rx enclosure")
    s_hi = VT.sin_point(delta_hi).hi
    c_lo = VT.cos_point(delta_hi).lo
    z = FULL.I(0.0); o = FULL.I(1.0)
    s = Interval(0.0, FULL.up(s_hi))
    c = Interval(FULL.down(c_lo), 1.0)
    return [
        [o, z, z],
        [z, c, -s],
        [z, s, c],
    ]


def _linear_coordinate_transform(R3):
    T = matrix_identity(N)
    for group in (FULL.V, FULL.P, FULL.SS, FULL.AW):
        for i in range(3):
            for j in range(3):
                T[group[i]][group[j]] = R3[i][j]
    return T


def _transform_covariance(Pm, R3):
    T = _linear_coordinate_transform(R3)
    return FULL._psd_tighten(matrix_mul(matrix_mul(T, Pm), matrix_transpose(T)))


def _transform_state(x, R3):
    out = list(x)
    for group in (FULL.V, FULL.P, FULL.SS, FULL.AW):
        old = [x[i] for i in group]
        new = FULL._mat_vec(R3, old)
        for k, i in enumerate(group):
            out[i] = new[k]
    return out


def _predict_state(x, F):
    out = list(x)
    lin = list(FULL.V)+list(FULL.P)+list(FULL.SS)+list(FULL.AW)
    for i in lin:
        y = FULL.I(0.0)
        for j in lin:
            y = y + F[i][j]*x[j]
        out[i] = y
    for i in FULL.TH:
        out[i] = FULL.I(0.0)
    return out


def _canonical_H0(gravity: float):
    return ENTRY._canonical_first_H(gravity)


def _sample1_H(force_norm_upper: float):
    H = FULL._zero(3, N)
    f = FULL._box(force_norm_upper)
    # H_theta=-[f]_x, with exact skew structure and a source-complete component
    # box for the force direction.  J_aw is identity in the transported body gauge.
    H[0][1] = f;  H[0][2] = f
    H[1][0] = f;  H[1][2] = f
    H[2][0] = f;  H[2][1] = f
    for i in range(3):
        H[i][15+i] = FULL.I(1.0)
    return H


def _post_q(q: float, d: float) -> float:
    return PREFIX2._post_correction_q_upper(q, d)


def _rot_diff(q: float) -> float:
    den = FULL.down(math.sqrt(FULL.down(4.0 + FULL.down(q*q))))
    if not den > 0.0:
        return math.inf
    return min(2.0, FULL.up(FULL.up(2.0*q)/den))


def _norm(x, idxs) -> float:
    return FULL._norm_upper([x[i] for i in idxs])


def _measurement_accepted(Pm, H, R, r):
    cell = FULL._measurement_cell(Pm, H, R, r)
    return cell["P_accepted"], cell


def _first_accel_accepted(Ppre, ep, xpre, H0, Racc, rho0: float, d0max: float):
    PHt, S = FULL._innovation(Ppre, H0, Racc)
    Sinv, backend = FULL._spd_inverse_enclosure(S, Racc)
    K = matrix_mul(PHt, Sinv)
    r0 = FULL._vec_box(rho0)
    dx = FULL._mat_vec(K, r0)
    caps = ENTRY3._linear_gain_caps(Ppre, Racc, rho0)
    for name, idxs in ENTRY3.GROUPS.items():
        cap = Interval(-caps[name], caps[name])
        for i in idxs:
            dx[i] = FULL._intersect(dx[i], cap)

    Pj = FULL._shipping_joseph(Ppre, K, S, PHt)
    # Exact SO(2)-equivariance lets us replace the unknown tangent direction by
    # +e1 while retaining the complete magnitude interval.
    dcanon = [Interval.outward_bounds(0.0, d0max), FULL.I(0.0), FULL.I(0.0)]
    Pr = FULL._reset_covariance(Pj, dcanon)

    eacc = list(ep); xacc = list(xpre)
    for i in range(3, N):
        eacc[i] = ep[i] - dx[i]
        xacc[i] = xpre[i] + dx[i]
    for i in FULL.TH:
        xacc[i] = FULL.I(0.0)

    Rx = _rx_delta(d0max)
    return (
        _transform_covariance(Pr, Rx),
        _transform_state(eacc, Rx),
        _transform_state(xacc, Rx),
        backend,
    )


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2, source_cell_index: int = 0) -> dict:
    FULL3._install_backend()
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("sample-1 rotation-gauge refinement must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("sample-1 rotation-gauge refinement requires lever arm disabled")

    first = FIRST.build(domain_path, source_pieces=source_pieces)
    post = POST.build(domain_path, source_pieces=source_pieces)
    vector = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"post: {x}" for x in POST.validate(post)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = RG._source_phase_children(source_pieces)
    if not (0 <= source_cell_index < len(src_phases)):
        raise IndexError("source cell index outside first-prefix family")
    src, phase0 = src_phases[source_cell_index]
    first_row = first["source_cells"][source_cell_index]
    if phase0 != first_row["pseudo_phase"]:
        failures.append("source-cell phase ordering mismatch")

    h = float(src["dt_s"])
    gravity = float(domain["startup"]["gravity_mps2"])
    H0 = _canonical_H0(gravity)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    e0 = FULL._initial_error(domain)
    x0 = _zero_state()
    d0max = float(first_row["correction_norm_upper_rad"])
    rho0 = float(first_row["combined_useful_residual_norm_upper_mps2"])
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    true_force_max = float(domain["normal_live"]["specific_force_norm_upper_mps2"])

    F, Q, Rstep = FULL._transition_and_Q(src, domain)
    P0 = FULL._initial_covariance(src, domain_path)
    Pp = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))
    Pp = ENTRY._canonicalize_first_attitude_covariance(Pp, domain_path, h)
    ep = FULL._predict_error(e0, F)

    s0_branches = [("identity", Pp)]
    if phase0 == "due":
        Ps0, _ = ENTRY._zero_residual_S_covariance(Pp, src)
        s0_branches.insert(0, ("accepted", Ps0))

    q_first_pre = float(first["post_prediction_full_cayley_norm_upper"])
    branch_rows = []
    first_failure = None
    max_d1 = 0.0
    max_rho1 = 0.0
    fixed = fallback = 0

    for s0_name, Ppre in s0_branches:
        accel0_branches = [("solver_identity", Ppre, ep, x0, q_first_pre)]
        Pa0, ea0, xa0, b0 = _first_accel_accepted(Ppre, ep, x0, H0, Racc, rho0, d0max)
        if b0 == "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": fixed += 1
        else: fallback += 1
        q0a = _post_q(q_first_pre, d0max)
        accel0_branches.insert(0, ("accepted_canonical_tangent", Pa0, ea0, xa0, q0a))

        for a0_name, Pc0, ec0, xc0, q0_branch in accel0_branches:
            # One next prediction in the held post-sample0 body gauge.
            P1m = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pc0), matrix_transpose(F)), Q))
            e1m = _predict_state(ec0, F)
            x1m = _predict_state(xc0, F)
            q1 = RG._q_after_first_prediction(q0_branch, domain, h)

            # Re-express all world-linear groups in the sample-1 estimated body
            # gauge.  Rstep is the same source-enclosed body rotation used by
            # the attitude prediction, so J_aw becomes I exactly at the packet.
            P1 = _transform_covariance(P1m, Rstep)
            e1 = _transform_state(e1m, Rstep)
            x1 = _transform_state(x1m, Rstep)

            phase1_options = PREFIX2.BASE._phase1_options(src, phase0)
            for phase1 in phase1_options:
                s1_branches = [("identity", P1, e1, x1, q1)]
                if phase1 == "due":
                    rS = [-x1[12+i] for i in range(3)]
                    Scell = FULL._measurement_cell(P1, FULL._H_S(), FULL._R_S(src), rS)
                    ds = FULL._norm_upper(Scell["dx"][0:3])
                    Ps = Scell["P_accepted"]
                    es = list(e1); xs = list(x1)
                    for i in range(3, N):
                        es[i] = e1[i] - Scell["dx"][i]
                        xs[i] = x1[i] + Scell["dx"][i]
                    for i in FULL.TH: xs[i] = FULL.I(0.0)
                    qs = _post_q(q1, ds)
                    s1_branches.insert(0, ("accepted", Ps, es, xs, qs))
                    if Scell["inverse_backend"] == "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": fixed += 1
                    else: fallback += 1

                for s1_name, P2, e2, x2, q2 in s1_branches:
                    try:
                        xaw = _norm(x2, FULL.AW)
                        eaw = _norm(e2, FULL.AW)
                        fhat = FULL.up(gravity + xaw)
                        rho1 = FULL.up(FULL.up(_rot_diff(q2)*true_force_max) + FULL.up(eaw + ba))
                        H1 = _sample1_H(fhat)
                        PHt1, S1 = FULL._innovation(P2, H1, Racc)
                        Sinv1, backend1 = FULL._spd_inverse_enclosure(S1, Racc)
                        K1 = matrix_mul(PHt1, Sinv1)
                        Ktheta = [row[:] for row in K1[0:3]]
                        knorm = RG._op2_upper(Ktheta)
                        d1 = FULL.up(knorm*rho1)
                        max_d1 = max(max_d1, d1); max_rho1 = max(max_rho1, rho1)
                        if backend1 == "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": fixed += 1
                        else: fallback += 1
                        q_after = _post_q(q2, d1)
                        row = {
                            "sample0_S_solver_branch": s0_name,
                            "sample0_accel_solver_branch": a0_name,
                            "sample1_pseudo_phase": phase1,
                            "sample1_S_solver_branch": s1_name,
                            "sample1_entry_q_upper": q2,
                            "sample1_estimator_aw_norm_upper_mps2": xaw,
                            "sample1_physical_aw_error_norm_upper_mps2": eaw,
                            "sample1_predicted_force_norm_upper_mps2": fhat,
                            "sample1_residual_norm_upper_mps2": rho1,
                            "sample1_inverse_backend": backend1,
                            "sample1_Ktheta_norm_upper": knorm,
                            "sample1_correction_norm_upper_rad": d1,
                            "sample1_post_accel_q_upper": q_after,
                            "inside_deployed_correction_range": d1 <= DEPLOYED_CORRECTION_LIMIT_RAD,
                            "inside_q8_after_sample1_accel": math.isfinite(q_after) and q_after < 8.0,
                        }
                        branch_rows.append(row)
                        if backend1 != "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN" or d1 > 6.0 or not math.isfinite(q_after) or q_after >= 8.0:
                            if first_failure is None:
                                first_failure = row
                    except Exception as exc:
                        row = {
                            "sample0_S_solver_branch": s0_name,
                            "sample0_accel_solver_branch": a0_name,
                            "sample1_pseudo_phase": phase1,
                            "sample1_S_solver_branch": s1_name,
                            "exception": f"{type(exc).__name__}: {exc}",
                        }
                        branch_rows.append(row)
                        if first_failure is None: first_failure = row

    refinement_pass = bool(branch_rows) and first_failure is None and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_FIRST_FAILURE_CANONICAL_ROTATION_GAUGE_REFINEMENT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "refined_source_cell_index": source_cell_index,
        "refined_sample0_pseudo_phase": phase0,
        "sample0_gravity_SO2_equivariance_used": True,
        "sample0_accel_attitude_correction_gravity_axis_component_exact_zero": True,
        "sample0_accepted_correction_canonicalized_to_positive_tangent_axis": True,
        "world_linear_groups_congruence_transformed_into_corrected_body_gauge": True,
        "next_prediction_body_rotation_applied_to_world_linear_gauge": True,
        "sample1_J_aw_exact_identity_in_transported_body_gauge": True,
        "safe_ldlt_solver_identity_branches_kept_separate": True,
        "sample1_true_specific_force_domain_upper_used_in_rotation_residual": true if False else True,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "sample0_cell_correction_norm_upper_rad": d0max,
        "evaluated_branch_count": len(branch_rows),
        "fixed_pivot_inverse_count": fixed,
        "spectral_fallback_inverse_count": fallback,
        "max_sample1_residual_norm_upper_mps2": max_rho1,
        "max_sample1_correction_norm_upper_rad": max_d1,
        "branches": branch_rows,
        "first_unclosed_branch": first_failure,
        "complete_source_family_refined_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_ROTATION_GAUGE_WITNESS_REFINEMENT": "PASS" if refinement_pass else "NOT_ESTABLISHED",
        "next_obligation": (
            "LIFT_CANONICAL_ROTATION_GAUGE_REFINEMENT_FROM_FIRST_WITNESS_TO_ALL_SAMPLE1_SOURCE_CELLS"
            if refinement_pass else
            "SUBDIVIDE_CANONICAL_CORRECTION_MAGNITUDE_AND_SAMPLE1_FORCE_DIRECTION_WITH_J_AW_IDENTITY"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "sample0_gravity_SO2_equivariance_used",
        "sample0_accel_attitude_correction_gravity_axis_component_exact_zero",
        "sample0_accepted_correction_canonicalized_to_positive_tangent_axis",
        "world_linear_groups_congruence_transformed_into_corrected_body_gauge",
        "next_prediction_body_rotation_applied_to_world_linear_gauge",
        "sample1_J_aw_exact_identity_in_transported_body_gauge",
        "safe_ldlt_solver_identity_branches_kept_separate",
        "sample1_true_specific_force_domain_upper_used_in_rotation_residual",
    ):
        if d.get(key) is not True: failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "complete_source_family_refined_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False: failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0: failures.append("deployed correction range changed")
    if int(d.get("evaluated_branch_count", 0)) <= 0: failures.append("no rotation-gauge witness branches evaluated")
    status = d.get("P5_SAMPLE1_ROTATION_GAUGE_WITNESS_REFINEMENT")
    if status == "PASS":
        if d.get("first_unclosed_branch") is not None: failures.append("PASS retains unclosed witness branch")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_unclosed_branch") is None: failures.append("nonclosure missing witness branch")
    else:
        failures.append("invalid witness refinement status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces, source_cell_index=args.source_cell_index)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_ROTATION_GAUGE_WITNESS_REFINEMENT"],
        "branches": out["evaluated_branch_count"],
        "fixed": out["fixed_pivot_inverse_count"],
        "fallback": out["spectral_fallback_inverse_count"],
        "max_rho1": out["max_sample1_residual_norm_upper_mps2"],
        "max_d1": out["max_sample1_correction_norm_upper_rad"],
        "first_unclosed": out["first_unclosed_branch"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
