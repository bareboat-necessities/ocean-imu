#!/usr/bin/env python3
"""Complete-word signed-Joseph design backend for OU-III P4.

This is the first executable backend following the final P4 proof calculus.  It
is intentionally a DESIGN diagnostic, not theorem promotion: one exact P2
source node, one outer-Cayley cover cell and one canonical PE vector realization
are evaluated before the same calculation is expanded to every phased P2 path
and source vector cell.

Unlike the abandoned interval-gain route, no interval matrix K is formed.
Every accepted update keeps the same P,H,R,r cell through

    x  = S^-1 r,
    dx = P H^T x,
    dW = r^T x - eta^T R^-1 eta.

The state map is interval-AD in the word-entry coordinates.  For each accepted
operation the derivative enclosures of r and eta are used to form the signed
directional quadratic family

    J_r^T S^-1 J_r - J_eta^T R^-1 J_eta,

and these forms are summed over the complete word before any scalarization.
The direct correction ``P H^T(S^-1 r)`` is used to propagate the nonlinear
state to the next operation.  Prediction dissipation from Q is nonnegative and
is omitted in this first lower bound; reset is an exact congruence.

The endpoint test compares the accumulated signed word form to the entry
information metric P0^-1.  A positive generalized margin mu would imply the
homogeneous design inequality W+ <= (1-mu) W.  Promotion remains blocked until
all outer-ball cells, all H/A source-vector cells, optional accepted branches,
and all phased P2 paths are outward validated.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval, matrix_add, matrix_mul, matrix_transpose,
    symmetric_positive_definite_ldlt,
)
import ou3_interval_ad as AD
import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_joint_joseph as JJ
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_p5_full_h_prefix_cells_v2 as H2
import ou3_verified_spd_inverse as VINV
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero(n: int, m: int):
    return [[I(0.0) for _ in range(m)] for _ in range(n)]


def _box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(math.nextafter(-a, -math.inf), math.nextafter(a, math.inf))


def _ad_matvec(A, x):
    n = x[0].n
    out = []
    for row in A:
        z = AD.constant(0.0, n)
        for a, b in zip(row, x):
            z = z + AD.constant(a, n) * b
        out.append(z)
    return out


def _mat_sub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _form(J, W):
    return matrix_mul(matrix_mul(matrix_transpose(J), W), J)


def _form_add(A, B):
    return matrix_add(A, B)


def _diag_inverse(R):
    n = len(R)
    X = _zero(n, n)
    for i in range(n):
        if not R[i][i].lo > 0.0:
            raise RuntimeError("R diagonal lost positivity")
        X[i][i] = Interval.outward_bounds(1.0 / R[i][i].hi, 1.0 / R[i][i].lo)
        for j in range(n):
            if i != j and not (R[i][j].lo == 0.0 and R[i][j].hi == 0.0):
                raise RuntimeError("design backend currently requires diagonal R")
    return X


def _signed_form(residual, eta, Sinv, R):
    Jr = AD.jacobian(residual)
    Je = AD.jacobian(eta)
    return _mat_sub(_form(Jr, Sinv), _form(Je, _diag_inverse(R)))


def _minus_mu(A, M, mu: float):
    n = len(A)
    q = Interval.outward_bounds(float(mu), float(mu))
    return [[A[i][j] - q * M[i][j] for j in range(n)] for i in range(n)]


def _generalized_margin(A, M):
    ok, _ = symmetric_positive_definite_ldlt(A)
    if not ok:
        return 0.0
    # Find a robust bracket from diagonal ratios, then binary search using only
    # validated interval LDLT.  The result is a strict lower generalized margin.
    hi = math.inf
    for i in range(len(A)):
        if M[i][i].hi > 0.0:
            hi = min(hi, A[i][i].lo / M[i][i].hi)
    if not math.isfinite(hi) or hi <= 0.0:
        return 0.0
    lo = 0.0
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        good, _ = symmetric_positive_definite_ldlt(_minus_mu(A, M, mid))
        if good:
            lo = mid
        else:
            hi = mid
    return math.nextafter(lo, -math.inf)


def _initial_ad(mode: str, domain: dict, cbox):
    CAND._configure_mode(mode)
    e, _ba, _pos = CAND._initial_error(mode, domain)
    vals = [Interval(a, b) for a, b in cbox] + list(e[3:])
    n = 18 if mode == "H" else 21
    if len(vals) != n:
        raise RuntimeError("initial AD state dimension mismatch")
    return [AD.independent(vals[i], i, n) for i in range(n)]


def _prediction(mode: str, z, F):
    n = len(z)
    out = list(z)
    # Homogeneous proof map: no deterministic gyro disturbance is inserted
    # here.  Its certified effect belongs in the later affine b term.
    c = list(z[:3])
    bg = list(z[3:6])
    Rstep = [[F[i][j] for j in range(3)] for i in range(3)]
    Bstep = [[F[i][3 + j] for j in range(3)] for i in range(3)]
    transported = _ad_matvec(Rstep, c)
    db = _ad_matvec(Bstep, bg)
    out[:3] = AD.deployed_correct_cayley(transported, db)
    for i in range(3, n):
        y = AD.constant(0.0, n)
        for j in range(3, n):
            y = y + AD.constant(F[i][j], n) * z[j]
        out[i] = y
    return out


def _canonical_vectors(domain: dict):
    live = domain["normal_live"]
    # Point realization at the declared information floors.  This is why the
    # file is design-only; final promotion must subdivide the full source vector
    # family, including magnitudes/orientation/covariance correlation.
    f = float(live["specific_force_norm_lower_mps2"])
    m = float(live["magnetic_vector_norm_lower_uT"])
    s = float(live["vector_sine_separation_lower"])
    c = math.sqrt(max(0.0, 1.0 - s * s))
    return [I(0.0), I(0.0), I(f)], [I(m * s), I(0.0), I(m * c)]


def _H_acc(mode: str, force, n: int):
    fx, fy, fz = force
    M = _zero(3, n)
    M[0][1] = fz; M[0][2] = -fy
    M[1][0] = -fz; M[1][2] = fx
    M[2][0] = fy; M[2][1] = -fx
    for i in range(3):
        M[i][15 + i] = I(1.0)
    if mode == "A":
        for i in range(3):
            M[i][18 + i] = I(1.0)
    return M


def _H_mag(mag, n: int):
    mx, my, mz = mag
    M = _zero(3, n)
    M[0][1] = mz; M[0][2] = -my
    M[1][0] = -mz; M[1][2] = mx
    M[2][0] = my; M[2][1] = -mx
    return M


def _H_S(n: int):
    M = _zero(3, n)
    for i in range(3):
        M[i][12 + i] = I(1.0)
    return M


def _exact_acc_residual(mode: str, z, force):
    n = len(z)
    R = AD.rotation_from_cayley(z[:3])
    f = [AD.constant(x, n) for x in force]
    aw = list(z[15:18])
    Rf = AD.matvec(R, f)
    Raw = AD.matvec(R, aw)
    r = [Rf[i] - f[i] + Raw[i] for i in range(3)]
    if mode == "A":
        r = [r[i] + z[18 + i] for i in range(3)]
    return r


def _exact_mag_residual(z, mag):
    n = len(z)
    R = AD.rotation_from_cayley(z[:3])
    m = [AD.constant(x, n) for x in mag]
    Rm = AD.matvec(R, m)
    return [Rm[i] - m[i] for i in range(3)]


def _linear_residual(Hm, z):
    return _ad_matvec(Hm, z)


def _eta(exact, linear):
    return [a - b for a, b in zip(exact, linear)]


def _ad_joint_update(Pm, z, Hm, Rm, residual, eta):
    PHt, S = JJ.innovation(Pm, Hm, Rm)
    Sinv, meta = JJ.verified_inverse(S)
    x = _ad_matvec(Sinv, residual)
    dx = _ad_matvec(PHt, x)
    signed = _signed_form(residual, eta, Sinv, Rm)

    # Covariance update uses the same verified S inverse but never forms K.
    Pj = JJ.posterior_covariance(Pm, PHt, Sinv, psd_tighten=H._psd_tighten)
    d = [-q for q in dx[:3]]
    out = list(z)
    out[:3] = AD.deployed_correct_cayley_right(z[:3], d)
    for i in range(3, len(z)):
        out[i] = z[i] - dx[i]
    Pout = H._reset_covariance(Pj, [q.val for q in dx[:3]])
    return Pout, out, signed, {
        "inverse_q_inf_upper": meta["neumann_q_inf_upper"],
        "correction_theta_norm_upper": AD._norm_upper([q.val for q in dx[:3]]),
        "K_interval_matrix_materialized": False,
    }


def _schedule(domain_path: Path, samples: int, h: float):
    words = WORDS.build(domain_path)
    wc = words["word_contract"]
    tr = wc["translation_recurrence"]
    q = int(tr["spread_index_q_W"])
    dS = float(tr["pseudo_gap_min_s"])
    sep = max(1, int(math.floor(q * dS / h)))
    last = samples - 1
    Ssteps = [last - 3 * sep, last - 2 * sep, last - sep, last]
    gap = VECTOR.build()["operating_envelope"]["packet_gap_s"]
    vg = max(1, int(math.ceil(float(gap[1]) / h)))
    Vsteps = [max(0, last - vg), last]
    if min(Ssteps) < 0:
        raise RuntimeError("mandatory four-S schedule does not fit word")
    return {"S_steps": Ssteps, "vector_steps": Vsteps}


def _mode(mode: str, path: Path, domain: dict, source_node_index: int, cbox):
    CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    nodes = NODES.build()
    src = NODES.h18_source_cell(source_node_index, nodes)
    F, Q, _Rstep, _ba = CAND._transition_and_Q(mode, src, domain)
    Pm = CAND._initial_covariance(mode, src, path)
    P0inv, p0meta = VINV.inverse_enclosure(Pm, symmetric_certified=True, spd_certified=True)
    z = _initial_ad(mode, domain, cbox)
    force, mag = _canonical_vectors(domain)
    Ha, Hm, Hs = _H_acc(mode, force, n), _H_mag(mag, n), _H_S(n)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    h = float(src["dt_s"])
    samples = int(WORDS.build(path)["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = _schedule(path, samples, h)
    D = _zero(n, n)
    ops = []

    for k in range(samples):
        Pm = H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
        z = _prediction(mode, z, F)

        if k in schedule["S_steps"]:
            r = [z[12 + i] for i in range(3)]
            eta = [AD.constant(0.0, n) for _ in range(3)]
            Pm, z, dform, meta = _ad_joint_update(Pm, z, Hs, RS, r, eta)
            D = _form_add(D, dform)
            ops.append({"step": k, "kind": "S", **meta})

        if k in schedule["vector_steps"]:
            r = _exact_acc_residual(mode, z, force)
            eta = _eta(r, _linear_residual(Ha, z))
            Pm, z, dform, meta = _ad_joint_update(Pm, z, Ha, Racc, r, eta)
            D = _form_add(D, dform)
            ops.append({"step": k, "kind": "acc", **meta})

            r = _exact_mag_residual(z, mag)
            eta = _eta(r, _linear_residual(Hm, z))
            Pm, z, dform, meta = _ad_joint_update(Pm, z, Hm, Rmag, r, eta)
            D = _form_add(D, dform)
            ops.append({"step": k, "kind": "mag", **meta})

    mu = _generalized_margin(D, P0inv)
    return {
        "dimension": n,
        "source_node_index": source_node_index,
        "entry_cayley_box": cbox,
        "samples": samples,
        "schedule": schedule,
        "operations": ops,
        "operation_count": len(ops),
        "entry_metric_inverse_q_inf_upper": p0meta["neumann_q_inf_upper"],
        "signed_word_generalized_margin_design": mu,
        "rho_homogeneous_design_upper": math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0,
        "signed_word_form_interval_ldlt_positive": _generalized_margin(D, P0inv) > 0.0,
        "K_interval_matrix_materialized": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          ball_inflation: float = 1.5) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    words = WORDS.build(path)
    failures = [f"word: {x}" for x in WORDS.validate(words)]
    q = 2.0 * math.tan(0.80 / 2.0)
    covers = CAND._ball_box_cover(q, max_box_norm_factor=ball_inflation)
    cbox = covers[0]
    modes = {}
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode(mode, path, domain, int(source_node_index), cbox)
        except Exception as exc:
            failures.append(f"{mode}: {type(exc).__name__}: {exc}")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_JOINT_JOSEPH_COMPLETE_WORD_DISSIPATION_DESIGN",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "design_only_not_theorem_promotion": True,
        "P1_changed": False,
        "P3_delta_used_as_physical_basin": False,
        "joint_P_H_R_r_used_through_innovation_solve": True,
        "K_interval_matrix_materialized": False,
        "signed_Joseph_identity_used": "r^T S^-1 r - eta^T R^-1 eta",
        "directional_forms_accumulated_before_scalarization": True,
        "prediction_process_dissipation_credited": False,
        "homogeneous_disturbance_set_to_zero_for_rho_design": True,
        "affine_b_term_established_here": False,
        "outer_angle_rad": 0.80,
        "outer_ball_cover_total": len(covers),
        "outer_ball_cells_checked": 1 if modes else 0,
        "source_node_index": int(source_node_index),
        "phased_P2_paths_checked": False,
        "optional_accepted_branch_family_checked": False,
        "full_source_vector_family_checked": False,
        "modes": modes,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d: dict):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "design_only_not_theorem_promotion",
        "joint_P_H_R_r_used_through_innovation_solve",
        "directional_forms_accumulated_before_scalarization",
        "homogeneous_disturbance_set_to_zero_for_rho_design",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "P1_changed",
        "P3_delta_used_as_physical_basin", "K_interval_matrix_materialized",
        "prediction_process_dissipation_credited", "affine_b_term_established_here",
        "phased_P2_paths_checked", "optional_accepted_branch_family_checked",
        "full_source_vector_family_checked", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED", "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_node_index=a.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu_design": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho_design": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "ops": d.get("modes", {}).get(m, {}).get("operation_count"),
        } for m in ("H", "A")},
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
