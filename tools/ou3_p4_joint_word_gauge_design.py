#!/usr/bin/env python3
"""Co-gauged source design point for complete-word OU-III P4 dissipation.

The first joint-word trace showed that the first mandatory accelerometer
innovation was destroyed by independently hulling two objects that are
source-correlated: the anisotropic attitude covariance yaw axis and the gravity
/specific-force direction.  Its artificial S_xy interval was about +/-1.16,
large enough to defeat any verified inverse although every physical innovation
matrix is SPD.

This diagnostic keeps those objects in one orthogonal gauge.  It chooses the
admissible zero-body-rate source realization and the gravity-aligned handoff
frame d=e3.  In that frame the exact goLive attitude seed is

    diag(sigma_tilt^2, sigma_tilt^2, sigma_yaw^2),

and zero-rate prediction transports both the covariance yaw axis and gravity
without a relative rotation.  The zero-rate attitude/gyro-bias process block is
inserted from its exact integrated source formula.  Translation and active-bias
blocks still come from the existing validated source producer.

This is a DESIGN POINT ONLY.  It proves no source-complete P4 claim.  Its role
is to test whether the requested joint Joseph complete-word calculus has a
usable positive margin once the known orientation-hull artifact is removed.
The promoted proof must replace this point by a finite gauge/cone cover of all
admissible body motion, specific-force deviation, magnetic separation and
phased P2 tuner paths.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_full_process_ucc as PROCESS
import ou3_implementation_word_language as WORDS
import ou3_interval_ad as AD
import ou3_p4_candidate_full_word as CAND
import ou3_p4_joint_word_dissipation_design as D
import ou3_p4_joint_word_postprediction_design as POST
import ou3_p4_source_node_cells as NODES
import ou3_p5_full_h_prefix_cells as H
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = POST.DEFAULT_DOMAIN
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero_rate_transition_process(mode: str, src: dict, domain: dict):
    """Source-admissible omega=0 transition with exact correlated Q_ab block."""
    CAND._configure_mode(mode)
    F, Q, _Rstep, ba_meta = CAND._transition_and_Q(mode, src, domain)
    h = float(src["dt_s"])

    R0, B0 = H._rotation_step_box(0.0, h)
    for i in range(3):
        for j in range(3):
            F[i][j] = R0[i][j]
            F[i][3 + j] = B0[i][j]

    pc = PROCESS.build()["source_constants"]
    qg = float(pc["gyro_noise_density_rad_sqrt_s"]) ** 2
    qb = float(pc["gyro_bias_rw_variance_density"])
    qtt = qg * h + qb * h ** 3 / 3.0
    qtb = -qb * h ** 2 / 2.0
    qbb = qb * h
    for i in range(6):
        for j in range(6):
            Q[i][j] = I(0.0)
    for ax in range(3):
        Q[ax][ax] = Interval.outward_bounds(qtt, qtt)
        Q[ax][3 + ax] = Interval.outward_bounds(qtb, qtb)
        Q[3 + ax][ax] = Q[ax][3 + ax]
        Q[3 + ax][3 + ax] = Interval.outward_bounds(qbb, qbb)
    return F, H._psd_tighten(Q), ba_meta


def _gauged_golive_covariance(mode: str, src: dict, path: Path):
    """Exact H goLive covariance structure in d=e3 gauge; A adds source ba seed."""
    CAND._configure_mode(mode)
    seed = GOLIVE.build(path)["goLive_H_covariance_seed"]
    n = 18 if mode == "H" else 21
    P = D._zero(n, n)
    a = seed["attitude_covariance_seed"]
    vt = float(a["tilt_variance"])
    vy = float(a["gauged_yaw_variance"])
    for i, v in enumerate((vt, vt, vy)):
        P[i][i] = Interval.outward_bounds(v, v)

    # Bootstrap withholds MEKF time_update(), so gyro-bias covariance remains
    # the isotropic constructor seed at goLive.  Use the source value already
    # parsed by the prefix backend rather than the stale startup-RW upper box.
    pbg = float(H._source_pb0())
    for i in range(3, 6):
        P[i][i] = Interval.outward_bounds(pbg, pbg)

    for i in range(6, 9):
        v = float(seed["P_vv_variance_per_axis"])
        P[i][i] = Interval.outward_bounds(v, v)
    for i in range(9, 12):
        v = float(seed["P_pp_variance_per_axis"])
        P[i][i] = Interval.outward_bounds(v, v)
    for i in range(12, 15):
        v = float(seed["P_SS_variance_per_axis"])
        P[i][i] = Interval.outward_bounds(v, v)
    saw = src["sigma_aw_mps2"]
    aw = Interval.outward_bounds(saw.lo * saw.lo, saw.hi * saw.hi)
    for i in range(15, 18):
        P[i][i] = aw

    if mode == "A":
        s = CAND._source_sigma_bacc0()
        v = Interval.outward_bounds(s * s, s * s)
        for i in range(18, 21):
            P[i][i] = v
    return H._psd_tighten(P)


def _mode(mode: str, path: Path, domain: dict, source_node_index: int, cbox):
    CAND._configure_mode(mode)
    n = 18 if mode == "H" else 21
    src = NODES.h18_source_cell(source_node_index, NODES.build())
    F, Q, _ba_meta = _zero_rate_transition_process(mode, src, domain)
    Ppre = _gauged_golive_covariance(mode, src, path)
    Pm = H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Ppre), matrix_transpose(F)), Q))

    zpre = D._initial_ad(mode, domain, cbox)
    z = POST._rebase_postprediction(D._prediction(mode, zpre, F))
    Mupper, metric_meta = POST._entry_information_upper(mode)

    force, mag = D._canonical_vectors(domain)
    Ha, Hm, Hs = D._H_acc(mode, force, n), D._H_mag(mag, n), D._H_S(n)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    h = float(src["dt_s"])
    words = WORDS.build(path)
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = D._schedule(path, samples, h)
    signed_word = D._zero(n, n)
    ops = []

    for k in range(samples):
        if k in schedule["S_steps"]:
            r = [z[12 + i] for i in range(3)]
            eta = [AD.constant(0.0, n) for _ in range(3)]
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Hs, RS, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "S", **meta})
        if k in schedule["vector_steps"]:
            r = D._exact_acc_residual(mode, z, force)
            eta = D._eta(r, D._linear_residual(Ha, z))
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Ha, Racc, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "acc", **meta})

            r = D._exact_mag_residual(z, mag)
            eta = D._eta(r, D._linear_residual(Hm, z))
            Pm, z, dform, meta = D._ad_joint_update(Pm, z, Hm, Rmag, r, eta)
            signed_word = D._form_add(signed_word, dform)
            ops.append({"step": k, "kind": "mag", **meta})

        Pm = H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
        z = D._prediction(mode, z, F)

    mu = D._generalized_margin(signed_word, Mupper)
    return {
        "dimension": n,
        "source_node_index": source_node_index,
        "source_realization_body_rate_rad_s": 0.0,
        "gravity_yaw_covariance_gauge": "d=e3",
        "entry_cayley_box": cbox,
        "samples": samples,
        "schedule": schedule,
        "operations": ops,
        "operation_count": len(ops),
        "entry_information_metric_upper": metric_meta,
        "signed_word_generalized_margin_design": mu,
        "rho_homogeneous_design_upper": math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0,
        "K_interval_matrix_materialized": False,
    }


def build(path: Path = DEFAULT_DOMAIN, source_node_index: int = 0):
    path = Path(path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    q = 2.0 * math.tan(0.80 / 2.0)
    covers = CAND._ball_box_cover(q, max_box_norm_factor=1.5)
    cbox = covers[0]
    modes = {}
    failures = []
    for mode in ("H", "A"):
        try:
            modes[mode] = _mode(mode, path, domain, source_node_index, cbox)
        except Exception as exc:
            failures.append(f"{mode}: {type(exc).__name__}: {exc}")
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CO_GAUGED_JOINT_JOSEPH_WORD_DESIGN_POINT",
        "design_point_only": True,
        "source_complete": False,
        "source_realization_admissible": True,
        "trajectory_replay_used": False,
        "P1_changed": False,
        "P3_delta_used_as_physical_basin": False,
        "joint_P_H_R_r_used": True,
        "K_interval_matrix_materialized": False,
        "orientation_cartesian_covariance_hull_used": False,
        "outer_angle_rad": 0.80,
        "outer_cover_cells_total": len(covers),
        "outer_cover_cells_checked": 1 if modes else 0,
        "modes": modes,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.source_node_index)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "modes": {m: {
            "mu": d.get("modes", {}).get(m, {}).get("signed_word_generalized_margin_design"),
            "rho": d.get("modes", {}).get(m, {}).get("rho_homogeneous_design_upper"),
            "ops": d.get("modes", {}).get(m, {}).get("operation_count"),
            "last_op": (d.get("modes", {}).get(m, {}).get("operations") or [None])[-1],
        } for m in ("H", "A")},
        "failures": d["failures"],
    }, indent=2, sort_keys=True))
    return 0 if not d["failures"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
