#!/usr/bin/env python3
"""Outward P5 source-cell propagation for the later H-word prefixes.

This is the first numerical layer after the exact effective-vector-input
lemmas.  It deliberately carries one *joint* source tuple through every scalar
bound used for a correction cell instead of multiplying independently selected
global extrema.  Each tuple contains

  (tau, sigma_aw, R_S, q-cell, P envelope, H, R, S, K, r, d_eff).

The source parameter partition is the same validated partition used by the P3
matrix certificate.  The attitude Cayley partition is the outward annular
partition already certified by P5.  For every product cell this module computes
source-correlated innovation/gain envelopes for the shipping S, accelerometer,
and magnetometer corrections.  The configured magnetometer uses the exact
radial-null / tangent-only identity and the accelerometer uses the exact
``a_w`` effective input; no standalone vector-eta penalty is reintroduced.

This stage also makes an important fail-closed distinction.  The P3 source-cell
producer currently exports a directional Loewner envelope, not the full signed
matrix-valued covariance cell.  Such an envelope is sufficient for outward
S/K norm bounds, but it is *not* sufficient to recover the signed ``a^T c``
needed by the exact Cayley composition.  Therefore this producer propagates the
complete scalar P/H/R/S/K/r/d_eff enclosure and records the first place where
full matrix-direction correlation is still required.  It never replaces that
missing correlation by ``-|a||c|`` and never promotes P5 from a norm-only cell.

The purpose is twofold: (1) turn the previous prose obligation into numerical
cells with explicit worst tuples and margins, and (2) identify whether the next
subdivision must target covariance magnitude, effective-vector geometry, or the
signed K*r direction.  A later full matrix interval cell producer can replace
the directional P envelope without changing this interface.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_s_exact_prefix as FIRSTS
import ou3_p5_first_s_state_prefix_certificate as FIRSTSTATE
import ou3_p5_mag_information_reduction as MAGINFO
import ou3_source_reachable_matrix_p3 as P3CELL
import ou3_vector_uco_certificate as VECTOR
from ou3_interval import Interval

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _sqrt_up(x: float) -> float:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("finite nonnegative square-root input required")
    return up(math.sqrt(x))


def _source_cells(domain: dict) -> tuple[list[dict], dict]:
    """Regenerate the P3 H source partition and retain every joint tuple.

    ``mode_cell`` is intentionally called with the same x/sigma/R_S cell.  No
    extremum from another cell is substituted into its covariance envelope.
    """
    live = domain["normal_live"]
    vector = VECTOR.build()
    process = PROCESS.build()
    vf = VECTOR.validate(vector)
    pf = PROCESS.validate(process)
    if vf or pf:
        raise RuntimeError(f"source-cell prerequisites failed: vector={vf}, process={pf}")

    sched = P3CELL.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    xlo, xhi = h / tau_hi, h / tau_lo
    edges = P3CELL.geom_edges(xlo, xhi, 24)
    if xlo < P3CELL.BRANCH_X < xhi:
        edges = sorted(set(edges + [P3CELL.BRANCH_X]))
    xcells: list[tuple[Interval, float]] = []
    for x in P3CELL.interval_cells(edges):
        xcells.extend(P3CELL.split_x_cell(x))
    sigmas = P3CELL.interval_cells(P3CELL.geom_edges(0.05, 6.0, 5))
    rs_lo, rs_hi = map(float, sched["R_S_applied_invariant"])
    rss = P3CELL.interval_cells(P3CELL.geom_edges(rs_lo, rs_hi, 8))
    alpha6 = P3CELL.pos(P3CELL.vector_alpha6(live, vector), "declared alpha6")

    out: list[dict] = []
    index = 0
    for x, rho_t in xcells:
        for sigma in sigmas:
            for rs in rss:
                c = P3CELL.mode_cell(
                    "H", x, rho_t, sigma, rs, live, vector, process, sched, alpha6
                )
                out.append({
                    "index": index,
                    "x_h_over_tau": list(c["x_h_over_tau"]),
                    "tau_s": list(c["tau_s"]),
                    "sigma_aw_mps2": list(c["sigma_aw_mps2"]),
                    "R_S_filter_std": list(c["R_S_filter_std"]),
                    "P_directional_diagonal_upper": list(c["Sigma_diagonal_upper"]),
                    "P_scaled_lambda_max_upper": float(c["Sigma_scaled_lambda_max_upper"]),
                    "P_relative_Riccati_injection_margin_lower": float(c["relative_Riccati_injection_margin_lower"]),
                    "word_horizon_s_upper": float(c["word_horizon_s_upper"]),
                })
                index += 1
    return out, {
        "x_cells": len(xcells),
        "sigma_cells": len(sigmas),
        "R_S_cells": len(rss),
        "joint_source_cells": len(out),
        "source_schedule": sched,
        "vector": vector,
    }


def _mag_cell(src: dict, qrow: dict, live: dict, vc: dict) -> dict:
    p = src["P_directional_diagonal_upper"]
    ptheta = max(map(float, p[0:3]))
    mlo = float(live["magnetic_vector_norm_lower_uT"])
    mhi = float(live["magnetic_vector_norm_upper_uT"])
    rm = up(float(vc["mag_measurement_std_uT"]) ** 2)
    qlo, qhi = map(float, qrow["q_interval"])

    s_lo = down(rm)
    s_hi = up(rm + up(mhi * mhi * ptheta))
    ktheta_hi = _sqrt_up(ptheta / rm)
    geff_lo = float(qrow["effective_tangent_gain_lower"])
    geff_hi = float(qrow["effective_tangent_gain_upper"])
    deff_hi = up(geff_hi * qhi)
    r_tangent_hi = up(mhi * deff_hi)
    dtheta_hi = up(ktheta_hi * r_tangent_hi)

    # Same source tuple and same q cell: the useful term is lower-bounded with
    # its own innovation upper and the exact tangent-only d_eff coercivity.
    useful_per_cperp2 = down((mlo * mlo * geff_lo * geff_lo) / s_hi)
    tangent_penalty_ratio = float(
        qrow["effective_vs_linear_tangent_penalty_information_ratio_upper"]
    )
    return {
        "P_theta_lambda_upper": ptheta,
        "H_theta_operator_norm_interval_uT": [down(mlo), up(mhi)],
        "R_variance": rm,
        "S_innovation_lambda_interval": [s_lo, s_hi],
        "K_theta_operator_norm_upper": ktheta_hi,
        "d_eff_norm_upper": deff_hi,
        "r_tangent_norm_upper": r_tangent_hi,
        "attitude_injection_norm_upper": dtheta_hi,
        "effective_tangent_information_per_cperp2_lower": useful_per_cperp2,
        "linear_cayley_tangent_penalty_ratio_upper": tangent_penalty_ratio,
        "radial_K_action_exact_zero": True,
        "radial_Joseph_information_exact_zero": True,
        "standalone_eta_penalty_used": False,
        "q_interval": [qlo, qhi],
    }


def _acc_cell(src: dict, qrow: dict, live: dict, vc: dict, aw_radius: float) -> dict:
    p = src["P_directional_diagonal_upper"]
    ptheta = max(map(float, p[0:3]))
    paw = max(map(float, p[15:18]))
    flo = float(live["specific_force_norm_lower_mps2"])
    fhi = float(live["specific_force_norm_upper_mps2"])
    ra = up(float(vc["acc_measurement_std_mps2"]) ** 2)
    qlo, qhi = map(float, qrow["q_interval"])

    # For the same source cell, Cauchy on the theta/a_w covariance gives
    # H P H^T <= (|f| sqrt(Ptheta)+sqrt(Paw))^2.  This keeps the two blocks in
    # one cell rather than selecting Ptheta and Paw from unrelated cells.
    root = up(fhi * _sqrt_up(ptheta) + _sqrt_up(paw))
    s_hi = up(ra + up(root * root))
    ktheta_hi = _sqrt_up(ptheta / ra)

    a_att = float(qrow["acc_effective_aw_attitude_eta_per_vector_norm_upper"])
    a_lat = float(qrow["acc_effective_aw_latent_cross_gain_upper"])
    eeta_hi = up(up(a_att * fhi) + up(a_lat * aw_radius))
    return {
        "P_theta_lambda_upper": ptheta,
        "P_aw_lambda_upper": paw,
        "H_att_operator_norm_interval_mps2": [down(flo), up(fhi)],
        "H_aw_operator_norm": 1.0,
        "R_variance": ra,
        "S_innovation_lambda_interval": [down(ra), s_hi],
        "K_theta_operator_norm_upper": ktheta_hi,
        "effective_aw_input_norm_upper_mps2": eeta_hi,
        "effective_aw_attitude_coefficient": a_att,
        "effective_aw_latent_coefficient": a_lat,
        "standalone_eta_penalty_used": False,
        "exact_effective_input_identity_used": True,
        "q_interval": [qlo, qhi],
    }


def _S_cell(src: dict, S_radius: float) -> dict:
    p = src["P_directional_diagonal_upper"]
    ptheta = max(map(float, p[0:3]))
    pS = max(map(float, p[12:15]))
    rs_lo, rs_hi = map(float, src["R_S_filter_std"])
    if not rs_lo > 0.0:
        raise RuntimeError("positive R_S standard deviation required")
    rlo2 = down(rs_lo * rs_lo)
    rhi2 = up(rs_hi * rs_hi)
    # Functional-covariance extremum: ||C_thetaS(PSS+R)^-1|| <=
    # sqrt(Ptheta)/(2 sqrt(R)) for every PSD joint covariance.  It keeps the
    # complete S->attitude gain and is finite on each R_S source cell.
    kthetaS = up(_sqrt_up(ptheta) / (2.0 * rs_lo))
    return {
        "P_theta_lambda_upper": ptheta,
        "P_S_lambda_upper": pS,
        "H_S_operator_norm": 1.0,
        "R_std_interval": [rs_lo, rs_hi],
        "R_variance_interval": [rlo2, rhi2],
        "S_innovation_lambda_interval": [rlo2, up(pS + rhi2)],
        "K_thetaS_operator_norm_upper": kthetaS,
        "r_S_state_norm_upper": float(S_radius),
        "attitude_injection_norm_upper": up(kthetaS * float(S_radius)),
        "eta_exact_zero": True,
        "full_S_to_attitude_gain_retained": True,
    }


def _worse(cur: dict | None, value: float, row: dict) -> dict:
    if cur is None or value > float(cur["value"]):
        return {"value": float(value), **row}
    return cur


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("joint-prefix source domain must not be trajectory fitted")

    first = FIRSTS.build(domain_path)
    firststate = FIRSTSTATE.build(domain_path)
    veff = VEFF.build(domain_path)
    maginfo = MAGINFO.build(domain_path)
    failures = [f"first-S: {x}" for x in FIRSTS.validate(first)]
    failures += [f"first-S-state: {x}" for x in FIRSTSTATE.validate(firststate)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"mag-information: {x}" for x in MAGINFO.validate(maginfo)]

    source_cells, meta = _source_cells(domain)
    qcells = list(maginfo["annular_information_cells"])
    veff_by_index = {int(x["index"]): x for x in veff["annular_effective_input_cells"]}
    if len(qcells) != len(veff_by_index):
        failures.append("effective and magnetic q partitions differ")

    live = domain["normal_live"]
    vc = meta["vector"]["configured_measurement_bounds"]
    aw_radius = float(domain["startup"]["physical_handoff_coordinate_bounds"]["latent_acceleration_error_norm_upper_mps2"])
    S_radius = float(firststate["first_due_S_error_norm_upper_m_s"])

    worst_mag_d = None
    worst_acc_e = None
    worst_S_d = None
    min_mag_info = None
    product_cells = 0
    finite = True

    # Do not serialize tens of thousands of cells.  The emitted witnesses retain
    # the complete source tuple and q interval that attain each extremum.
    for src in source_cells:
        Scell = _S_cell(src, S_radius)
        worst_S_d = _worse(
            worst_S_d, Scell["attitude_injection_norm_upper"],
            {"source_cell": {k: src[k] for k in ("index", "tau_s", "sigma_aw_mps2", "R_S_filter_std")}, "S_cell": Scell},
        )
        for qrow in qcells:
            idx = int(qrow["index"])
            vr = veff_by_index[idx]
            # Merge the two exact q-cell payloads without taking data from a
            # different annulus.
            qjoint = dict(qrow)
            qjoint.update({
                "acc_effective_aw_attitude_eta_per_vector_norm_upper": vr["acc_effective_aw_attitude_eta_per_vector_norm_upper"],
                "acc_effective_aw_latent_cross_gain_upper": vr["acc_effective_aw_latent_cross_gain_upper"],
            })
            m = _mag_cell(src, qjoint, live, vc)
            a = _acc_cell(src, qjoint, live, vc, aw_radius)
            product_cells += 1
            values = (
                m["S_innovation_lambda_interval"][1], m["K_theta_operator_norm_upper"],
                m["d_eff_norm_upper"], m["attitude_injection_norm_upper"],
                a["S_innovation_lambda_interval"][1], a["K_theta_operator_norm_upper"],
                a["effective_aw_input_norm_upper_mps2"],
            )
            finite = finite and all(math.isfinite(float(x)) and float(x) >= 0.0 for x in values)
            witness = {
                "source_cell": {k: src[k] for k in ("index", "tau_s", "sigma_aw_mps2", "R_S_filter_std")},
                "q_interval": list(qjoint["q_interval"]),
            }
            worst_mag_d = _worse(worst_mag_d, m["attitude_injection_norm_upper"], {**witness, "mag_cell": m})
            worst_acc_e = _worse(worst_acc_e, a["effective_aw_input_norm_upper_mps2"], {**witness, "acc_cell": a})
            info = float(m["effective_tangent_information_per_cperp2_lower"])
            if min_mag_info is None or info < float(min_mag_info["value"]):
                min_mag_info = {"value": info, **witness, "mag_cell": m}

    if not source_cells or not qcells or product_cells <= 0:
        failures.append("joint prefix cell partition is empty")
    if not finite:
        failures.append("joint prefix cell arithmetic emitted nonfinite bound")
    if min_mag_info is None or not float(min_mag_info["value"]) > 0.0:
        failures.append("tangent-only magnetometer information lost strict positivity")

    # The current P3 source cells expose directional Loewner envelopes only.
    # They do not expose a signed interval matrix P, hence K*r direction cannot
    # be formed without inventing correlations.  Refuse to call a norm-only
    # replacement a signed Cayley cell.
    full_signed_matrix_cells_available = False
    signed_dot_replaced_by_independent_norm_product = False
    first_unclosed = "FULL_MATRIX_P_H_R_K_R_DIRECTION_CELL_NOT_YET_PROPAGATED"

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_JOINT_SOURCE_PREFIX_CELL_ENCLOSURE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "outward_rounded_scalar_arithmetic": True,
        "joint_source_partition_reused_from_P3": True,
        "independent_global_extrema_product_used": False,
        "standalone_vector_eta_penalty_used": False,
        "magnetometer_radial_K_action_exact_zero": True,
        "magnetometer_radial_Joseph_information_exact_zero": True,
        "accelerometer_effective_aw_input_used": True,
        "full_S_to_attitude_gain_retained": True,
        "source_cell_partition": {k: meta[k] for k in ("x_cells", "sigma_cells", "R_S_cells", "joint_source_cells")},
        "attitude_q_cell_count": len(qcells),
        "joint_prefix_product_cell_count": product_cells,
        "word_samples_upper": int(math.ceil(float(domain["normal_live"]["vector_pe_recurrence_window_s"]) / float(meta["source_schedule"]["dt_s"]))) + 2,
        "first_due_S_seed": {
            "widened_cayley_norm_upper": float(first["widened_prefix_cayley_norm_upper"]),
            "first_due_S_state_norm_upper": S_radius,
            "exact_first_S_prefix_status": first["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
        },
        "numerical_extrema": {
            "worst_magnetometer_attitude_injection": worst_mag_d,
            "worst_accelerometer_effective_aw_input": worst_acc_e,
            "worst_S_attitude_injection": worst_S_d,
            "minimum_tangent_only_magnetometer_information": min_mag_info,
        },
        "P_H_R_K_S_r_d_eff_scalar_cell_enclosure_complete": not failures,
        "full_signed_matrix_covariance_cells_available": full_signed_matrix_cells_available,
        "signed_a_dot_c_replaced_by_independent_abs_product": signed_dot_replaced_by_independent_norm_product,
        "signed_cayley_prefix_composition_closed": False,
        "P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE": "NOT_ESTABLISHED",
        "first_unclosed_numerical_obligation": first_unclosed,
        "next_obligation": (
            "replace each directional P envelope by a forward-propagated outward full matrix covariance cell from the exact goLive/first-S source state; compute H,R,S,K,r,d_eff in that same cell, then feed the resulting signed K*r correction vector directly to ou3_p5_signed_cayley_cell without an abs-norm denominator"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "outward_rounded_scalar_arithmetic",
        "joint_source_partition_reused_from_P3", "magnetometer_radial_K_action_exact_zero",
        "magnetometer_radial_Joseph_information_exact_zero", "accelerometer_effective_aw_input_used",
        "full_S_to_attitude_gain_retained", "P_H_R_K_S_r_d_eff_scalar_cell_enclosure_complete",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "independent_global_extrema_product_used",
        "standalone_vector_eta_penalty_used", "signed_a_dot_c_replaced_by_independent_abs_product",
        "full_signed_matrix_covariance_cells_available", "signed_cayley_prefix_composition_closed",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    part = d.get("source_cell_partition", {})
    if int(part.get("joint_source_cells", 0)) <= 0 or int(d.get("joint_prefix_product_cell_count", 0)) <= 0:
        failures.append("joint source-cell counts are empty")
    ext = d.get("numerical_extrema", {})
    mi = ext.get("minimum_tangent_only_magnetometer_information", {})
    if not (isinstance(mi.get("value"), (int, float)) and math.isfinite(float(mi["value"])) and float(mi["value"]) > 0.0):
        failures.append("minimum tangent-only magnetometer information is not strict")
    if d.get("P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("full matrix prefix certificate promoted without signed P/K*r cells")
    if d.get("first_unclosed_numerical_obligation") != "FULL_MATRIX_P_H_R_K_R_DIRECTION_CELL_NOT_YET_PROPAGATED":
        failures.append("unexpected next joint-prefix obligation")
    if not failures and d.get("P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE") != "PASS":
        failures.append("joint scalar cell certificate did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "scalar_cells": out["P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE"],
        "full_matrix": out["P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE"],
        "partition": out["source_cell_partition"],
        "product_cells": out["joint_prefix_product_cell_count"],
        "extrema": out["numerical_extrema"],
        "next": out["first_unclosed_numerical_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
