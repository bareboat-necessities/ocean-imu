#!/usr/bin/env python3
"""Fast value-only companion to the non-promoting physical P4 point probe.

The slow reference evaluator uses outward AD for both values and full 18/21 by
18/21 Jacobians at every event and repeats source-only Joseph solves for every
finite scale.  That is appropriate for universal enclosure development but is
wasteful for a point feasibility sweep whose source P/H/R/K history is fixed.

This companion preserves the exact same point source payload and metric path:

* :mod:`ou3_p4_physical_finite_map_feasibility` still builds the reset-normalized
  Riccati path and performs one complete zero-state outward-AD parity pass;
* every source-only Joseph gain K and every P_after^-1 are precomputed once;
* finite values use the same deployed quaternion branch, exact Cayley
  composition, exact physical S/vector/accelerometer residuals, and the same
  A21 bias projection at the declared radius (currently 0.4 m/s^2);
* every due S update therefore still uses its actual applied anisotropic R_S.

No theorem enclosure is obtained here.  This remains a point
feasibility/falsification diagnostic and cannot promote P4.
The reset-deleted covariance path is a rebuilt diagnostic path, not an
attached finite gauge representative; see ou3_p4_retained_word_attachment.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import ou3_p4_physical_finite_map_feasibility as BASE


def _skew(v: np.ndarray) -> np.ndarray:
    x, y, z = map(float, v)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]], dtype=float)


def _rotation_from_cayley(c: np.ndarray) -> np.ndarray:
    c = np.asarray(c, dtype=float)
    C = _skew(c)
    den = 4.0 + float(c @ c)
    return np.eye(3) + (4.0 / den) * C + (2.0 / den) * (C @ C)


def _rotation_to_cayley(R: np.ndarray) -> np.ndarray:
    R = np.asarray(R, dtype=float)
    den = 1.0 + float(np.trace(R))
    if not den > 1.0e-14:
        raise RuntimeError("finite point map reached Cayley antipode")
    vee = np.array(
        [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]],
        dtype=float,
    )
    return 2.0 * vee / den


def _deployed_quaternion_rotation(d: np.ndarray) -> np.ndarray:
    """Point-value shipping quat_from_delta_theta followed by normalization."""
    d = np.asarray(d, dtype=float)
    theta2 = float(d @ d)
    theta = math.sqrt(max(0.0, theta2))
    if theta < 1.0e-2:
        theta4 = theta2 * theta2
        w = 1.0 - theta2 / 8.0 + theta4 / 384.0
        k = 0.5 - theta2 / 48.0 + theta4 / 3840.0
    else:
        w = math.cos(0.5 * theta)
        k = math.sin(0.5 * theta) / theta if theta > 0.0 else 0.5
    q = np.concatenate(([w], k * d)).astype(float)
    qn = float(np.linalg.norm(q))
    if not qn > 0.0:
        raise RuntimeError("deployed correction quaternion lost norm")
    q /= qn
    w, x, y, z = map(float, q)
    return np.array(
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
            [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
            [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)],
        ],
        dtype=float,
    )


def _prediction_value(ev: dict, state: np.ndarray, n: int) -> np.ndarray:
    # All non-attitude mean channels are exactly linear in shipping.  Reuse the
    # captured source F for bg/v/p/S/a_w/b_a and replace only the finite attitude
    # row by the exact relative-rotation composition.
    out = ev["linear_shipping"][:n, :n] @ state
    E = _rotation_from_cayley(state[:3])
    dbg = state[BASE.OFF_BG:BASE.OFF_BG + 3]
    d_shadow = (-ev["omega"] + dbg) * ev["h"]
    d_nom = -ev["omega"] * ev["h"]
    Rs = _deployed_quaternion_rotation(d_shadow)
    Rn = _deployed_quaternion_rotation(d_nom)
    out[:3] = _rotation_to_cayley(Rs @ E @ Rn.T)
    return out


def _physical_residual(ev: dict, state: np.ndarray, n: int) -> np.ndarray:
    geom = BASE._event_geometry(ev, n)
    if geom["kind"] == "S_zero":
        return state[BASE.OFF_S:BASE.OFF_S + 3].copy()
    E = _rotation_from_cayley(state[:3])
    if geom["kind"] == "vector":
        m = geom["m_body"]
        return (E - np.eye(3)) @ m
    f = geom["f_hat"]
    Rhat = geom["R_hat"]
    y = (E - np.eye(3)) @ f + E @ Rhat @ state[BASE.OFF_AW:BASE.OFF_AW + 3]
    if n == 21:
        y = y + state[BASE.OFF_BA:BASE.OFF_BA + 3]
    return y


def _joseph_value(ev: dict, state: np.ndarray, n: int, K: np.ndarray,
                  projection_limit: float) -> tuple[np.ndarray, str]:
    y = _physical_residual(ev, state, n)
    d = K @ y
    out = state - d
    # source_joseph_event uses deployed_correct_cayley_right(c,-dtheta), i.e.
    # E_plus=E Q(dtheta)^-1.
    E = _rotation_from_cayley(state[:3])
    out[:3] = _rotation_to_cayley(E @ _deployed_quaternion_rotation(-d[:3]))

    branch = "not_applicable"
    if n == 21:
        b = out[BASE.OFF_BA:BASE.OFF_BA + 3]
        r = float(np.linalg.norm(b))
        if r < projection_limit:
            branch = "inactive"
        elif r > projection_limit:
            out[BASE.OFF_BA:BASE.OFF_BA + 3] = (projection_limit / r) * b
            branch = "active"
        else:
            branch = "clarke_hull"
    return out, branch


def prepare_kernels(payload: dict, linear: dict) -> list[dict]:
    """Precompute scale-independent gain and energy kernels exactly once."""
    n = payload["mode_dim"]
    kernels = []
    for ev, lp in zip(payload["events"], linear["path"]):
        row = {"Qafter": BASE._spd_inverse(lp["Pafter"]), "K": None}
        if ev["type"] not in (BASE.EV_PRED, BASE.EV_FLOOR):
            H = ev["H"][:, :n]
            S = BASE._sym(H @ lp["Pbefore"] @ H.T + ev["R"])
            row["K"] = np.linalg.solve(S, H @ lp["Pbefore"]).T
        kernels.append(row)
    return kernels


def propagate_case(payload: dict, linear: dict, kernels: list[dict], domain: dict,
                   input_path: Path, projection_limit: float, scale: float) -> dict:
    n = payload["mode_dim"]
    state = scale * linear["direction"]
    V0 = float(state @ linear["Q0"] @ state)
    if not (V0 > 0.0 and math.isfinite(V0)):
        raise RuntimeError("invalid point initial energy")
    hs = BASE.parse_hs(input_path)
    max_ratio = 1.0
    first_domain_failure = None
    projection_counts = {"inactive": 0, "active": 0, "clarke_hull": 0, "not_applicable": 0}

    for idx, (ev, kernel) in enumerate(zip(payload["events"], kernels)):
        if ev["type"] == BASE.EV_PRED:
            state = _prediction_value(ev, state, n)
        elif ev["type"] == BASE.EV_FLOOR:
            state = state.copy()
        else:
            K = kernel["K"]
            if K is None:
                raise RuntimeError("measurement event lost precomputed Joseph gain")
            state, branch = _joseph_value(ev, state, n, K, projection_limit)
            projection_counts[branch] = projection_counts.get(branch, 0) + 1

        V = float(state @ kernel["Qafter"] @ state)
        max_ratio = max(max_ratio, V / V0)
        if first_domain_failure is None:
            ok, why = BASE.domain_status(state, n, domain, hs)
            if not ok:
                first_domain_failure = {
                    "event_index": idx, "event": ev["name"], "reason": why, "time": ev["time"]
                }

    Vend = float(state @ linear["QN"] @ state)
    return {
        "scale": float(scale), "absolute_scale": abs(float(scale)),
        "V0": V0, "VN": Vend, "rho_endpoint": Vend / V0,
        "max_prefix_ratio": max_ratio,
        "domain_retained": first_domain_failure is None,
        "first_domain_failure": first_domain_failure,
        "max_interval_width": 0.0,
        "projection_counts": projection_counts,
        "endpoint_group_norms": {
            "theta": BASE.group_norm(state, 0, 3), "bg": BASE.group_norm(state, 3, 6),
            "v": BASE.group_norm(state, 6, 9), "p": BASE.group_norm(state, 9, 12),
            "S": BASE.group_norm(state, 12, 15), "aw": BASE.group_norm(state, 15, 18),
            "ba": BASE.group_norm(state, 18, 21) if n == 21 else 0.0,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload", type=Path, required=True)
    ap.add_argument("--domain", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    payload = BASE.read_payload(args.payload)
    domain = json.loads(args.domain.read_text(encoding="utf-8"))
    projection_limit = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    linear = BASE.build_reset_normalized_linear_path(payload, projection_limit)

    # Keep the expensive theorem-facing AD machinery as a once-per-word parity
    # test.  Only the repeated finite scale sweep is value-only.
    parity = BASE.zero_state_parity(payload, linear, projection_limit)
    if not parity["pass"]:
        raise RuntimeError(f"physical finite-map zero-state parity failed: {parity}")
    if linear["minimum_measurement_R_eigenvalue"] <= 0.0:
        raise RuntimeError("physical point payload has non-positive measurement R")
    if linear["counts"]["prediction"] != 600 or linear["counts"]["accelerometer"] != 600:
        raise RuntimeError(f"physical point payload is not a 600-sample complete word: {linear['counts']}")
    if linear["actual_RS_std_ratio_max_error"] > 2e-5:
        raise RuntimeError("physical point payload lost actual R_S anisotropy")

    kernels = prepare_kernels(payload, linear)
    limit, limit_detail = BASE.scale_limit(linear["direction"], payload["mode_dim"], domain, args.input)
    scales = BASE.choose_scales(limit)
    cases = []
    for s in scales:
        cases.append(propagate_case(payload, linear, kernels, domain, args.input, projection_limit, -s))
        cases.append(propagate_case(payload, linear, kernels, domain, args.input, projection_limit, +s))

    valid = [c for c in cases if c["domain_retained"]]
    crossings = [c for c in valid if c["rho_endpoint"] >= 1.0]
    report = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_PHYSICAL_FINITE_MAP_POINT_FAST_V1",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "physical_true_minus_estimated_map": True,
        "reset_gauge_normalized_representative": False,
        "same_single_shipping_observer_source_payload": True,
        "second_estimator_or_Riccati_history_used": True,
        "diagnostic_Riccati_rebuilt_without_resets": True,
        "finite_gauge_attachment_verified": False,
        "source_only_Joseph_gains_precomputed_once": True,
        "finite_scale_sweep_value_only_same_deployed_equations": True,
        "zero_state_full_AD_parity_still_required": True,
        "all_due_S_updates_with_actual_RS_retained": True,
        "packet_count_remainder_budget_used": False,
        "state_elimination_used": False,
        "trajectory_replay_promoted_to_theorem_source": False,
        "P4_promoted": False,
        "P5_may_start": False,
        "mode": "H18" if payload["mode_dim"] == 18 else "A21",
        "dimension": payload["mode_dim"],
        "word_t0": payload["t0"], "word_t1": payload["t1"],
        "event_counts": linear["counts"],
        "reset_normalized_linear_rho": linear["rho_linear"],
        "reset_normalized_linear_distance_to_one": 1.0 - linear["rho_linear"],
        "zero_state_parity": parity,
        "actual_RS_std_ratio_max_error": linear["actual_RS_std_ratio_max_error"],
        "minimum_measurement_R_eigenvalue": linear["minimum_measurement_R_eigenvalue"],
        "direction": linear["direction"].tolist(),
        "direction_information_energy": float(linear["direction"] @ linear["Q0"] @ linear["direction"]),
        "declared_scale_limit": limit,
        "declared_scale_limit_detail": limit_detail,
        "tested_absolute_scales": scales,
        "cases": cases,
        "valid_domain_cases": len(valid),
        "strict_endpoint_contraction_on_all_domain_retained_tested_cases": bool(valid) and not crossings,
        "first_endpoint_rho_ge_one_absolute_scale": min((c["absolute_scale"] for c in crossings), default=None),
        "worst_endpoint_rho_on_domain_retained_cases": max((c["rho_endpoint"] for c in valid), default=None),
        "worst_prefix_ratio_on_domain_retained_cases": max((c["max_prefix_ratio"] for c in valid), default=None),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "mode": report["mode"],
        "linear_rho": report["reset_normalized_linear_rho"],
        "zero_parity": report["zero_state_parity"],
        "scale_limit": report["declared_scale_limit"],
        "valid_cases": report["valid_domain_cases"],
        "worst_endpoint_rho": report["worst_endpoint_rho_on_domain_retained_cases"],
        "first_endpoint_crossing": report["first_endpoint_rho_ge_one_absolute_scale"],
        "worst_prefix_ratio": report["worst_prefix_ratio_on_domain_retained_cases"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
