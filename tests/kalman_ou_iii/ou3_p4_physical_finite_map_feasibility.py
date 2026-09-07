#!/usr/bin/env python3
"""Non-promoting physical finite-map feasibility on one complete SEA3 point word.

The source/event payload is emitted by the same shipping observer used by the
linear complete-word scan.  This evaluator does not run another estimator,
scheduler, tuner or Riccati recursion.  It chooses the reset-gauge-normalized
representative justified by the P3 reset-congruence bridge: prediction/floor
source data and every H/R measurement cell are retained, while covariance
resets are represented by G=I.  The resulting Riccati path is rebuilt from the
same source F/Q/H/R sequence and every due S operation retains the actual
applied R_S from that source word.

The finite state is the paper's true-minus-estimated physical Cayley error.  It
is propagated by the existing theorem-facing physical prediction and Joseph
event functions.  A21 uses homogeneous source true residual bias zero and the
projection radius read from the current proof domain; the full delta-b_a error
coordinate is retained.

Before any finite ratio is reported, the zero-state physical Jacobian cocycle
must reproduce the reset-normalized linear cocycle.  This is a point
feasibility/falsification diagnostic only.  It never promotes P4 or substitutes
for source-uniform complete-SEA3 enclosure.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[2]
STABILITY = REPO / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_point
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_prediction as PRED

NX = 21
OFF_BG = 3
OFF_V = 6
OFF_P = 9
OFF_S = 12
OFF_AW = 15
OFF_BA = 18

EV_PRED = 1
EV_FLOOR = 2
EV_S = 3
EV_ACC = 4
EV_VECTOR = 5
EVENT_NAMES = {
    EV_PRED: "prediction",
    EV_FLOOR: "aw_floor",
    EV_S: "S_zero",
    EV_ACC: "accelerometer",
    EV_VECTOR: "vector",
}

HEADER = struct.Struct("<8sIIIdd")
REC_PREFIX = struct.Struct("<IIdfffff")
REC_FLOATS = 6 + 4 * NX * NX + 3 * NX + 9 + 9
REC_BYTES = REC_PREFIX.size + 4 * REC_FLOATS


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _mid(x: Interval) -> float:
    return 0.5 * (float(x.lo) + float(x.hi))


def _max_width(xs) -> float:
    return max((float(x.hi) - float(x.lo) for x in xs), default=0.0)


def _mat_mid(A) -> np.ndarray:
    return np.asarray([[_mid(x) for x in row] for row in A], dtype=float)


def _sym(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def _spd_inverse(P: np.ndarray) -> np.ndarray:
    P = _sym(P)
    ev = np.linalg.eigvalsh(P)
    if not (np.all(np.isfinite(ev)) and ev[0] > 0.0):
        raise RuntimeError(f"point covariance is not SPD; lambda_min={ev[0] if len(ev) else None}")
    return _sym(np.linalg.inv(P))


def _read_matrix(vals: np.ndarray, offset: int, rows: int, cols: int):
    count = rows * cols
    return vals[offset:offset + count].reshape(rows, cols).astype(float), offset + count


def read_payload(path: Path) -> dict:
    data = path.read_bytes()
    if len(data) < HEADER.size:
        raise RuntimeError("physical payload is truncated")
    magic, version, nx, mode_dim, t0, t1 = HEADER.unpack_from(data, 0)
    if magic != b"OU3PHY1\0" or version != 1 or nx != NX or mode_dim not in (18, 21):
        raise RuntimeError("unexpected physical payload schema")
    tail = len(data) - HEADER.size
    if tail % REC_BYTES != 0:
        raise RuntimeError(f"physical payload record alignment failed: tail={tail} rec={REC_BYTES}")
    events = []
    off = HEADER.size
    while off < len(data):
        event_type, event_dim, t, h, tau_aw, tau_ba, rs_scalar, pseudo = REC_PREFIX.unpack_from(data, off)
        off += REC_PREFIX.size
        if event_type not in EVENT_NAMES or event_dim != mode_dim:
            raise RuntimeError("physical payload event type/dimension drifted")
        vals = np.frombuffer(data, dtype="<f4", count=REC_FLOATS, offset=off).astype(float)
        off += 4 * REC_FLOATS
        k = 0
        omega, k = _read_matrix(vals, k, 3, 1)
        dtheta, k = _read_matrix(vals, k, 3, 1)
        Pbefore, k = _read_matrix(vals, k, NX, NX)
        Pafter, k = _read_matrix(vals, k, NX, NX)
        linear, k = _read_matrix(vals, k, NX, NX)
        Q, k = _read_matrix(vals, k, NX, NX)
        H, k = _read_matrix(vals, k, 3, NX)
        R, k = _read_matrix(vals, k, 3, 3)
        aux, k = _read_matrix(vals, k, 3, 3)
        if k != REC_FLOATS:
            raise RuntimeError("physical payload parser offset drifted")
        events.append({
            "type": int(event_type), "name": EVENT_NAMES[int(event_type)],
            "time": float(t), "h": float(h), "tau_aw": float(tau_aw),
            "tau_ba": float(tau_ba), "rs_scalar": float(rs_scalar),
            "pseudo_period": float(pseudo), "omega": omega[:, 0],
            "dtheta": dtheta[:, 0], "Pbefore_shipping": _sym(Pbefore),
            "Pafter_shipping": _sym(Pafter), "linear_shipping": linear,
            "Q": _sym(Q), "H": H, "R": _sym(R), "aux": _sym(aux),
        })
    if not events or events[0]["type"] != EV_PRED:
        raise RuntimeError("selected physical payload must begin with prediction")
    return {"mode_dim": int(mode_dim), "t0": float(t0), "t1": float(t1), "events": events}


def _positive_part(A: np.ndarray) -> np.ndarray:
    A = _sym(A)
    w, V = np.linalg.eigh(A)
    return _sym(V @ np.diag(np.maximum(w, 0.0)) @ V.T)


def _joseph(P: np.ndarray, H: np.ndarray, R: np.ndarray):
    S = _sym(H @ P @ H.T + R)
    K = np.linalg.solve(S, H @ P).T
    A = np.eye(P.shape[0]) - K @ H
    Pn = _sym(A @ P @ A.T + K @ R @ K.T)
    return Pn, K, A


def _event_geometry(event: dict, n: int):
    H = event["H"][:, :n]
    if event["type"] == EV_S:
        return {"kind": "S_zero"}
    if event["type"] == EV_ACC:
        f = np.array([H[1, 2], H[2, 0], H[0, 1]], dtype=float)
        Rhat = H[:, OFF_AW:OFF_AW + 3].copy()
        return {"kind": "accelerometer", "f_hat": f, "R_hat": Rhat}
    if event["type"] == EV_VECTOR:
        m = np.array([H[1, 2], H[2, 0], H[0, 1]], dtype=float)
        return {"kind": "vector", "m_body": m}
    raise ValueError("geometry requested for non-measurement event")


def _point_matrix(A: np.ndarray):
    return matrix_point([[float(x) for x in row] for row in A])


def _point_vector(x: np.ndarray):
    return [I(float(v)) for v in x]


def _generated_H(event: dict, n: int, P: np.ndarray, projection_limit: float):
    zero = [I(0.0) for _ in range(n)]
    geom = _event_geometry(event, n)
    kwargs = {}
    if geom["kind"] == "accelerometer":
        kwargs.update(
            f_hat=_point_vector(geom["f_hat"]),
            R_hat=_point_matrix(geom["R_hat"]),
        )
    elif geom["kind"] == "vector":
        kwargs.update(m_body=_point_vector(geom["m_body"]))
    if geom["kind"] == "S_zero":
        kwargs["R_provenance"] = EVENTS.ACTUAL_RS_PROVENANCE
    if n == 21:
        kwargs.update(
            bias_true=[I(0.0), I(0.0), I(0.0)],
            bias_projection_limit=projection_limit,
        )
    out = EVENTS.source_joseph_event(
        "H" if n == 18 else "A", zero, _point_matrix(P), _point_matrix(event["R"]),
        geom["kind"], **kwargs,
    )
    return out


def build_reset_normalized_linear_path(payload: dict, projection_limit: float) -> dict:
    n = payload["mode_dim"]
    events = payload["events"]
    P = events[0]["Pbefore_shipping"][:n, :n].copy()
    P0 = P.copy()
    M = np.eye(n)
    path = []
    counts = {name: 0 for name in EVENT_NAMES.values()}
    rs_ratio_max_error = 0.0
    min_R_eig = math.inf

    for ev in events:
        Pbefore = P.copy()
        if ev["type"] == EV_PRED:
            A = ev["linear_shipping"][:n, :n]
            Q = ev["Q"][:n, :n]
            P = _sym(A @ P @ A.T + Q)
            C = A
        elif ev["type"] == EV_FLOOR:
            target = ev["aux"]
            D = _positive_part(target - P[OFF_AW:OFF_AW + 3, OFF_AW:OFF_AW + 3])
            P = P.copy()
            P[OFF_AW:OFF_AW + 3, OFF_AW:OFF_AW + 3] += D
            P = _sym(P)
            C = np.eye(n)
        else:
            H = ev["H"][:, :n]
            R = ev["R"]
            min_R_eig = min(min_R_eig, float(np.linalg.eigvalsh(R)[0]))
            P, _K, C = _joseph(P, H, R)
            if ev["type"] == EV_S:
                std = np.sqrt(np.maximum(np.diag(R), 0.0))
                if std[2] > 0:
                    rs_ratio_max_error = max(
                        rs_ratio_max_error,
                        abs(float(std[0] / std[2]) - 0.72),
                        abs(float(std[1] / std[2]) - 0.72),
                    )
        M = C @ M
        counts[ev["name"]] += 1
        path.append({"Pbefore": Pbefore, "Pafter": P.copy(), "C": C.copy()})

    Q0 = _spd_inverse(P0)
    QN = _spd_inverse(P)
    Aend = _sym(M.T @ QN @ M)
    L = np.linalg.cholesky(Q0)
    Linv = np.linalg.inv(L)
    B = _sym(Linv @ Aend @ Linv.T)
    w, V = np.linalg.eigh(B)
    y = V[:, -1]
    direction = np.linalg.solve(L.T, y)
    energy = float(direction @ Q0 @ direction)
    direction /= math.sqrt(energy)
    rho = float(w[-1])
    return {
        "P0": P0, "PN": P, "Q0": Q0, "QN": QN, "M": M,
        "path": path, "rho_linear": rho, "direction": direction,
        "counts": counts, "actual_RS_std_ratio_max_error": rs_ratio_max_error,
        "minimum_measurement_R_eigenvalue": min_R_eig,
    }


def zero_state_parity(payload: dict, linear: dict, projection_limit: float) -> dict:
    n = payload["mode_dim"]
    mode = "H" if n == 18 else "A"
    Jprod = np.eye(n)
    max_event_rel = 0.0
    max_H_rel = 0.0
    max_width = 0.0
    projection_branches = set()

    for ev, lp in zip(payload["events"], linear["path"]):
        if ev["type"] == EV_PRED:
            out = PRED.prediction_event(
                mode, [I(0.0) for _ in range(n)], _point_vector(ev["omega"]),
                I(ev["h"]), I(ev["tau_aw"]),
                tau_ba=I(ev["tau_ba"]) if n == 21 else None,
            )
            J = _mat_mid(out["J_state"])
            max_width = max(max_width, max(_max_width(row) for row in out["J_state"]))
        elif ev["type"] == EV_FLOOR:
            J = np.eye(n)
        else:
            out = _generated_H(ev, n, lp["Pbefore"], projection_limit)
            J = _mat_mid(out["J_state"])
            Hgen = _mat_mid(out["H"])
            Hsrc = ev["H"][:, :n]
            max_H_rel = max(max_H_rel, float(np.linalg.norm(Hgen-Hsrc) / max(1.0, np.linalg.norm(Hsrc))))
            max_width = max(max_width, max(_max_width(row) for row in out["J_state"]))
            projection_branches.add(out["bias_projection_branch"])
        C = lp["C"]
        max_event_rel = max(max_event_rel, float(np.linalg.norm(J-C) / max(1.0, np.linalg.norm(C))))
        Jprod = J @ Jprod

    word_rel = float(np.linalg.norm(Jprod-linear["M"]) / max(1.0, np.linalg.norm(linear["M"])))
    return {
        "max_event_relative_difference": max_event_rel,
        "word_relative_difference": word_rel,
        "max_generated_H_relative_difference": max_H_rel,
        "max_interval_width": max_width,
        "projection_branches": sorted(projection_branches),
        "pass": max_event_rel <= 2.0e-4 and word_rel <= 3.0e-4 and max_H_rel <= 2.0e-5,
    }


def parse_hs(path: Path) -> float | None:
    import re
    m = re.search(r"_H([0-9]+(?:\.[0-9]+)?)_", path.name)
    return float(m.group(1)) if m else None


def group_norm(x: np.ndarray, a: int, b: int) -> float:
    return float(np.linalg.norm(x[a:b]))


def scale_limit(direction: np.ndarray, n: int, domain: dict, input_path: Path) -> tuple[float, dict]:
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    angle_deg = max(float(x) for x in domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"])
    cayley_cap = 2.0 * math.tan(0.5 * math.radians(angle_deg))
    caps = {
        "theta_cayley": (group_norm(direction, 0, 3), cayley_cap),
        "bg": (group_norm(direction, 3, 6), float(handoff["gyro_bias_error_norm_upper_rad_s"])),
        "v": (group_norm(direction, 6, 9), float(handoff["velocity_error_norm_upper_mps"])),
        "p": (group_norm(direction, 9, 12), float(handoff["position_error_norm_upper_m"])),
        "S": (group_norm(direction, 12, 15), float(handoff["integral_displacement_error_norm_upper_m_s"])),
        "aw": (group_norm(direction, 15, 18), float(handoff["latent_acceleration_error_norm_upper_mps2"])),
    }
    if n == 21:
        caps["ba"] = (group_norm(direction, 18, 21), float(handoff["accelerometer_bias_error_norm_upper_mps2"]))
    limits = {k: cap/g for k,(g,cap) in caps.items() if g > 0.0}
    hs = parse_hs(input_path)
    if hs is not None:
        for j, x in enumerate(direction[9:12]):
            if abs(float(x)) > 0.0:
                limits[f"p_component_{j}"] = 0.5 * hs / abs(float(x))
    limiter = min(limits, key=limits.get)
    return float(limits[limiter]), {
        "limiting_constraint": limiter, "all_limits": limits,
        "attitude_candidate_deg": angle_deg, "cayley_norm_cap": cayley_cap, "Hs_m": hs,
    }


def domain_status(x: np.ndarray, n: int, domain: dict, hs: float | None) -> tuple[bool, str | None]:
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    angle_deg = max(float(v) for v in domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"])
    cayley_cap = 2.0 * math.tan(0.5 * math.radians(angle_deg))
    checks = [
        (group_norm(x,0,3) <= cayley_cap*(1+1e-9), "theta"),
        (group_norm(x,3,6) <= float(handoff["gyro_bias_error_norm_upper_rad_s"])*(1+1e-9), "bg"),
        (group_norm(x,6,9) <= float(handoff["velocity_error_norm_upper_mps"])*(1+1e-9), "v"),
        (group_norm(x,9,12) <= float(handoff["position_error_norm_upper_m"])*(1+1e-9), "p"),
        (group_norm(x,12,15) <= float(handoff["integral_displacement_error_norm_upper_m_s"])*(1+1e-9), "S"),
        (group_norm(x,15,18) <= float(handoff["latent_acceleration_error_norm_upper_mps2"])*(1+1e-9), "aw"),
    ]
    if n == 21:
        checks.append((group_norm(x,18,21) <= float(handoff["accelerometer_bias_error_norm_upper_mps2"])*(1+1e-9), "ba"))
    if hs is not None:
        checks.append((bool(np.all(np.abs(x[9:12]) <= 0.5*hs*(1+1e-9))), "p_component_Hs"))
    for ok, name in checks:
        if not ok:
            return False, name
    return True, None


def propagate_case(payload: dict, linear: dict, domain: dict, input_path: Path,
                   projection_limit: float, scale: float) -> dict:
    n = payload["mode_dim"]
    mode = "H" if n == 18 else "A"
    state = scale * linear["direction"]
    V0 = float(state @ linear["Q0"] @ state)
    if not (V0 > 0.0 and math.isfinite(V0)):
        raise RuntimeError("invalid point initial energy")
    hs = parse_hs(input_path)
    max_ratio = 1.0
    first_domain_failure = None
    max_interval_width = 0.0
    projection_counts = {"inactive":0, "active":0, "clarke_hull":0, "not_applicable":0}

    for idx, (ev, lp) in enumerate(zip(payload["events"], linear["path"])):
        if ev["type"] == EV_PRED:
            out = PRED.prediction_event(
                mode, _point_vector(state), _point_vector(ev["omega"]), I(ev["h"]), I(ev["tau_aw"]),
                tau_ba=I(ev["tau_ba"]) if n == 21 else None,
            )
            state_i = out["state_out"]
        elif ev["type"] == EV_FLOOR:
            state_i = _point_vector(state)
        else:
            geom = _event_geometry(ev, n)
            kwargs = {}
            if geom["kind"] == "accelerometer":
                kwargs.update(f_hat=_point_vector(geom["f_hat"]), R_hat=_point_matrix(geom["R_hat"]))
            elif geom["kind"] == "vector":
                kwargs.update(m_body=_point_vector(geom["m_body"]))
            if geom["kind"] == "S_zero":
                kwargs["R_provenance"] = EVENTS.ACTUAL_RS_PROVENANCE
            if n == 21:
                kwargs.update(
                    bias_true=[I(0.0), I(0.0), I(0.0)],
                    bias_projection_limit=projection_limit,
                )
            out = EVENTS.source_joseph_event(
                mode, _point_vector(state), _point_matrix(lp["Pbefore"]), _point_matrix(ev["R"]),
                geom["kind"], **kwargs,
            )
            state_i = out["state_out"]
            projection_counts[out["bias_projection_branch"]] = projection_counts.get(out["bias_projection_branch"],0)+1
        max_interval_width = max(max_interval_width, _max_width(state_i))
        state = np.asarray([_mid(x) for x in state_i], dtype=float)
        Qp = _spd_inverse(lp["Pafter"])
        V = float(state @ Qp @ state)
        max_ratio = max(max_ratio, V/V0)
        if first_domain_failure is None:
            ok, why = domain_status(state, n, domain, hs)
            if not ok:
                first_domain_failure = {"event_index": idx, "event": ev["name"], "reason": why, "time": ev["time"]}

    Vend = float(state @ linear["QN"] @ state)
    return {
        "scale": float(scale), "absolute_scale": abs(float(scale)),
        "V0": V0, "VN": Vend, "rho_endpoint": Vend/V0,
        "max_prefix_ratio": max_ratio,
        "domain_retained": first_domain_failure is None,
        "first_domain_failure": first_domain_failure,
        "max_interval_width": max_interval_width,
        "projection_counts": projection_counts,
        "endpoint_group_norms": {
            "theta": group_norm(state,0,3), "bg": group_norm(state,3,6),
            "v": group_norm(state,6,9), "p": group_norm(state,9,12),
            "S": group_norm(state,12,15), "aw": group_norm(state,15,18),
            "ba": group_norm(state,18,21) if n==21 else 0.0,
        },
    }


def choose_scales(limit: float) -> list[float]:
    base = [0.125, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 16.0, 32.0]
    xs = [x for x in base if x <= limit*(1+1e-12)]
    if not xs:
        xs = [0.125*limit, 0.25*limit, 0.5*limit]
    if limit > 0 and (not xs or abs(xs[-1]-limit) > 1e-9*max(1.0,limit)):
        xs.append(limit)
    return sorted(set(float(x) for x in xs if x > 0.0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload", type=Path, required=True)
    ap.add_argument("--domain", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    payload = read_payload(args.payload)
    domain = json.loads(args.domain.read_text(encoding="utf-8"))
    projection_limit = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    linear = build_reset_normalized_linear_path(payload, projection_limit)
    parity = zero_state_parity(payload, linear, projection_limit)
    if not parity["pass"]:
        raise RuntimeError(f"physical finite-map zero-state parity failed: {parity}")
    if linear["minimum_measurement_R_eigenvalue"] <= 0.0:
        raise RuntimeError("physical point payload has non-positive measurement R")
    if linear["counts"]["prediction"] != 600 or linear["counts"]["accelerometer"] != 600:
        raise RuntimeError(f"physical point payload is not a 600-sample complete word: {linear['counts']}")
    if linear["actual_RS_std_ratio_max_error"] > 2e-5:
        raise RuntimeError("physical point payload lost actual R_S anisotropy")

    limit, limit_detail = scale_limit(linear["direction"], payload["mode_dim"], domain, args.input)
    scales = choose_scales(limit)
    cases = []
    for s in scales:
        cases.append(propagate_case(payload, linear, domain, args.input, projection_limit, -s))
        cases.append(propagate_case(payload, linear, domain, args.input, projection_limit, +s))

    valid = [c for c in cases if c["domain_retained"]]
    endpoint_crossings = [c for c in valid if c["rho_endpoint"] >= 1.0]
    report = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_PHYSICAL_FINITE_MAP_POINT_V1",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "physical_true_minus_estimated_map": True,
        "reset_gauge_normalized_representative": True,
        "same_single_shipping_observer_source_payload": True,
        "second_estimator_or_Riccati_history_used": False,
        "all_due_S_updates_with_actual_RS_retained": True,
        "packet_count_remainder_budget_used": False,
        "state_elimination_used": False,
        "trajectory_replay_promoted_to_theorem_source": False,
        "P4_promoted": False,
        "P5_may_start": False,
        "mode": "H18" if payload["mode_dim"]==18 else "A21",
        "dimension": payload["mode_dim"],
        "word_t0": payload["t0"], "word_t1": payload["t1"],
        "event_counts": linear["counts"],
        "reset_normalized_linear_rho": linear["rho_linear"],
        "reset_normalized_linear_distance_to_one": 1.0-linear["rho_linear"],
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
        "strict_endpoint_contraction_on_all_domain_retained_tested_cases": bool(valid) and not endpoint_crossings,
        "first_endpoint_rho_ge_one_absolute_scale": min((c["absolute_scale"] for c in endpoint_crossings), default=None),
        "worst_endpoint_rho_on_domain_retained_cases": max((c["rho_endpoint"] for c in valid), default=None),
        "worst_prefix_ratio_on_domain_retained_cases": max((c["max_prefix_ratio"] for c in valid), default=None),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
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
