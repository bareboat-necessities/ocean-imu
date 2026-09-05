#!/usr/bin/env python3
"""Stable source-history-free SEA3 finite-word translation lower.

This is the arithmetic refinement of ``ou3_sea3_finite_word_translation``.
The theorem is unchanged: tau may take any value in the complete current SEA3
invariant independently at every 5 ms step.  The global x=h/tau box is split
only into instantaneous arithmetic cells.  No predecessor/source history is
created.

For each step and each x cell we form a validated prediction lower.  A single
common matrix is then obtained by a generalized Loewner comparison against the
lowest-x cell.  The common matrix is valid for every x cell and is the only
state propagated to the next step.  Thus the cell index is forgotten after
every sample (history depth zero), while the full 4x4 integrated-OU geometry is
retained instead of replacing it with the one-step scalar process floor.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_finite_word_translation as BASE
import ou3_sea3_riccati_tube as TUBE_BASE
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_vector_uco_certificate as VECTOR

SCHEMA = 2
QUALIFICATION = "OU3_SEA3_SOURCE_HISTORY_FREE_FINITE_WORD_TRANSLATION_LOWER_V2"
USEFUL_GATE = BASE.USEFUL_GATE
WORD_HORIZON_S = BASE.WORD_HORIZON_S
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN


def _minus_matrix(A, alpha: float, B):
    q = BASE.I(float(alpha))
    return BASE._sym([
        [A[i][j] - q * B[i][j] for j in range(4)]
        for i in range(4)
    ])


def _relative_floor(A, B) -> float:
    """Largest certified alpha in [0,1] with A >= alpha B."""
    if not symmetric_positive_definite_ldlt(A)[0]:
        return 0.0
    if not symmetric_positive_definite_ldlt(B)[0]:
        raise RuntimeError("reference common lower is not SPD")
    if symmetric_positive_definite_ldlt(_minus_matrix(A, 1.0, B))[0]:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(58):
        mid = 0.5 * (lo + hi)
        if symmetric_positive_definite_ldlt(_minus_matrix(A, mid, B))[0]:
            lo = mid
        else:
            hi = mid
    return BASE.down(lo)


def _cell_data(global_x: Interval, sigma2: float):
    # Reuse the canonical factored Q splitter.  The returned leaves cover the
    # whole current source interval and are not connected across time.
    leaves = TUBE.split_x_cell(global_x)
    cells = []
    for x, _rho in leaves:
        F = BASE._transition(x)
        Q = BASE._scale(TUBE.step_scaled_q(x), sigma2)
        cells.append((x, F, BASE._transpose(F), Q))
    if not cells:
        raise RuntimeError("empty instantaneous x cover")
    return cells


def _prediction_lower(P, cell):
    _x, F, Ft, Q = cell
    A = BASE._sym(BASE._add(BASE._mul(BASE._mul(F, P), Ft), Q))
    L, eps, route = BASE._common_point_lower(A)
    return L, eps, route


def _common_prediction(P, cells):
    lowers = []
    max_eps = 0.0
    relative_shaves = 0
    for cell in cells:
        L, eps, route = _prediction_lower(P, cell)
        if not symmetric_positive_definite_ldlt(L)[0]:
            x = cell[0]
            raise RuntimeError(f"instantaneous prediction lower lost SPD on x cell {x.as_list()}")
        lowers.append(L)
        max_eps = max(max_eps, eps)
        relative_shaves += int(route == "RELATIVE_DIAGONAL")

    # Lowest-x cell is the natural low-process reference.  We do not assume it
    # is Loewner-minimal: every other cell is explicitly compared to it.
    ref = lowers[0]
    alpha = 1.0
    limiting_cell = 0
    for i, L in enumerate(lowers[1:], 1):
        a = _relative_floor(L, ref)
        if a < alpha:
            alpha = a
            limiting_cell = i
    if not alpha > 0.0:
        raise RuntimeError("instantaneous SEA3 cells have no positive common prediction lower")
    common = BASE._scale(ref, BASE.down(alpha))
    if not symmetric_positive_definite_ldlt(common)[0]:
        raise RuntimeError("scaled common prediction lower lost SPD")
    return common, max_eps, relative_shaves, alpha, limiting_cell


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("finite-word translation proof may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"SEA3 dynamic source prerequisite failed: {df}")
    tube = BASE._load_tube(path, tube_path)

    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = float(rates["dt_s"])
    n = int(round(WORD_HORIZON_S / h))
    if n < 1 or abs(n * h - WORD_HORIZON_S) > 5.0e-6:
        raise RuntimeError("finite word is not represented by deployed samples")

    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_lo = float(inv["sigma_aw_filter_mps2"][0])
    if not (tau_lo > 0.0 and tau_hi >= tau_lo and sigma_lo > 0.0):
        raise RuntimeError("SEA3 dynamic invariant lost positive translation source")
    global_x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    if not (0.0 < global_x.lo <= global_x.hi < 0.05):
        raise RuntimeError("global SEA3 h/tau interval left audited range")
    cells = _cell_data(global_x, BASE.down(sigma_lo * sigma_lo))

    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    if vf:
        raise RuntimeError(f"vector measurement prerequisite failed: {vf}")
    R_aw = BASE.down(float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"]) ** 2)
    axis_factor = min(map(float, TUBE_BASE._axis_factors()))
    rs_lo = float(inv["R_S_applied"][0])
    rS = BASE.down(rs_lo * axis_factor)
    R_S_z = BASE.down(BASE.down(rS * rS) / BASE.up(h ** 6))
    if not (R_aw > 0.0 and R_S_z > 0.0):
        raise RuntimeError("strongest translation measurement variance lost positivity")

    P = BASE._zero()
    max_prediction_eps = 0.0
    max_measurement_eps = 0.0
    min_common_alpha = 1.0
    limiting_step = 0
    limiting_cell = 0
    relative_prediction_shaves = 0
    relative_measurement_shaves = 0

    for step in range(1, n + 1):
        P, eps, nshave, alpha, icell = _common_prediction(P, cells)
        max_prediction_eps = max(max_prediction_eps, eps)
        relative_prediction_shaves += nshave
        if alpha < min_common_alpha:
            min_common_alpha = alpha
            limiting_step = step
            limiting_cell = icell

        # Applying every possible same-sample direct a_w and S=0 observation,
        # at their smallest declared R, gives a covariance no larger than any
        # shipping subset and is therefore a valid lower recursion.
        P, eps, route = BASE._measurement_lower(P, 3, R_aw)
        if not symmetric_positive_definite_ldlt(P)[0]:
            raise RuntimeError(f"a_w measurement lower lost SPD at step {step}")
        max_measurement_eps = max(max_measurement_eps, eps)
        relative_measurement_shaves += int(route == "RELATIVE_DIAGONAL")

        P, eps, route = BASE._measurement_lower(P, 2, R_S_z)
        if not symmetric_positive_definite_ldlt(P)[0]:
            raise RuntimeError(f"S measurement lower lost SPD at step {step}")
        max_measurement_eps = max(max_measurement_eps, eps)
        relative_measurement_shaves += int(route == "RELATIVE_DIAGONAL")

    pdiag = list(map(float, tube["modes"]["H"]["Pbar_diagonal_variance_upper"]))
    upper = [pdiag[6], pdiag[9], pdiag[12], pdiag[15]]
    dscale = [h, h * h, h * h * h, 1.0]
    upper_z = [BASE.up(upper[i] / BASE.down(dscale[i] * dscale[i])) for i in range(4)]
    rho = BASE._generalized_rho(P, upper_z)
    gate_pass = rho >= USEFUL_GATE

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_history_frontier_consumed": False,
        "SEA3_dynamic_source_consumed": True,
        "instantaneous_source_cell_partition_consumed": True,
        "history_depth": 0,
        "instantaneous_x_cells": len(cells),
        "global_source_interval_reselected_independently_each_step": True,
        "cell_index_forgotten_after_each_step": True,
        "tau_path_correlation_assumed": False,
        "sigma_path_correlation_assumed": False,
        "R_S_path_correlation_assumed": False,
        "full_4x4_translation_matrix_retained": True,
        "determinant_or_e3_scalarization_used": False,
        "fixed_physical_scaling": "z=[v/h,p/h^2,S/h^3,a_w]",
        "strongest_accelerometer_measurement_each_sample_for_lower": True,
        "strongest_S_zero_measurement_each_sample_for_lower": True,
        "nuisance_states_conditioned_known_for_translation_lower": True,
        "recursive_natural_interval_Riccati_subtraction_used": False,
        "deterministic_common_Loewner_lower_propagated": True,
        "word_horizon_s": WORD_HORIZON_S,
        "prediction_steps": n,
        "dt_s": h,
        "global_x_h_over_tau": global_x.as_list(),
        "sigma_process_lower_mps2": sigma_lo,
        "R_aw_lower": R_aw,
        "R_S_z_lower": R_S_z,
        "translation_covariance_upper": upper,
        "translation_covariance_upper_z": upper_z,
        "endpoint_common_lower_diagonal_z": [P[i][i].lo for i in range(4)],
        "relative_word_injection_floor_lower": rho,
        "useful_gate": USEFUL_GATE,
        "useful_margin_pass": gate_pass,
        "numerical_profile": {
            "max_prediction_cell_shave": max_prediction_eps,
            "max_measurement_shave": max_measurement_eps,
            "relative_prediction_shaves": relative_prediction_shaves,
            "relative_measurement_shaves": relative_measurement_shaves,
            "minimum_cross_cell_common_scale": min_common_alpha,
            "minimum_cross_cell_common_scale_step": limiting_step,
            "minimum_cross_cell_common_scale_cell": limiting_cell,
        },
        "theorem_identity": {
            "finite_word_concavity": "R_W(P)-D R_W(P)[P] >= R_W(0)",
            "selected_process_lower": "R_W(0) >= L_translation",
            "comparison": "L_translation >= delta_translation * Pbar_translation",
        },
        "pass": gate_pass,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "instantaneous_source_cell_partition_consumed",
        "global_source_interval_reselected_independently_each_step",
        "cell_index_forgotten_after_each_step",
        "full_4x4_translation_matrix_retained",
        "strongest_accelerometer_measurement_each_sample_for_lower",
        "strongest_S_zero_measurement_each_sample_for_lower",
        "nuisance_states_conditioned_known_for_translation_lower",
        "deterministic_common_Loewner_lower_propagated",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "old_P2_history_frontier_consumed",
        "tau_path_correlation_assumed", "sigma_path_correlation_assumed",
        "R_S_path_correlation_assumed", "determinant_or_e3_scalarization_used",
        "recursive_natural_interval_Riccati_subtraction_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("history_depth") != 0:
        f.append("instantaneous arithmetic cells acquired source history")
    if not isinstance(d.get("instantaneous_x_cells"), int) or d.get("instantaneous_x_cells", 0) <= 0:
        f.append("invalid instantaneous x cover")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("translation useful gate changed")
    rho = d.get("relative_word_injection_floor_lower")
    if not isinstance(rho, (int, float)) or not math.isfinite(float(rho)) or float(rho) < 0.0:
        f.append("translation finite-word floor is not finite nonnegative")
    expected = isinstance(rho, (int, float)) and float(rho) >= USEFUL_GATE
    if d.get("useful_margin_pass") is not expected or d.get("pass") is not expected:
        f.append("translation pass flag disagrees with finite-word margin")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "x_cells": d["instantaneous_x_cells"],
        "global_x_h_over_tau": d["global_x_h_over_tau"],
        "steps": d["prediction_steps"],
        "translation_delta": d["relative_word_injection_floor_lower"],
        "useful_gate": d["useful_gate"],
        "pass": d["pass"],
        "endpoint_diag": d["endpoint_common_lower_diagonal_z"],
        "profile": d["numerical_profile"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
