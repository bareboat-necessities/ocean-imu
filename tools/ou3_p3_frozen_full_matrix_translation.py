#!/usr/bin/env python3
"""Rigorous frozen-source full-matrix translation P3 headroom certificate.

This producer answers one narrow but decisive question before the source-varying
P3 fixed point is built: if one of the retained 800 physical tuner cells is
held fixed, does the *matrix* integrated-OU process covariance remain useful
after the strongest translation measurements that can occur over a 1.25 s
suffix?

The earlier determinant/e3 route reduces the complete [v,p,S,a_w] covariance to
one scalar before measurement attenuation.  That loses the chain geometry.  We
instead propagate the exact shipping one-step normalized process covariance as a
4x4 interval matrix in

    y = [v/(sigma h), p/(sigma h^2), S/(sigma h^3), a_w/sigma].

For a frozen source cell, x=h/tau is constant.  The exact normalized transition
is composed for 250 deployed 5-ms steps.  After every prediction we apply both
of the strongest translation measurements that could possibly reduce the
selected process covariance:

* an accepted accelerometer a_w measurement every sample, with the smallest
  normalized R_a allowed by the source cell; and
* an S=0 pseudo measurement every sample, even though the shipping cadence is
  never faster than one per sample, with the strongest horizontal R_S axis.

Conditioning attitude / bias nuisance states as known and accepting packets that
may actually be rejected can only decrease this selected-process covariance.
Thus the result is a lower covariance for the frozen translation subsystem.  We
then test directly, with outward-rounded LDL^T,

    D_lo P_selected D_lo - 1e-18 diag(Sigma_upper) >> 0

on every x subdivision of every physical source cell.  The covariance upper is
the retained finite-memory cell upper; no replay values enter.

This is intentionally NOT source-complete P3: tau/sigma/R_S are frozen inside
the 1.25-s suffix.  Passing all 800 cells proves that the remaining obstacle is
source variation / invariant-metric composition, not static process excitation
or measurement attenuation.  It cannot promote P3/P4/P5.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p3_scaled_process as SCALED
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
HORIZON_S = 1.25
X_SUBCELLS = 16
SERIES_ORDER = 12


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _mat_zero(n: int = 4):
    z = I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _transpose(A):
    return [list(x) for x in zip(*A)]


def _mul(A, B):
    rows, inner, cols = len(A), len(B), len(B[0])
    if len(A[0]) != inner:
        raise ValueError("matrix shape mismatch")
    out = _mat_zero(rows)
    if cols != rows:
        out = [[I(0.0) for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            s = I(0.0)
            for k in range(inner):
                s = s + A[i][k] * B[k][j]
            out[i][j] = s
    return out


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _sym(A):
    out = [[A[i][j] for j in range(len(A))] for i in range(len(A))]
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            lo = min(out[i][j].lo, out[j][i].lo)
            hi = max(out[i][j].hi, out[j][i].hi)
            x = Interval(lo, hi)
            out[i][j] = x
            out[j][i] = x
    return out


def _ipow(x: Interval, n: int) -> Interval:
    y = I(1.0)
    for _ in range(int(n)):
        y = y * x
    return y


def _normalized_integral_series(x: Interval, m: int) -> Interval:
    """Enclose sum_{n>=0} (-1)^n x^n/(n+m)! for m=1,2,3."""
    if m not in (1, 2, 3) or not (0.0 < x.lo <= x.hi < 0.05):
        raise ValueError("transition series outside audited small-x range")
    y = I(0.0)
    for n in range(SERIES_ORDER + 1):
        c = Fraction((-1) ** n, math.factorial(n + m))
        cf = I(float(c))
        y = y + cf * _ipow(x, n)
    first = x.hi ** (SERIES_ORDER + 1) / math.factorial(SERIES_ORDER + 1 + m)
    ratio = x.hi / float(SERIES_ORDER + m + 2)
    if not ratio < 1.0:
        raise RuntimeError("transition-series tail ratio is not contractive")
    tail = math.nextafter(first / (1.0 - ratio), math.inf)
    return Interval(math.nextafter(y.lo - tail, -math.inf), math.nextafter(y.hi + tail, math.inf))


def _transition(x: Interval):
    k1 = _normalized_integral_series(x, 1)
    k2 = _normalized_integral_series(x, 2)
    k3 = _normalized_integral_series(x, 3)
    e = VT.exp_interval(-x)
    return [
        [I(1.0), I(0.0), I(0.0), k1],
        [I(1.0), I(1.0), I(0.0), k2],
        [I(0.5), I(1.0), I(1.0), k3],
        [I(0.0), I(0.0), I(0.0), e],
    ]


def _scaled_Q(x: Interval):
    if x.hi < SCALED.BRANCH_X:
        B = SCALED.small_normalized_matrix(x)
    elif x.lo >= SCALED.BRANCH_X and x.hi <= SCALED.NEAR_EXACT_SERIES_MAX_X:
        B = SCALED.near_exact_normalized_matrix(x)
    else:
        raise ValueError("x cell crosses scaled-process branch; split before evaluation")
    return _sym([[x * B[i][j] for j in range(4)] for i in range(4)])


def _measurement_update(P, coordinate: int, R: float):
    den = P[coordinate][coordinate] + I(R)
    if den.lo <= 0.0:
        raise RuntimeError("measurement innovation denominator lost positivity")
    col = [P[i][coordinate] for i in range(4)]
    out = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            out[i][j] = P[i][j] - (col[i] * col[j]) / den
    return _sym(out)


def _split_x(lo: float, hi: float, count: int):
    cuts = {float(lo), float(hi)}
    if lo < SCALED.BRANCH_X < hi:
        cuts.add(SCALED.BRANCH_X)
    base = sorted(cuts)
    out = []
    total_log = math.log(hi / lo)
    for a, b in zip(base, base[1:]):
        n = max(1, int(math.ceil(count * math.log(b / a) / total_log))) if total_log > 0 else 1
        ratio = (b / a) ** (1.0 / n)
        edges = [a]
        for _ in range(n - 1):
            edges.append(edges[-1] * ratio)
        edges.append(b)
        for k in range(n):
            aa = math.nextafter(edges[k], -math.inf) if k else a
            bb = math.nextafter(edges[k + 1], math.inf) if k + 1 < n else b
            if b == SCALED.BRANCH_X and k + 1 == n:
                bb = math.nextafter(SCALED.BRANCH_X, -math.inf)
            if a == SCALED.BRANCH_X and k == 0:
                aa = SCALED.BRANCH_X
            out.append(Interval(aa, bb))
    return out


def _gate_matrix(Pnorm, sigma_lower: float, h: float, upper: list[float], gate: float):
    powers = [h, h * h, h * h * h, 1.0]
    d = [sigma_lower * z for z in powers]
    A = [[I(d[i]) * Pnorm[i][j] * I(d[j]) for j in range(4)] for i in range(4)]
    for i in range(4):
        A[i][i] = A[i][i] - I(gate) * I(float(upper[i]))
    return _sym(A)


def _one_subcell(x: Interval, steps: int, R_aw_norm: float, R_S_norm: float,
                 sigma_lower: float, h: float, upper: list[float], gate: float) -> dict:
    F = _transition(x)
    Ft = _transpose(F)
    Q = _scaled_Q(x)
    P = _mat_zero()
    first_spd_loss = None
    for step in range(1, steps + 1):
        P = _sym(_add(_mul(_mul(F, P), Ft), Q))
        P = _measurement_update(P, 3, R_aw_norm)
        P = _measurement_update(P, 2, R_S_norm)
        if step in (1, steps) and not symmetric_positive_definite_ldlt(P)[0]:
            first_spd_loss = step
            break
    gate_ok = False
    if first_spd_loss is None:
        gate_ok = bool(symmetric_positive_definite_ldlt(
            _gate_matrix(P, sigma_lower, h, upper, gate)
        )[0])
    return {
        "x_h_over_tau": x.as_list(),
        "posterior_interval_spd": first_spd_loss is None,
        "first_interval_spd_loss_step": first_spd_loss,
        "useful_gate_spd": gate_ok,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("frozen full-matrix diagnostic must not be trajectory fitted")
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("frozen translation certificate requires zero lever arm")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        raise RuntimeError("frozen translation certificate requires transparent vibration guard")

    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"source nodes invalid: {nf}")
    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    steps = int(math.floor(HORIZON_S / h + 1e-12))
    if steps * h > HORIZON_S + 1e-12 or steps < 1:
        raise RuntimeError("invalid deployed horizon discretization")
    gate = float(BASE.MIN_USEFUL_DELTA)
    if gate != 1.0e-18:
        raise RuntimeError("frozen matrix certificate requires unchanged 1e-18 useful gate")
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    vc = BASE.VECTOR.build()["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    axis_factor = min(BASE.source_rs_axis_std_factors())

    rows = []
    gate_pass = 0
    interval_closed = 0
    worst_fail = None
    for node in nodes_payload["nodes"]:
        tau = Interval(*map(float, node["tau_s"]))
        sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
        rs = Interval(*map(float, node["R_S_filter_std"]))
        upper, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)
        upper = [float(x) for x in upper]
        if float(timing["word_horizon_s_lower"]) + 1e-15 < HORIZON_S:
            raise RuntimeError("chosen frozen suffix exceeds certified P3 word")
        R_aw_norm = BASE.down(ra / BASE.up(sigma.hi * sigma.hi))
        rS = BASE.down(rs.lo * axis_factor)
        denomS = BASE.up(sigma.hi * sigma.hi * h ** 6)
        R_S_norm = BASE.down(BASE.down(rS * rS) / denomS)
        if not (R_aw_norm > 0.0 and R_S_norm > 0.0):
            raise RuntimeError("normalized measurement noise lost positivity")
        xlo = BASE.down(h / tau.hi)
        xhi = BASE.up(h / tau.lo)
        subrows = []
        for xcell in _split_x(xlo, xhi, X_SUBCELLS):
            subrows.append(_one_subcell(
                xcell, steps, R_aw_norm, R_S_norm,
                sigma.lo, h, upper, gate,
            ))
        closed = all(r["posterior_interval_spd"] for r in subrows)
        passed = closed and all(r["useful_gate_spd"] for r in subrows)
        interval_closed += int(closed)
        gate_pass += int(passed)
        row = {
            "index": int(node["index"]),
            "tau_index": int(node["tau_index"]),
            "sigma_raw_index": int(node["sigma_raw_index"]),
            "R_S_index": int(node["R_S_index"]),
            "tau_s": tau.as_list(),
            "sigma_filter_committed_mps2": sigma.as_list(),
            "R_S_filter_std": rs.as_list(),
            "Sigma_translation_diagonal_upper": upper,
            "R_aw_normalized_lower": R_aw_norm,
            "R_S_normalized_lower": R_S_norm,
            "x_subcells": len(subrows),
            "interval_propagation_closed": closed,
            "useful_gate_pass": passed,
            "subcells": subrows,
        }
        if not passed and worst_fail is None:
            worst_fail = row
        rows.append(row)

    all_pass = gate_pass == NODES.EXPECTED_STATES
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FROZEN_SOURCE_FULL_MATRIX_TRANSLATION_SELECTED_PROCESS_CERTIFICATE",
        "source_partition_only": True,
        "frozen_source_inside_suffix": True,
        "source_complete_time_variation_covered": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "full_4x4_translation_matrix_retained": True,
        "determinant_e3_scalarization_used": False,
        "accepted_accelerometer_every_sample_for_lower_bound": True,
        "S_zero_every_sample_for_lower_bound": True,
        "nuisance_states_conditioned_known": True,
        "measurement_model_stronger_than_every_shipping_translation_branch": True,
        "validated_interval_arithmetic": True,
        "ordinary_floating_eigensolver_used": False,
        "useful_gate": gate,
        "horizon_s": HORIZON_S,
        "prediction_steps": steps,
        "x_subdivision_target": X_SUBCELLS,
        "source_cells_scanned": len(rows),
        "interval_propagation_closed_cells": interval_closed,
        "useful_gate_pass_cells": gate_pass,
        "all_frozen_source_cells_useful": all_pass,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "first_nonpassing_cell": worst_fail,
        "rows": rows,
        "next_obligation": (
            "if all frozen cells pass, replace the frozen x/sigma/R_S cell by the retained staged/committed source automaton and propagate a source-indexed matrix lower/invariant across every finite-speed segment; do not promote this frozen diagnostic as source-complete P3"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FROZEN_SOURCE_FULL_MATRIX_TRANSLATION_SELECTED_PROCESS_CERTIFICATE":
        f.append("wrong qualification")
    for key in (
        "source_partition_only", "frozen_source_inside_suffix", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch", "full_4x4_translation_matrix_retained",
        "accepted_accelerometer_every_sample_for_lower_bound", "S_zero_every_sample_for_lower_bound",
        "nuisance_states_conditioned_known", "measurement_model_stronger_than_every_shipping_translation_branch",
        "validated_interval_arithmetic",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_complete_time_variation_covered", "trajectory_replay_used", "filter_changed",
        "declared_domain_changed", "determinant_e3_scalarization_used",
        "ordinary_floating_eigensolver_used", "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", 0.0)) != 1.0e-18:
        f.append("useful gate changed")
    if int(d.get("source_cells_scanned", 0)) != NODES.EXPECTED_STATES:
        f.append("did not scan all 800 source cells")
    if not 0 <= int(d.get("interval_propagation_closed_cells", -1)) <= NODES.EXPECTED_STATES:
        f.append("invalid interval-closed count")
    if not 0 <= int(d.get("useful_gate_pass_cells", -1)) <= NODES.EXPECTED_STATES:
        f.append("invalid useful-cell count")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "horizon_s": d["horizon_s"],
        "prediction_steps": d["prediction_steps"],
        "source_cells_scanned": d["source_cells_scanned"],
        "interval_propagation_closed_cells": d["interval_propagation_closed_cells"],
        "useful_gate_pass_cells": d["useful_gate_pass_cells"],
        "all_frozen_source_cells_useful": d["all_frozen_source_cells_useful"],
        "first_nonpassing_cell": None if d["first_nonpassing_cell"] is None else {
            "index": d["first_nonpassing_cell"]["index"],
            "tau_index": d["first_nonpassing_cell"]["tau_index"],
            "sigma_raw_index": d["first_nonpassing_cell"]["sigma_raw_index"],
            "R_S_index": d["first_nonpassing_cell"]["R_S_index"],
            "closed": d["first_nonpassing_cell"]["interval_propagation_closed"],
            "useful": d["first_nonpassing_cell"]["useful_gate_pass"],
        },
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
