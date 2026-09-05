#!/usr/bin/env python3
"""Source-history-free finite-word translation lower for canonical OU-III P3.

The post-#489 one-sample P3 comparison is structurally microscopic because the
weakest integrated-OU direction carries powers of the 5 ms sample interval.
P3, however, is a recurrent finite-window theorem.  This producer therefore
keeps the complete one-axis [v,p,S,a_w] matrix over a full recurrent word.

The proof is deliberately stronger than a frozen-source scan.  At *every*
prediction step tau may be any value in the complete SEA3 dynamic invariant.
The same global interval x=h/tau is reevaluated independently at every step, so
no equality or path correlation between successive source values is assumed.
Sigma is replaced by its global positive lower endpoint and the strongest
possible accelerometer/S=0 measurements are applied at every sample.  Both
choices can only decrease the selected-process covariance and therefore give a
Loewner lower for every shipping word in the declared Normal-Live language.

Natural interval propagation of the subtractive Riccati formula is avoided.
After each prediction/measurement an interval matrix A is collapsed to a
single deterministic common Loewner lower L.  If C is the binary64 midpoint
and |A-C| <= E entrywise, a diagonally scaled row-sum bound gives

    A_real >= C - eps D^2

for every represented symmetric matrix.  This is the same validated
midpoint/radius idea that made the retained correlated-segment diagnostic
stable, but it is implemented here without importing or consuming any P2
source graph/history object.

The propagated coordinates are fixed physical coordinates

    z = [v/h, p/h^2, S/h^3, a_w],

so source motion never silently rescales the metric.  The endpoint lower is
compared directly with the canonical SEA3 covariance upper.  The resulting
relative floor is the translation contribution to the recurrent-word Riccati
injection; it is not a per-sample contraction claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p3_scaled_process as SCALED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_SOURCE_HISTORY_FREE_FINITE_WORD_TRANSLATION_LOWER"
WORD_HORIZON_S = 1.0
SERIES_ORDER = 14
USEFUL_GATE = 1.0e-18


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _zero(n: int = 4):
    return [[I(0.0) for _ in range(n)] for _ in range(n)]


def _transpose(A):
    return [list(row) for row in zip(*A)]


def _mul(A, B):
    if not A or not B or len(A[0]) != len(B):
        raise ValueError("matrix shape mismatch")
    out = [[I(0.0) for _ in range(len(B[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            s = I(0.0)
            for k in range(len(B)):
                s = s + A[i][k] * B[k][j]
            out[i][j] = s
    return out


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _scale(A, s: float):
    q = I(float(s))
    return [[q * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _sym(A):
    n = len(A)
    out = [[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
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
    """Enclose sum (-1)^n x^n/(n+m)! for m=1,2,3."""
    if m not in (1, 2, 3) or not (0.0 < x.lo <= x.hi < 0.05):
        raise ValueError("transition series outside audited small-x range")
    y = I(0.0)
    for n in range(SERIES_ORDER + 1):
        c = Fraction((-1) ** n, math.factorial(n + m))
        y = y + I(float(c)) * _ipow(x, n)
    first = x.hi ** (SERIES_ORDER + 1) / math.factorial(SERIES_ORDER + 1 + m)
    ratio = x.hi / float(SERIES_ORDER + m + 2)
    if not ratio < 1.0:
        raise RuntimeError("transition-series remainder ratio is not contractive")
    tail = up(first / (1.0 - ratio))
    return Interval(math.nextafter(y.lo - tail, -math.inf),
                    math.nextafter(y.hi + tail, math.inf))


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


def _scaled_process(x: Interval):
    """Shipping Q in z coordinates, before the common sigma^2 factor."""
    if x.hi < SCALED.BRANCH_X:
        B = SCALED.small_normalized_matrix(x)
        return _sym([[x * B[i][j] for j in range(4)] for i in range(4)])
    if x.lo >= SCALED.BRANCH_X and x.hi <= SCALED.NEAR_EXACT_SERIES_MAX_X:
        B = SCALED.near_exact_normalized_matrix(x)
        return _sym([[x * B[i][j] for j in range(4)] for i in range(4)])
    if x.lo < SCALED.BRANCH_X <= x.hi:
        left = Interval(x.lo, math.nextafter(SCALED.BRANCH_X, -math.inf))
        right = Interval(SCALED.BRANCH_X, x.hi)
        families = [_scaled_process(left), _scaled_process(right)]
        return [[Interval(min(A[i][j].lo for A in families),
                          max(A[i][j].hi for A in families))
                 for j in range(4)] for i in range(4)]
    raise ValueError("SEA3 x interval outside audited scaled-process range")


def _center_radius(A):
    A = _sym(A)
    n = len(A)
    center = [[0.0] * n for _ in range(n)]
    radius = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            a = A[i][j]
            c = float(a.lo + 0.5 * (a.hi - a.lo))
            c = min(max(c, a.lo), a.hi)
            center[i][j] = c
            radius[i][j] = up(max(c - a.lo, a.hi - c))
    for i in range(n):
        for j in range(i + 1, n):
            c = 0.5 * (center[i][j] + center[j][i])
            center[i][j] = center[j][i] = c
            r = up(max(radius[i][j], radius[j][i]))
            radius[i][j] = radius[j][i] = r
    return center, radius


def _common_point_lower(A):
    """One deterministic L with A_real >= L for every symmetric A_real in A."""
    center, radius = _center_radius(A)
    n = len(center)
    scales = []
    relative = True
    for i in range(n):
        if center[i][i] <= 0.0:
            relative = False
            break
        d = down(math.sqrt(center[i][i]))
        if not d > 0.0:
            relative = False
            break
        scales.append(d)

    if relative:
        eps = 0.0
        for i in range(n):
            row = 0.0
            for j in range(n):
                denom = down(scales[i] * scales[j])
                if not denom > 0.0:
                    relative = False
                    break
                row = up(row + up(radius[i][j] / denom))
            if not relative:
                break
            eps = max(eps, row)
        eps = up(eps)

    if relative:
        shifts = [up(eps * up(scales[i] * scales[i])) for i in range(n)]
        route = "RELATIVE_DIAGONAL"
    else:
        eps = 0.0
        for i in range(n):
            row = 0.0
            for j in range(n):
                row = up(row + radius[i][j])
            eps = max(eps, row)
        eps = up(eps)
        shifts = [eps] * n
        route = "ABSOLUTE_FALLBACK"

    L = [[I(center[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        L[i][i] = I(down(center[i][i] - shifts[i]))
    return _sym(L), eps, route


def _measurement_lower(L, coordinate: int, R_lower: float):
    """Scalar covariance update followed by a non-recursive common lower shave."""
    c = int(coordinate)
    den = L[c][c] + I(float(R_lower))
    if den.lo <= 0.0:
        raise RuntimeError("measurement lower innovation denominator lost positivity")
    col = [L[i][c] for i in range(4)]
    A = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            A[i][j] = L[i][j] - (col[i] * col[j]) / den
    return _common_point_lower(_sym(A))


def _minus_weight(P, rho: float, upper_z: list[float]):
    A = [[P[i][j] for j in range(4)] for i in range(4)]
    for i in range(4):
        A[i][i] = A[i][i] - I(float(rho)) * I(float(upper_z[i]))
    return _sym(A)


def _generalized_rho(P, upper_z: list[float]) -> float:
    if not symmetric_positive_definite_ldlt(P)[0]:
        return 0.0
    hi = min(P[i][i].lo / float(upper_z[i]) for i in range(4))
    if not (math.isfinite(hi) and hi > 0.0):
        return 0.0
    lo = 0.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        if symmetric_positive_definite_ldlt(_minus_weight(P, mid, upper_z))[0]:
            lo = mid
        else:
            hi = mid
    return down(lo)


def _load_tube(path: Path, tube_path: Path | None):
    tube = TUBE.build(path) if tube_path is None else json.loads(Path(tube_path).read_text(encoding="utf-8"))
    vf = TUBE.validate(tube)
    if vf:
        raise RuntimeError(f"SEA3 covariance tube validation failed: {vf}")
    return tube


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("finite-word translation proof may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"SEA3 dynamic source prerequisite failed: {df}")
    tube = _load_tube(path, tube_path)

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
    x = Interval.outward_bounds(down(h / tau_hi), up(h / tau_lo))
    if not (0.0 < x.lo <= x.hi < 0.05):
        raise RuntimeError("global SEA3 h/tau interval left audited transition range")

    F = _transition(x)
    Ft = _transpose(F)
    Q = _scale(_scaled_process(x), down(sigma_lo * sigma_lo))

    vc = tube["measurement_source"]
    R_aw = down(float(vc["acc_measurement_std_mps2"]) ** 2)
    axis_factor = min(map(float, tube["source_bounds"]["R_S_axis_std_factors"]))
    rs_lo = float(inv["R_S_applied"][0])
    rS = down(rs_lo * axis_factor)
    R_S_z = down(down(rS * rS) / up(h ** 6))
    if not (R_aw > 0.0 and R_S_z > 0.0):
        raise RuntimeError("strongest translation measurement variance lost positivity")

    P = _zero()
    max_prediction_eps = 0.0
    max_measurement_eps = 0.0
    relative_prediction_shaves = 0
    relative_measurement_shaves = 0
    for _step in range(1, n + 1):
        pred_interval = _sym(_add(_mul(_mul(F, P), Ft), Q))
        P, eps, route = _common_point_lower(pred_interval)
        max_prediction_eps = max(max_prediction_eps, eps)
        relative_prediction_shaves += int(route == "RELATIVE_DIAGONAL")

        P, eps, route = _measurement_lower(P, 3, R_aw)
        max_measurement_eps = max(max_measurement_eps, eps)
        relative_measurement_shaves += int(route == "RELATIVE_DIAGONAL")
        # Hypothetical S=0 every sample is stronger than every shipping pseudo
        # schedule and therefore produces a valid lower covariance.
        P, eps, route = _measurement_lower(P, 2, R_S_z)
        max_measurement_eps = max(max_measurement_eps, eps)
        relative_measurement_shaves += int(route == "RELATIVE_DIAGONAL")

    upper = list(map(float, tube["covariance_upper"]["translation_diagonal_variance_upper"]))
    dscale = [h, h * h, h * h * h, 1.0]
    upper_z = [up(upper[i] / down(dscale[i] * dscale[i])) for i in range(4)]
    rho = _generalized_rho(P, upper_z)
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
        "current_source_cell_partition_consumed": False,
        "global_source_interval_reselected_independently_each_step": True,
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
        "global_x_h_over_tau": x.as_list(),
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
            "max_prediction_relative_or_absolute_shave": max_prediction_eps,
            "max_measurement_relative_or_absolute_shave": max_measurement_eps,
            "relative_prediction_shaves": relative_prediction_shaves,
            "relative_measurement_shaves": relative_measurement_shaves,
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
        "global_source_interval_reselected_independently_each_step",
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
        "current_source_cell_partition_consumed", "tau_path_correlation_assumed",
        "sigma_path_correlation_assumed", "R_S_path_correlation_assumed",
        "determinant_or_e3_scalarization_used",
        "recursive_natural_interval_Riccati_subtraction_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("translation useful gate changed")
    rho = d.get("relative_word_injection_floor_lower")
    if not isinstance(rho, (int, float)) or not math.isfinite(float(rho)) or float(rho) < 0.0:
        f.append("translation finite-word floor is not finite nonnegative")
    expected = isinstance(rho, (int, float)) and float(rho) >= USEFUL_GATE
    if d.get("useful_margin_pass") is not expected or d.get("pass") is not expected:
        f.append("translation pass flag disagrees with finite-word margin")
    if abs(float(d.get("word_horizon_s", 0.0)) - WORD_HORIZON_S) > 1e-12:
        f.append("translation word horizon changed")
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
        "global_x_h_over_tau": d["global_x_h_over_tau"],
        "steps": d["prediction_steps"],
        "translation_delta": d["relative_word_injection_floor_lower"],
        "useful_gate": d["useful_gate"],
        "pass": d["pass"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
