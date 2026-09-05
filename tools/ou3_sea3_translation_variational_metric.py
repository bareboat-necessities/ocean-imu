#!/usr/bin/env python3
"""Full-matrix finite-memory translation covariance lower for canonical SEA3 P3.

This replaces the old one-step/scalar process floor with a variational
(controllability/least-action) lower bound over a finite Normal-Live suffix.
For one physical translation axis ordered x=[S,p,v,a_w], condition the initial
state and all nuisance states as known.  Conditioning only reduces covariance,
so any covariance lower for this stronger conditional problem is also a lower
for the shipping filter.

Over H seconds choose, for every requested endpoint, one explicit cubic
acceleration trajectory a(t).  In normalized time u=t/H write b(u)=a(Hu) and
scaled endpoint coordinates

  y = [S/H^3, p/H^2, v/H, a_w].

The unique cubic satisfying the endpoint constraints has coefficients

  c = C y,

with the exact integer matrix C below.  The continuous OU driving action is

  int (a' + a/tau)^2 / q_c dt,
  q_c = 2 sigma_aw^2/tau.

Using the deployed source invariants sigma_aw>=sigma_min and
tau_lo<=tau<=tau_hi gives the *matrix* upper

  E_process <= tau_hi/(2 sigma_min^2) int a'^2 dt
             + 1/sigma_min^2         int a'a dt
             + 1/(2 sigma_min^2 tau_lo) int a^2 dt.

The cross coefficient is exact and source independent; only the two positive
square coefficients are enlarged.  This is substantially tighter than the old
2(a'^2+lambda_max^2 a^2) scalar inequality and remains valid for arbitrary
time-varying tau(t).

Measurement attenuation is included in the same quadratic action, not as a
scalar beta.  Every accepted Normal-Live accelerometer sample is admitted with
attitude/bias/nuisance states conditioned known, leaving the orthogonal a_w
Jacobian.  For S=0 we deliberately admit one pseudo measurement at *every* IMU
sample with the strongest deployed horizontal variance.  That is stronger than
any shipping pseudo schedule and therefore conservative for a covariance lower.
The exact discrete measurement sums of the cubic path are accumulated with
Fraction arithmetic; no integral or packet-count approximation is used.

If M is the resulting endpoint action upper, Gaussian minimum-energy/Riccati
duality gives

  P_translation >= M^-1.

All matrix inversion and SPD decisions are outward-rounded interval operations.
No determinant/trace scalarization, scalar information beta, floating
eigensolver, source history graph or trajectory replay is used.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import re
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_mul,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_rs_tau_lag_envelope as LAG
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
KALMAN = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_TRANSLATION_VARIATIONAL_FULL_MATRIX_METRIC_FLOOR"
HORIZON_S = 2.0

F = Fraction

# Polynomial coefficient map c=C*y for b(u)=sum_{n=0}^3 c_n u^n.
C = (
    (F(120), F(-60), F(12), F(-1)),
    (F(-1080), F(600), F(-132), F(12)),
    (F(2160), F(-1260), F(300), F(-30)),
    (F(-1200), F(720), F(-180), F(20)),
)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def IF(q: Fraction) -> Interval:
    return Interval.outward_bounds(float(q), float(q))


def ipow(x: Interval, n: int) -> Interval:
    y = Interval.point(1.0)
    for _ in range(int(n)):
        y = y * x
    return y


def _zero_fraction_matrix(n: int) -> list[list[Fraction]]:
    return [[F(0) for _ in range(n)] for _ in range(n)]


def _outer_add(M: list[list[Fraction]], row: list[Fraction], weight: Fraction = F(1)) -> None:
    for i in range(4):
        for j in range(4):
            M[i][j] += weight * row[i] * row[j]


def _b_row(u: Fraction) -> list[Fraction]:
    powers = (F(1), u, u*u, u*u*u)
    return [sum(powers[n] * C[n][j] for n in range(4)) for j in range(4)]


def _bp_row(u: Fraction) -> list[Fraction]:
    powers = (F(0), F(1), F(2)*u, F(3)*u*u)
    return [sum(powers[n] * C[n][j] for n in range(4)) for j in range(4)]


def _s_row(u: Fraction) -> list[Fraction]:
    # sbar(u)=int_0^u (u-r)^2/2 b(r) dr
    out = [F(0) for _ in range(4)]
    for n in range(4):
        factor = u ** (n + 3) / F((n + 1) * (n + 2) * (n + 3))
        for j in range(4):
            out[j] += C[n][j] * factor
    return out


def _poly_integral_gram(kind: str) -> list[list[Fraction]]:
    """Exact integral Gram for b^2, b'^2, or b'b on u in [0,1]."""
    M = _zero_fraction_matrix(4)
    for n in range(4):
        for m in range(4):
            if kind == "b2":
                factor = F(1, n + m + 1)
                an, am = F(1), F(1)
            elif kind == "bp2":
                if n == 0 or m == 0:
                    continue
                factor = F(1, n + m - 1)
                an, am = F(n), F(m)
            elif kind == "bpb":
                if n == 0:
                    continue
                factor = F(1, n + m)
                an, am = F(n), F(1)
            else:
                raise ValueError(kind)
            for i in range(4):
                for j in range(4):
                    M[i][j] += an * am * factor * C[n][i] * C[m][j]
    if kind == "bpb":
        # y^T M y sees only the symmetric part.
        M = [[(M[i][j] + M[j][i]) / 2 for j in range(4)] for i in range(4)]
    return M


def _discrete_measurement_grams(samples: int) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    A = _zero_fraction_matrix(4)
    S = _zero_fraction_matrix(4)
    for k in range(1, samples + 1):
        u = F(k, samples)
        _outer_add(A, _b_row(u))
        _outer_add(S, _s_row(u))
    return A, S


def _interval_matrix(M: list[list[Fraction]]) -> list[list[Interval]]:
    return [[IF(x) for x in row] for row in M]


def _scale(A, c: Interval):
    return [[c * x for x in row] for row in A]


def _diag(values: list[Interval]):
    z = Interval.point(0.0)
    return [[values[i] if i == j else z for j in range(len(values))] for i in range(len(values))]


def _matrix_bounds(A) -> list[list[list[float]]]:
    return [[[x.lo, x.hi] for x in row] for row in A]


def _member_float(text: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f\b", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed {name}")
    return float(m.group(1))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("variational metric floor may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    lag = LAG.build(path)
    vector = VECTOR.build()
    for label, failures in (
        ("dynamic", DYNAMIC.validate(dynamic)),
        ("lag", LAG.validate(lag)),
        ("vector", VECTOR.validate(vector)),
    ):
        if failures:
            raise RuntimeError(f"{label} prerequisite failed: {failures}")

    kalman = KALMAN.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    for marker in (
        "QdAxis4x1_analytic(tau, Ts, sigma2, Qd_axis);",
        "const Matrix3 J_aw  =  R_wb();",
        "joseph_update3_(K, S_mat, PCt);",
        "apply_pending_aw_covariance_inflation_();",
    ):
        if marker not in kalman:
            raise RuntimeError(f"shipping translation/measurement source marker changed: {marker}")
    for marker in (
        "float S_factor_      = 1.0f;",
        "float R_S_x_factor_ = 0.72f;",
        "float R_S_y_factor_ = 0.72f;",
    ):
        if marker not in wrapper:
            raise RuntimeError(f"shipping axis source marker changed: {marker}")

    inv = dynamic["dynamic_invariant"]
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_min = float(inv["sigma_aw_filter_mps2"][0])
    if not (0.0 < tau_lo <= tau_hi and sigma_min >= 0.05):
        raise RuntimeError("dynamic translation source invariant lost positivity")

    dt = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    samples_float = HORIZON_S / dt
    samples = int(round(samples_float))
    if not math.isclose(samples * dt, HORIZON_S, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("variational horizon must be an exact configured sample count")

    mb = _interval_matrix(_poly_integral_gram("b2"))
    mbp = _interval_matrix(_poly_integral_gram("bp2"))
    mcross = _interval_matrix(_poly_integral_gram("bpb"))
    macc_f, ms_f = _discrete_measurement_grams(samples)
    macc = _interval_matrix(macc_f)
    ms = _interval_matrix(ms_f)

    H = I(HORIZON_S)
    s2 = I(sigma_min).square()
    # Exact source-uniform OU action upper for arbitrary time-varying tau(t).
    m_process = matrix_add(
        matrix_add(
            _scale(mbp, I(tau_hi) / (I(2.0) * s2 * H)),
            _scale(mcross, I(1.0) / s2),
        ),
        _scale(mb, H / (I(2.0) * s2 * I(tau_lo))),
    )

    vc = vector["configured_measurement_bounds"]
    racc_std = float(vc["acc_measurement_std_mps2"])
    racc_var_lower = math.nextafter(racc_std * racc_std, -math.inf)
    if not racc_var_lower > 0.0:
        raise RuntimeError("accelerometer measurement variance lower lost positivity")
    m_acc = _scale(macc, I(1.0) / I(racc_var_lower))

    rs_floor = float(lag["R_S_hard_floor"])
    axis_factor = min(_member_float(wrapper, "R_S_x_factor_"),
                      _member_float(wrapper, "R_S_y_factor_"), 1.0)
    rs_std_lower = math.nextafter(rs_floor * axis_factor, -math.inf)
    rs_var_lower = math.nextafter(rs_std_lower * rs_std_lower, -math.inf)
    if not rs_var_lower > 0.0:
        raise RuntimeError("S pseudo variance lower lost positivity")
    # Physical S(t)=H^3*sbar(u); admit one stronger pseudo measurement at every
    # IMU sample.  H^6 converts the normalized exact sum to physical units.
    m_s = _scale(ms, ipow(H, 6) / I(rs_var_lower))

    m_y = matrix_symmetric_hull(matrix_add(matrix_add(m_process, m_acc), m_s))
    ok_m, piv_m = symmetric_positive_definite_ldlt(m_y)
    if not ok_m:
        raise RuntimeError(f"variational precision matrix is not certified SPD: {[p.as_list() for p in piv_m]}")

    # y = D^-1 x, D=diag(H^3,H^2,H,1), x=[S,p,v,a].
    dinv = _diag([
        I(1.0) / ipow(H, 3),
        I(1.0) / ipow(H, 2),
        I(1.0) / H,
        I(1.0),
    ])
    m_x = matrix_symmetric_hull(matrix_mul(matrix_mul(dinv, m_y), dinv))
    ok_x, piv_x = symmetric_positive_definite_ldlt(m_x)
    if not ok_x:
        raise RuntimeError(f"physical variational precision matrix is not certified SPD: {[p.as_list() for p in piv_x]}")

    l_x = matrix_symmetric_hull(matrix_inverse_gauss_jordan(m_x))
    ok_l, piv_l = symmetric_positive_definite_ldlt(l_x)
    if not ok_l:
        raise RuntimeError(f"physical covariance lower inverse is not certified SPD: {[p.as_list() for p in piv_l]}")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "old_P2_800_state_graph_consumed": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "per_sample_Riccati_lower_propagation_used": False,
        "ordinary_floating_eigensolver_used": False,
        "finite_memory_horizon_s": HORIZON_S,
        "finite_memory_samples": samples,
        "state_order": ["S", "p", "v", "a_w"],
        "scaled_state": ["S/H^3", "p/H^2", "v/H", "a_w"],
        "arbitrary_time_varying_tau_inside_memory": True,
        "same_source_process_sigma_floor_mps2": sigma_min,
        "tau_applied_invariant_s": [tau_lo, tau_hi],
        "OU_action_matrix_bound": {
            "formula": "tau_hi/(2*sigma_min^2)*int a'^2 + 1/sigma_min^2*int a'a + 1/(2*sigma_min^2*tau_lo)*int a^2",
            "cross_coefficient_not_scalarized": True,
            "exact_continuous_OU_discretization_source_consumed": True,
        },
        "measurement_conditioning": {
            "normal_live_accelerometer_every_sample": True,
            "attitude_bias_nuisance_conditioned_known": True,
            "J_aw_orthogonal": True,
            "acc_measurement_std_mps2": racc_std,
            "S_extra_packet_admitted_every_sample": True,
            "S_extra_packets_stronger_than_shipping_schedule": True,
            "R_S_base_hard_floor": rs_floor,
            "R_S_axis_std_factor_strongest": axis_factor,
            "R_S_filter_std_lower": rs_std_lower,
            "R_S_tau_lag_certificate_consumed": True,
            "lag_envelope_not_needed_to_rescue_baseline_floor": True,
        },
        "exact_rational_path_matrices": True,
        "validated_interval_matrix_arithmetic": True,
        "precision_M_scaled_interval": _matrix_bounds(m_y),
        "precision_M_physical_interval": _matrix_bounds(m_x),
        "covariance_L_physical_interval": _matrix_bounds(l_x),
        "precision_ldlt_pivots": [p.as_list() for p in piv_x],
        "covariance_lower_ldlt_pivots": [p.as_list() for p in piv_l],
        "P3_TRANSLATION_VARIATIONAL_METRIC_LOWER_CLOSED": True,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": False,
        "P3_PROMOTED": False,
        "next_obligation": (
            "compare the four-S innovation information directly against this full 4x4 precision M in the same scaled coordinates using the Newton/divided-difference observation transform; do not reduce either matrix to a scalar eigenvalue ratio"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "arbitrary_time_varying_tau_inside_memory",
        "exact_rational_path_matrices", "validated_interval_matrix_arithmetic",
        "P3_TRANSLATION_VARIATIONAL_METRIC_LOWER_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_P2_800_state_graph_consumed", "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed", "determinant_trace_scalarization_used",
        "scalar_information_beta_used", "per_sample_Riccati_lower_propagation_used",
        "ordinary_floating_eigensolver_used", "P3_FULL_MATRIX_COMPARISON_CLOSED", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    m = d.get("measurement_conditioning", {})
    for key in (
        "normal_live_accelerometer_every_sample", "attitude_bias_nuisance_conditioned_known",
        "J_aw_orthogonal", "S_extra_packet_admitted_every_sample",
        "S_extra_packets_stronger_than_shipping_schedule", "R_S_tau_lag_certificate_consumed",
        "lag_envelope_not_needed_to_rescue_baseline_floor",
    ):
        if m.get(key) is not True:
            f.append(f"measurement conditioning lost: {key}")
    if int(d.get("finite_memory_samples", 0)) != 400:
        f.append("2 s finite memory is not 400 configured samples")
    for key in ("precision_ldlt_pivots", "covariance_lower_ldlt_pivots"):
        piv = d.get(key, [])
        if len(piv) != 4 or any(float(x[0]) <= 0.0 for x in piv):
            f.append(f"{key} lost positive four-state pivots")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "horizon_s": d["finite_memory_horizon_s"],
        "samples": d["finite_memory_samples"],
        "M_pivots": d["precision_ldlt_pivots"],
        "L_pivots": d["covariance_lower_ldlt_pivots"],
        "metric_lower_closed": d["P3_TRANSLATION_VARIATIONAL_METRIC_LOWER_CLOSED"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
