#!/usr/bin/env python3
"""Direct generalized-matrix P3 comparison for deployed OU-III.

This layer removes the last scalarization in the source-reachable P3 proof.
The conditioned backend already certifies the integrated-OU process family in
an exact rational C=R L^{-1} congruence basis.  Here the same congruence is
applied to both the word-noise lower matrix and covariance upper matrix and we
certify directly

    Omega_word - delta * Sigma_upper  >>  0

with outward-rounded interval LDL^T.

No theorem domain, filter gain, source schedule, or usefulness gate is changed.
The translation process matrix is retained as a 4x4 interval family for each
axis.  A scalar measurement-information norm is used only to obtain a uniform
posterior *matrix* lower bound; the process matrix is not collapsed to its
smallest eigenvalue.  The already-established directional covariance upper
bounds are retained as a diagonal Loewner upper matrix rather than collapsed to
one global maximum.  H and A are finally checked as assembled 18x18 and 21x21
matrices in the same conditioned coordinates.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_source_reachable_matrix_p3 as BASE
import ou3_source_reachable_matrix_p3_factored as FACTORED

DEFAULT_DOMAIN = FACTORED.DEFAULT_DOMAIN
SCHEMA = FACTORED.SCHEMA
MIN_USEFUL_DELTA = FACTORED.MIN_USEFUL_DELTA

_ORIGINAL_MODE_CELL = BASE.mode_cell


def _zero_matrix(n: int):
    z = Interval.point(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _diag_matrix(values: list[float]):
    A = _zero_matrix(len(values))
    for i, v in enumerate(values):
        A[i][i] = Interval.outward_bounds(v, v)
    return A


def _scale_matrix(A, s: float):
    q = Interval.outward_bounds(s, s)
    return [[q * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _subtract_delta(Omega, Sigma, delta: float):
    q = Interval.outward_bounds(delta, delta)
    return BASE.matrix_symmetric_hull([
        [Omega[i][j] - q * Sigma[i][j] for j in range(len(Omega))]
        for i in range(len(Omega))
    ])


def _spd_at_delta(Omega, Sigma, delta: float) -> bool:
    ok, _ = symmetric_positive_definite_ldlt(_subtract_delta(Omega, Sigma, delta))
    return ok


def _certified_generalized_delta(Omega, Sigma, gate: float) -> float:
    """Largest located certified delta lower bound, using log-domain bracketing."""
    gate = float(gate)
    if gate <= 0.0:
        raise ValueError("positive usefulness gate required")

    if _spd_at_delta(Omega, Sigma, gate):
        lo = gate
        hi = gate
        while hi < 1.0:
            trial = min(1.0, hi * 10.0)
            if not _spd_at_delta(Omega, Sigma, trial):
                hi = trial
                break
            lo = trial
            hi = trial
            if trial == 1.0:
                return BASE.down(lo)
        if hi == lo:
            return BASE.down(lo)
    else:
        hi = gate
        lo = gate
        for _ in range(80):
            lo /= 10.0
            if _spd_at_delta(Omega, Sigma, lo):
                break
        else:
            return 0.0

    for _ in range(44):
        mid = math.sqrt(lo * hi)
        if _spd_at_delta(Omega, Sigma, mid):
            lo = mid
        else:
            hi = mid
    return BASE.down(lo)


def _matrix_abs_row_sum_upper(A) -> float:
    best = 0.0
    for row in A:
        total = 0.0
        for a in row:
            total = BASE.up(total + a.abs_upper())
        best = max(best, total)
    return BASE.up(best)


def _measurement_beta_upper(mode: str, sigma: Interval, rs: Interval,
                            live: dict, vector: dict, process: dict,
                            sched: dict) -> float:
    h = sched["dt_s"]
    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qba_d = BASE.pos(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "active bias discrete process",
    )
    sS2 = (sigma.lo * h * h * h) ** 2
    sa2 = sigma.lo * sigma.lo
    betaS = BASE.up(sS2 / (rs.lo * rs.lo))
    betaAcc = BASE.up(
        (fhi * fhi * qtheta + sa2 + (qba_d if mode == "A" else 0.0)) / ra
    )
    betaMag = BASE.up((mhi * mhi * qtheta) / rm)
    return BASE.up(BASE.up(betaS + betaAcc) + betaMag)


def _translation_direct_blocks(x: Interval, sigma: Interval, raw: dict,
                               beta: float, sched: dict):
    """Return conditioned (Omega,Sigma), process factor and norm diagnostics.

    The process family is represented as Q_scaled = x * Shape(x).  We keep the
    already-conditioned Shape(x)=Q_scaled/x interval family, transform it by C,
    and use x.lo only after the congruence.  For every concrete source value in
    the cell, x*Shape(x) >= x.lo*Shape(x) in Loewner order because Shape(x) is
    PSD.  This avoids reintroducing the small-x interval dependency eliminated
    by the exact RL^{-1} process certificate.
    """
    shape = FACTORED.process_shape_q(x)
    shape_norm = _matrix_abs_row_sum_upper(shape)
    qnorm = BASE.up(x.hi * shape_norm)
    factor = BASE.down(1.0 / BASE.up(1.0 + BASE.up(beta * qnorm)))

    C = FACTORED._C_INTERVAL
    Ct = BASE.matrix_transpose(C)
    shape_tilde = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(C, shape), Ct)
    )
    omega_scale = BASE.down(factor * x.lo)
    Omega_tilde = _scale_matrix(shape_tilde, omega_scale)

    h = sched["dt_s"]
    scales = [sigma.lo * h, sigma.lo * h * h, sigma.lo * h * h * h, sigma.lo]
    upper = raw["Sigma_diagonal_upper"]
    physical = [upper[6], upper[9], upper[12], upper[15]]
    Sigma = _diag_matrix([
        BASE.up(physical[i] / (scales[i] * scales[i])) for i in range(4)
    ])
    Sigma_tilde = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(C, Sigma), Ct)
    )
    return Omega_tilde, Sigma_tilde, factor, qnorm


def _assemble_mode_matrices(mode: str, raw: dict, trans_Omega, trans_Sigma):
    n = 18 if mode == "H" else 21
    O = _zero_matrix(n)
    S = _zero_matrix(n)
    qpost = float(raw["post_measurement_scaled_Omega_lambda_min_lower"])
    scales2 = raw["comparison_scale_diagonal_squared"]
    upper = raw["Sigma_diagonal_upper"]

    trans_indices = set()
    for axis in range(3):
        idx = [6 + axis, 9 + axis, 12 + axis, 15 + axis]
        trans_indices.update(idx)
        for i in range(4):
            for j in range(4):
                O[idx[i]][idx[j]] = trans_Omega[i][j]
                S[idx[i]][idx[j]] = trans_Sigma[i][j]

    for i in range(n):
        if i in trans_indices:
            continue
        O[i][i] = Interval.outward_bounds(qpost, qpost)
        sui = BASE.up(upper[i] / scales2[i])
        S[i][i] = Interval.outward_bounds(sui, sui)
    return O, S


def mode_cell(mode: str, x: Interval, rho_trans: float, sigma: Interval,
              rs: Interval, live: dict, vector: dict, process: dict,
              sched: dict, alpha6: float) -> dict:
    raw = _ORIGINAL_MODE_CELL(
        mode, x, rho_trans, sigma, rs, live, vector, process, sched, alpha6
    )
    beta = _measurement_beta_upper(mode, sigma, rs, live, vector, process, sched)
    tO, tS, factor, qnorm = _translation_direct_blocks(x, sigma, raw, beta, sched)
    trans_delta = _certified_generalized_delta(tO, tS, MIN_USEFUL_DELTA)
    if trans_delta <= 0.0:
        raise RuntimeError(f"{mode} direct translation generalized inequality lost positivity")

    other = math.inf
    qpost = float(raw["post_measurement_scaled_Omega_lambda_min_lower"])
    trans_indices = {6+a for a in range(3)} | {9+a for a in range(3)} | {12+a for a in range(3)} | {15+a for a in range(3)}
    for i, (u, s2) in enumerate(zip(raw["Sigma_diagonal_upper"], raw["comparison_scale_diagonal_squared"])):
        if i in trans_indices:
            continue
        other = min(other, BASE.down(qpost / BASE.up(u / s2)))

    full_delta = BASE.down(min(trans_delta, other))
    if full_delta <= 0.0:
        raise RuntimeError(f"{mode} direct generalized matrix inequality lost positivity")
    fullO, fullS = _assemble_mode_matrices(mode, raw, tO, tS)
    full_ok = _spd_at_delta(fullO, fullS, full_delta)
    if not full_ok:
        raise RuntimeError(f"{mode} reported direct generalized delta did not re-certify")

    limiting = "translation_RL_inverse_block" if trans_delta <= other else "attitude_bias_or_active_ba_block"

    out = dict(raw)
    out["relative_Riccati_injection_margin_lower"] = full_delta
    out["direct_translation_generalized_margin_lower"] = trans_delta
    out["direct_nontranslation_margin_lower"] = other
    out["generalized_matrix_inequality"] = {
        "form": "Omega_word_RL_inverse_minus_delta_Sigma_upper_RL_inverse_is_SPD",
        "validated_interval_ldlt": True,
        "full_mode_dimension": len(fullO),
        "reported_delta_recertified": full_ok,
        "limiting_block": limiting,
        "measurement_information_beta_upper": beta,
        "translation_Qscaled_row_sum_norm_upper": qnorm,
        "translation_posterior_matrix_factor_lower": factor,
        "translation_process_representation": "x_lo*C*(Q_scaled/x)*C^T",
        "translation_congruence": "C=R*L_inverse applied to both Omega and Sigma",
        "full_margin_composition": "min(certified_translation_block, certified_nontranslation_diagonal_blocks), followed by full H/A LDLT recertification",
        "old_scalar_rho_over_max_scaled_upper_used": False,
    }
    return out


BASE.mode_cell = mode_cell
BASE._build_cached.cache_clear()


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = FACTORED.build(domain_path)
    out = dict(d)
    out["generalized_matrix_backend"] = "DIRECT_RL_INVERSE_CONGRUENCE_INTERVAL_LDLT"
    out["old_scalar_generalized_ratio_used"] = False
    return out


def validate(d: dict) -> list[str]:
    failures = FACTORED.validate(d)
    if d.get("generalized_matrix_backend") != "DIRECT_RL_INVERSE_CONGRUENCE_INTERVAL_LDLT":
        failures.append("direct generalized matrix backend missing")
    if d.get("old_scalar_generalized_ratio_used") is not False:
        failures.append("old scalar generalized ratio still active")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {}).get("matrix_comparison", {})
        g = m.get("generalized_matrix_inequality", {})
        if g.get("validated_interval_ldlt") is not True:
            failures.append(f"{mode} direct generalized inequality not LDLT certified")
        if g.get("reported_delta_recertified") is not True:
            failures.append(f"{mode} direct generalized delta not recertified")
        if g.get("old_scalar_rho_over_max_scaled_upper_used") is not False:
            failures.append(f"{mode} old scalar generalized ratio still active")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "backend": d["generalized_matrix_backend"],
        "cells": d["cell_partition"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
