#!/usr/bin/env python3
"""Binary32 arithmetic contract for the shipping OU-III prediction/Joseph/reset path.

This module separates finite-precision obligations that are genuinely closed by
explicit scalar-operation structure from the remaining solver/covariance
feedback obligations.  It is deliberately fail-closed: the existing projection
rounding certificate MUST NOT be relabelled as a certificate for the complete
Kalman event.

Shipping structure audited here:

* measurement innovation is factored by ``safe_ldlt3_`` using Eigen LDLT;
* gain rows are obtained by 3x3 ``ldlt.solve(rhs)`` calls;
* state correction is ``xext += K*r``;
* Joseph covariance is the scalar-loop identity
      P <- P - KCP - (KCP)^T + K S K^T;
* quaternion correction is injected on the left, normalized, and followed by
  the explicit covariance reset ``G=I+0.5*skew(dtheta)``;
* accelerometer-bias ball projection has its own already-closed binary32 value
  enclosure.

For explicit add/multiply/FMA/dot/reset operations we can use an absolute
rounding model valid through cancellation:

    |RN32(x)-x| <= u |x| + eta,

with u=2^-24 and eta=2^-150 (half the smallest binary32 subnormal spacing).
This does not require a relative-error assumption on a near-zero result.
Helpers below return conservative absolute local error envelopes for dot3,
state correction, Joseph scalar entries and the covariance reset once operand
magnitude bounds are supplied by a same-history event cell.

Two obligations remain before the complete filter finite-precision map is
closed:

1. prove/source-audit a backward-error enclosure for the concrete Eigen 3x3
   LDLT factor/solve path, including exclusion (or explicit modelling) of the
   ``safe_ldlt3_`` diagonal-bump fallback on the admitted source family;
2. recursively propagate the covariance rounding enclosure into the P used by
   the NEXT event's innovation/gain.  Bounding one Joseph event in isolation is
   insufficient because P is part of the future same-history coefficient cell.

Until both are attached to every selector prefix, ``gamma_n`` at the terminal
augmented LDLT may not claim complete binary32 implementation coverage.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_projection_binary32_enclosure as PROJECTION

REPO = Path(__file__).resolve().parents[2]
FILTER = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
SCHEMA = 1
QUALIFICATION = "OU3_P4_FULL_FILTER_BINARY32_EVENT_ARITHMETIC_CONTRACT_V1"
U = 2.0 ** -24
ETA = 2.0 ** -150


def _finite_nonnegative(x: float, name: str) -> float:
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError(f"{name} must be finite nonnegative")
    return x


def rn_abs_error(magnitude_upper: float) -> float:
    """One binary32 rounding, including the subnormal absolute-error floor."""
    m = _finite_nonnegative(magnitude_upper, "magnitude_upper")
    return math.nextafter(U * m + ETA, math.inf)


def dot3_abs_error(a_abs_upper: float, b_abs_upper: float) -> float:
    """Conservative error for a 3-term dot with rounded products/additions.

    We deliberately model three product roundings plus two additions.  FMA can
    only reduce this count.  Intermediate absolute sums are <=3*a*b.
    """
    a = _finite_nonnegative(a_abs_upper, "a_abs_upper")
    b = _finite_nonnegative(b_abs_upper, "b_abs_upper")
    prod = a * b
    e_prod = rn_abs_error(prod)
    # Each rounded product magnitude is bounded by prod+e_prod.
    p = prod + e_prod
    e_add1 = rn_abs_error(2.0 * p)
    e_add2 = rn_abs_error(3.0 * p + e_add1)
    return math.nextafter(3.0 * e_prod + e_add1 + e_add2, math.inf)


def state_correction_component_error(k_abs_upper: float, residual_abs_upper: float,
                                     state_abs_upper: float) -> float:
    """Local arithmetic error in one x_i <- x_i + dot3(K_i,r) component."""
    k = _finite_nonnegative(k_abs_upper, "k_abs_upper")
    r = _finite_nonnegative(residual_abs_upper, "residual_abs_upper")
    x = _finite_nonnegative(state_abs_upper, "state_abs_upper")
    dot_exact = 3.0 * k * r
    dot_err = dot3_abs_error(k, r)
    add_err = rn_abs_error(x + dot_exact + dot_err)
    return math.nextafter(dot_err + add_err, math.inf)


def joseph_entry_local_error(p_abs_upper: float, k_abs_upper: float,
                             pct_abs_upper: float, s_abs_upper: float) -> float:
    """Local scalar-loop error for one Joseph covariance entry.

    The exact identity is P_ij - dot3(K_i,PCt_j) - dot3(PCt_i,K_j)
    + K_i S K_j^T.  The last term is three dot3 contributions after a
    3x3-by-vector multiply, conservatively overcounted here.  This is ONLY the
    arithmetic of one event with exact operands; operand errors and next-event
    covariance feedback are separate.
    """
    p = _finite_nonnegative(p_abs_upper, "p_abs_upper")
    k = _finite_nonnegative(k_abs_upper, "k_abs_upper")
    c = _finite_nonnegative(pct_abs_upper, "pct_abs_upper")
    s = _finite_nonnegative(s_abs_upper, "s_abs_upper")
    e_kc = dot3_abs_error(k, c)
    # Each component of S*K_j has exact magnitude <=3*s*k.
    sk_mag = 3.0 * s * k
    e_sk = dot3_abs_error(s, k)
    e_ksk = dot3_abs_error(k, sk_mag + e_sk)
    exact_mag = p + 6.0 * k * c + 9.0 * k * k * s
    # Three signed scalar accumulations after the constituent dot products.
    e_acc = 3.0 * rn_abs_error(exact_mag + 2.0 * e_kc + e_sk + e_ksk)
    return math.nextafter(2.0 * e_kc + e_sk + e_ksk + e_acc, math.inf)


def reset_covariance_entry_local_error(p_abs_upper: float, dtheta_abs_upper: float) -> float:
    """Local error of the explicit left-error covariance reset scalar loops."""
    p = _finite_nonnegative(p_abs_upper, "p_abs_upper")
    d = _finite_nonnegative(dtheta_abs_upper, "dtheta_abs_upper")
    # Every G entry has |G_ij| <= 1 + d/2.  Paa uses G*P*G^T; cross rows G*P.
    g = 1.0 + 0.5 * d
    e_gp = dot3_abs_error(g, p)
    gp = 3.0 * g * p + e_gp
    e_gpg = dot3_abs_error(gp, g)
    # Final off-diagonal symmetrization adds and halves two already-rounded values.
    sym_mag = 3.0 * gp * g + e_gpg
    e_sym = rn_abs_error(2.0 * sym_mag) + rn_abs_error(sym_mag)
    return math.nextafter(e_gp + e_gpg + e_sym, math.inf)


def _source_audit() -> dict:
    f = FILTER.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    checks = {
        "safe_ldlt3_factor_present": "safe_ldlt3_" in f and "ldlt.compute(S)" in f,
        "gain_row_solve_present": "gain_from_ldlt3_" in f and "ldlt.solve(rhs)" in f,
        "state_correction_present": "xext.noalias() += K * r" in f,
        "joseph_scalar_identity_present": (
            "P ← P - KCP - (KCP)ᵀ + K S Kᵀ" in f
            or "P <- P - KCP" in f
            or "joseph_update3_" in f
        ),
        "left_quaternion_correction_present": "qref = corr * qref" in f and "qref.normalize()" in f,
        "explicit_reset_G_present": "Identity() + T(0.5)*skew(dtheta)" in c,
        "explicit_reset_scalar_dot_loops_present": "sum += G(i,k)*Paa_old(k,j)" in c,
        "safe_ldlt_diagonal_bump_branch_present": "S.diagonal().array() += bump" in f,
    }
    return checks


def build() -> dict:
    projection = PROJECTION.build()
    pf = PROJECTION.validate(projection)
    if pf:
        raise RuntimeError(f"projection binary32 prerequisite failed: {pf}")
    audit = _source_audit()
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "runtime_scalar_format": "IEEE754_binary32",
        "unit_roundoff": U,
        "half_min_subnormal_spacing": ETA,
        "source_audit": audit,
        "source_audit_closed": all(audit.values()),
        "cancellation_safe_absolute_rounding_primitive_available": True,
        "dot3_local_rounding_enclosure_available": True,
        "state_correction_local_rounding_enclosure_available": True,
        "joseph_scalar_loop_local_rounding_enclosure_available": True,
        "left_reset_scalar_loop_local_rounding_enclosure_available": True,
        "projection_binary32_enclosure_consumed": True,
        "projection_finite_precision_enclosure_closed": bool(
            projection["projection_finite_precision_enclosure_closed"]
        ),
        "eigen_3x3_LDLT_backward_error_enclosure_closed": False,
        "safe_ldlt_first_factorization_success_source_uniformly_proved": False,
        "safe_ldlt_diagonal_bump_modelled_in_exact_event_graph": False,
        "covariance_roundoff_recursively_attached_to_next_event_P": False,
        "future_gain_sensitivity_to_covariance_roundoff_closed": False,
        "full_prediction_arithmetic_enclosure_closed_here": False,
        "full_Kalman_reset_finite_precision_enclosure_closed": False,
        "terminal_gamma_n_may_claim_full_shipping_binary32": False,
        "packet_count_times_worst_roundoff_used": False,
        "filter_changed": False,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "next_obligation": (
            "close the source-uniform Eigen 3x3 LDLT/fallback branch and recursively propagate Joseph/reset "
            "covariance roundoff into the next event P; then suffix-propagate the resulting physical-state "
            "roundoff channels through every-event PrefixInput instead of a packet-count worst-case sum"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("runtime_scalar_format") != "IEEE754_binary32":
        f.append("runtime format changed")
    for k in (
        "source_audit_closed",
        "cancellation_safe_absolute_rounding_primitive_available",
        "dot3_local_rounding_enclosure_available",
        "state_correction_local_rounding_enclosure_available",
        "joseph_scalar_loop_local_rounding_enclosure_available",
        "left_reset_scalar_loop_local_rounding_enclosure_available",
        "projection_binary32_enclosure_consumed",
        "projection_finite_precision_enclosure_closed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "eigen_3x3_LDLT_backward_error_enclosure_closed",
        "safe_ldlt_first_factorization_success_source_uniformly_proved",
        "safe_ldlt_diagonal_bump_modelled_in_exact_event_graph",
        "covariance_roundoff_recursively_attached_to_next_event_P",
        "future_gain_sensitivity_to_covariance_roundoff_closed",
        "full_prediction_arithmetic_enclosure_closed_here",
        "full_Kalman_reset_finite_precision_enclosure_closed",
        "terminal_gamma_n_may_claim_full_shipping_binary32",
        "packet_count_times_worst_roundoff_used",
        "filter_changed",
        "P4_MOTION_PASS", "P4_PASS",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    # Smoke the arithmetic helpers on representative finite operands.
    vals = (
        rn_abs_error(1.0), dot3_abs_error(2.0, 3.0),
        state_correction_component_error(2.0, 3.0, 4.0),
        joseph_entry_local_error(4.0, 2.0, 3.0, 5.0),
        reset_covariance_entry_local_error(4.0, 0.4),
    )
    if not all(math.isfinite(x) and x > 0.0 for x in vals):
        f.append("local binary32 arithmetic helper produced invalid bound")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(); f = validate(d)
    d["validation_pass"] = not f; d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "local_scalar_arithmetic": d["joseph_scalar_loop_local_rounding_enclosure_available"],
        "ldlt": d["eigen_3x3_LDLT_backward_error_enclosure_closed"],
        "recursive_covariance": d["covariance_roundoff_recursively_attached_to_next_event_P"],
        "full_binary32": d["full_Kalman_reset_finite_precision_enclosure_closed"],
        "P4": d["P4_PASS"], "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
