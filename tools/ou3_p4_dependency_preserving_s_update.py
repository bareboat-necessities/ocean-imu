#!/usr/bin/env python3
"""Dependency-preserving homogeneous S=0 update for OU-III P4.

For the shipping S=0 pseudo-measurement in homogeneous physical error
coordinates, H is the exact selector of state entries 12:15 and r=H z.  The
accepted update is

    z+ = z - P H^T (H P H^T + R)^-1 r.

Writing S_i = H P H^T + R gives the exact measured-state identity

    H z+ = r - (S_i-R) S_i^-1 r = R S_i^-1 r.

Evaluating ``z_S - dx_S`` with independent interval boxes destroys this shared
residual dependence when the S gain is near identity.  In the current word it
turns a strongly contracting sample-198 update into a fictitious roughly
[-614,614] sample-199 residual.  This module replaces only the three measured
S-state values/derivatives by the algebraically identical ``R S_i^-1 r`` form.

Everything else remains unchanged: the same joint P,H,R,r innovation solve,
the same direct attitude correction, signed Joseph form, posterior covariance,
and reset congruence are used.  No K interval matrix is formed.

Scope is intentionally explicit.  The identity below is the homogeneous map
with physical S_true=0.  In the physical true-minus-estimate error system,
nonzero S_true is an external forcing term and must be charged to the eventual
complete-word affine b.  This helper therefore cannot by itself promote P4/P5.
"""
from __future__ import annotations

import json

import ou3_interval_ad as AD
import ou3_p4_joint_joseph as JJ
import ou3_p4_joint_word_dissipation_design as D

SCHEMA = 1
_ORIGINAL = D._ad_joint_update


def _is_point(x, value: float) -> bool:
    v = float(value)
    return float(x.lo) == v and float(x.hi) == v


def is_exact_S_selector(Hm, state_dimension: int) -> bool:
    if len(Hm) != 3 or any(len(row) != int(state_dimension) for row in Hm):
        return False
    if int(state_dimension) not in (18, 21):
        return False
    for i in range(3):
        for j in range(int(state_dimension)):
            target = 1.0 if j == 12 + i else 0.0
            if not _is_point(Hm[i][j], target):
                return False
    return True


def dependency_preserving_joint_update(Pm, z, Hm, Rm, residual, eta):
    """Route exact S selectors through H z+ = R S^-1 r; delegate all others."""
    n = len(z)
    if not is_exact_S_selector(Hm, n):
        return _ORIGINAL(Pm, z, Hm, Rm, residual, eta)

    PHt, Sinnov = JJ.innovation(Pm, Hm, Rm)
    Sinv, meta = JJ.verified_inverse(Sinnov)
    sol = D._ad_matvec(Sinv, residual)
    dx = D._ad_matvec(PHt, sol)
    signed = D._signed_form(residual, eta, Sinv, Rm)

    Pj = JJ.posterior_covariance(Pm, PHt, Sinv, psd_tighten=D.H._psd_tighten)
    out = list(z)
    d = [-q for q in dx[:3]]
    out[:3] = AD.deployed_correct_cayley_right(z[:3], d)

    # Unmeasured state rows retain the direct source map.  The three measured
    # S rows use the exact same-cell cancellation identity instead of natural
    # interval subtraction of two nearly equal, highly dependent boxes.
    for i in range(3, n):
        out[i] = z[i] - dx[i]
    rplus = D._ad_matvec(Rm, sol)
    for i in range(3):
        out[12 + i] = rplus[i]

    Pout = D.H._reset_covariance(Pj, [q.val for q in dx[:3]])
    return Pout, out, signed, {
        "inverse_q_inf_upper": meta.get("neumann_q_inf_upper"),
        "correction_theta_norm_upper": AD._norm_upper([q.val for q in dx[:3]]),
        "K_interval_matrix_materialized": False,
        "dependency_preserving_S_selector_identity_used": True,
        "measured_state_identity": "H z_plus = R (H P H^T + R)^-1 r",
        "homogeneous_S_true_zero_only": True,
        "physical_S_true_forcing_must_enter_affine_b": True,
    }


def install() -> None:
    D._ad_joint_update = dependency_preserving_joint_update


def _self_test() -> dict:
    # Structural guard: exact selector must be recognized; a single altered
    # element must be rejected so accel/mag updates cannot silently use it.
    n = 18
    Z = D._zero(3, n)
    for i in range(3):
        Z[i][12 + i] = D.I(1.0)
    exact = is_exact_S_selector(Z, n)
    Zbad = [[x for x in row] for row in Z]
    Zbad[0][0] = D.I(1.0)
    rejected_bad = not is_exact_S_selector(Zbad, n)
    return {
        "schema": SCHEMA,
        "pass": bool(exact and rejected_bad),
        "exact_selector_recognized": bool(exact),
        "nonselector_rejected": bool(rejected_bad),
        "identity": "H z_plus = R Sinnov^-1 r",
        "homogeneous_S_true_zero_only": True,
        "physical_S_true_forcing_must_enter_affine_b": True,
        "K_interval_matrix_materialized": False,
        "P4_PROMOTED_HERE": False,
        "P5_PROMOTED_HERE": False,
    }


# Install on import for the co-gauged proof backends that explicitly import this
# module through ou3_p4_joint_word_gauge_design_v2.
install()


if __name__ == "__main__":
    d = _self_test()
    print(json.dumps(d, indent=2, sort_keys=True))
    raise SystemExit(0 if d["pass"] else 2)
