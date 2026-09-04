#!/usr/bin/env python3
"""Endpoint-referenced SEA3 translation covariance ceiling for canonical P3.

The finite-memory S=0 observation argument reconstructs the error state at the
*end* of the covariance word.  If ``ell`` is the lag of a selected S packet,

    S(T-ell) = S(T) - ell p(T) + 0.5 ell^2 v(T) + disturbance.

Therefore the existing forward-time integrator observation matrix differs only
by the sign similarity ``diag(1,-1,1)`` on ``[S,p,v]``.  Once that similarity
is applied, the reconstructed ``[v,p,S]`` covariance already lives at the word
endpoint and must not be propagated by another full word.  Doing so pays the
integrator powers twice and is not a property of the shipping Riccati filter.

This helper deliberately keeps every #489 source-side choice intact:

* the global SEA3 dynamic invariant, not a P2 history graph;
* the progress-preserving pseudo-update recurrence bound;
* the same nuisance/process and R_S upper bounds;
* the same full-word diagonal process/noise dominator; and
* arbitrary time variation of the adaptive source inside the memory window.

Only the covariance reference time is corrected.  No filter constant, theorem
domain, usefulness gate, or source/history assumption changes.
"""
from __future__ import annotations

import ou3_sea3_riccati_tube as BASE


QUALIFICATION = "OU3_SEA3_ENDPOINT_REFERENCED_TRANSLATION_COVARIANCE"


def global_translation_upper(dynamic: dict, live: dict,
                             axis_factors: list[float]) -> tuple[list[float], dict]:
    """Return the #489 global ceiling with endpoint-referenced S observability."""
    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = BASE.pos(rates["dt_s"], "dt")
    tau_lo = BASE.pos(inv["tau_applied_s"][0], "tau lower")
    sigma_hi = BASE.pos(inv["sigma_aw_filter_mps2"][1], "sigma upper")
    rs_hi = BASE.pos(inv["R_S_applied"][1], "R_S upper")
    cadence_hi = BASE.pos(inv["pseudo_update_period_s"][1], "pseudo cadence upper")
    Tpe = BASE.pos(live["vector_pe_recurrence_window_s"], "PE recurrence")

    gap = BASE.up(cadence_hi + h)
    spacing = BASE.up(max(Tpe, 2.0 * gap))
    Tobs = BASE.up(2.0 * spacing + gap)
    Tword = BASE.up(Tobs + Tpe)

    Binv = BASE.integrator_inverse(gap, spacing)
    qc_hi = BASE.up(2.0 * sigma_hi * sigma_hi / tau_lo)
    s_nuis = BASE.up(sigma_hi * sigma_hi * (Tobs ** 3 / 6.0) ** 2)
    s_proc = BASE.up(qc_hi * Tobs ** 7 / 252.0)
    rmax = BASE.up((rs_hi * max(axis_factors)) ** 2)
    rstack = BASE.up(3.0 * (rmax + s_nuis + s_proc))
    R = [
        [BASE.I(rstack if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]

    # BASE.integrator_inverse is for rows [1,+ell,0.5 ell^2] in [S,p,v].
    # Endpoint lags require [1,-ell,0.5 ell^2], which is exactly the sign
    # similarity below.  No subsequent Phi(Tword) propagation is admissible:
    # this covariance already refers to x(T).
    Cspv_forward_sign = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(Binv, R), BASE.matrix_transpose(Binv))
    )
    sign = [
        [
            BASE.I(-1.0 if i == j == 1 else (1.0 if i == j else 0.0))
            for j in range(3)
        ]
        for i in range(3)
    ]
    Cspv_endpoint = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(
            BASE.matrix_mul(sign, Cspv_forward_sign),
            BASE.matrix_transpose(sign),
        )
    )
    order = (2, 1, 0)  # [S,p,v] -> [v,p,S]
    Cvps_endpoint = [
        [Cspv_endpoint[order[i]][order[j]] for j in range(3)]
        for i in range(3)
    ]
    u = BASE.diagonal_dominator(Cvps_endpoint)

    # Retain #489's full-word process/noise Loewner dominator unchanged.  This
    # is intentionally conservative; the present theorem correction removes
    # only the artificial second deterministic propagation of the estimator.
    variances = [
        BASE.up(sigma_hi * sigma_hi * Tword * Tword + qc_hi * Tword ** 3 / 3.0),
        BASE.up(sigma_hi * sigma_hi * Tword ** 4 / 4.0 + qc_hi * Tword ** 5 / 20.0),
        BASE.up(sigma_hi * sigma_hi * Tword ** 6 / 36.0 + qc_hi * Tword ** 7 / 252.0),
        BASE.up(sigma_hi * sigma_hi),
    ]
    roots = [BASE.math.sqrt(v) for v in variances]
    total = BASE.up(sum(roots))
    noise = [BASE.up(r * total) for r in roots]
    upper = [BASE.up(u[i] + noise[i]) for i in range(3)] + [noise[3]]

    return upper, {
        "pseudo_gap_s_upper": gap,
        "observation_window_s_upper": Tobs,
        "covariance_memory_window_s_upper": Tword,
        "q_c_global_upper": qc_hi,
        "source_motion_inside_window_allowed": True,
        "translation_reference": "word_endpoint",
        "endpoint_referenced_observability": True,
        "endpoint_p_sign_similarity_applied": True,
        "forward_propagation_after_endpoint_reconstruction": False,
        "full_word_process_noise_dominator_retained": True,
        "qualification": QUALIFICATION,
    }


def install(base=BASE) -> None:
    """Install the endpoint theorem into the single canonical BASE build."""
    if base is not BASE:
        raise RuntimeError("endpoint covariance patch must target canonical SEA3 tube module")
    base._global_translation_upper = global_translation_upper
