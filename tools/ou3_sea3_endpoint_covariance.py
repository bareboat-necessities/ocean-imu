#!/usr/bin/env python3
"""Endpoint-referenced SEA3 translation covariance ceiling for canonical P3.

The finite-memory S=0 observation argument reconstructs the error state at the
*end* of the covariance word.  If ``ell`` is the lag of a selected S packet,

    S(T-ell) = S(T) - ell p(T) + 0.5 ell^2 v(T) + disturbance.

Therefore the existing forward-time integrator observation matrix differs only
by the sign similarity ``diag(1,-1,1)`` on ``[S,p,v]``.  Once that similarity
is applied, the reconstructed ``[v,p,S]`` covariance already lives at the word
endpoint and must not be propagated by another full word.

The S-observation and vector-PE memories are also endpoint-referenced memories,
not serial phases of an execution.  A uniform S firing-gap upper ``g`` lets us
select one firing from each of the disjoint lag windows

    [0,g], [2g,3g], [4g,5g].

Those windows are separated by at least ``g`` for every legal choice.  The
inverse of the resulting quadratic Vandermonde matrix is evaluated here from
its exact Lagrange formula.  This matters for interval arithmetic: generic
Gauss-Jordan elimination forgets that repeated occurrences of the same lag are
dependent and can create a false zero-crossing pivot even though the three lag
windows are strictly ordered.  The closed-form rational expression preserves
that ordering and is outward evaluated using the repository Interval layer.

Consequently the translation observation memory is only ``5g``; it does not
need the 1 s vector PE recurrence as artificial spacing between S packets.  The
vector pair and the three S packets may live in overlapping backward windows
ending at the same Riccati endpoint, so the common covariance memory is
``max(5g,T_PE)``, not their sum.  Cross-block covariance is still paid by the
existing trace/diagonal Loewner dominators.

This helper deliberately keeps the #489 source-side architecture intact:

* the global SEA3 dynamic invariant, not a P2 history graph;
* the progress-preserving pseudo-update recurrence bound;
* the same nuisance/process and R_S upper bounds;
* the same full-memory diagonal process/noise dominator; and
* arbitrary time variation of the adaptive source inside the memory window.

No filter constant, theorem domain, usefulness gate, or source-history
assumption changes.
"""
from __future__ import annotations

import ou3_sea3_riccati_tube as BASE


QUALIFICATION = "OU3_SEA3_ENDPOINT_REFERENCED_TRANSLATION_COVARIANCE"


def _endpoint_window_integrator_inverse(gap: float):
    """Outward inverse for rows ``[1,t,t^2/2]`` on the three guaranteed windows.

    Let V have rows ``[1,t,t^2]``.  Its inverse is the coefficient matrix of
    the three Lagrange polynomials.  Since the physical observation matrix is
    ``B = V diag(1,1,1/2)``, ``B^-1 = diag(1,1,2) V^-1``.  Evaluating that
    identity directly avoids the dependency loss of interval elimination.
    """
    g = BASE.pos(gap, "gap")
    t0 = BASE.Interval.outward_bounds(0.0, g)
    t1 = BASE.Interval.outward_bounds(2.0 * g, 3.0 * g)
    t2 = BASE.Interval.outward_bounds(4.0 * g, 5.0 * g)

    d0 = (t0 - t1) * (t0 - t2)
    d1 = (t1 - t0) * (t1 - t2)
    d2 = (t2 - t0) * (t2 - t1)
    for i, d in enumerate((d0, d1, d2)):
        if d.lo <= 0.0 <= d.hi:
            raise RuntimeError(f"endpoint Vandermonde denominator {i} lost separation")

    return [
        [
            (t1 * t2) / d0,
            (t0 * t2) / d1,
            (t0 * t1) / d2,
        ],
        [
            -(t1 + t2) / d0,
            -(t0 + t2) / d1,
            -(t0 + t1) / d2,
        ],
        [
            BASE.I(2.0) / d0,
            BASE.I(2.0) / d1,
            BASE.I(2.0) / d2,
        ],
    ]


def global_translation_upper(dynamic: dict, live: dict,
                             axis_factors: list[float]) -> tuple[list[float], dict]:
    """Return the source-uniform endpoint covariance ceiling without serial memory."""
    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = BASE.pos(rates["dt_s"], "dt")
    tau_lo = BASE.pos(inv["tau_applied_s"][0], "tau lower")
    sigma_hi = BASE.pos(inv["sigma_aw_filter_mps2"][1], "sigma upper")
    rs_hi = BASE.pos(inv["R_S_applied"][1], "R_S upper")
    cadence_hi = BASE.pos(inv["pseudo_update_period_s"][1], "pseudo cadence upper")
    Tpe = BASE.pos(live["vector_pe_recurrence_window_s"], "PE recurrence")

    # The scheduler can miss a newly reached deadline by at most one configured
    # sample.  With max firing gap g, every lag interval of width g contains an
    # S packet.  Choosing three width-g intervals at 0, 2g and 4g guarantees
    # pairwise lag separation >=g without importing the unrelated vector-PE
    # recurrence into the translation observability matrix.
    gap = BASE.up(cadence_hi + h)
    spacing = BASE.up(2.0 * gap)
    Tobs = BASE.up(2.0 * spacing + gap)  # 5g

    # Both certificates reconstruct the *same endpoint*.  Their backward
    # evidence windows can overlap; serializing them paid a non-existent extra
    # T_PE interval.  The common full-state covariance memory is their union.
    Tword = BASE.up(max(Tobs, Tpe))

    Binv = _endpoint_window_integrator_inverse(gap)
    qc_hi = BASE.up(2.0 * sigma_hi * sigma_hi / tau_lo)
    s_nuis = BASE.up(sigma_hi * sigma_hi * (Tobs ** 3 / 6.0) ** 2)
    s_proc = BASE.up(qc_hi * Tobs ** 7 / 252.0)
    rmax = BASE.up((rs_hi * max(axis_factors)) ** 2)
    rstack = BASE.up(3.0 * (rmax + s_nuis + s_proc))
    R = [
        [BASE.I(rstack if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]

    # The closed-form inverse above is for rows [1,+ell,0.5 ell^2] in
    # [S,p,v].  Endpoint lags require [1,-ell,0.5 ell^2], which is exactly the
    # sign similarity below.  No subsequent Phi(Tword) propagation is
    # admissible: this covariance already refers to x(T).
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

    # Retain #489's source-uniform process/noise Loewner dominator, now over
    # the actual common backward memory rather than a serial concatenation.
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
        "S_observation_window_spacing_s": spacing,
        "observation_window_s_upper": Tobs,
        "vector_PE_window_s_upper": Tpe,
        "covariance_memory_window_s_upper": Tword,
        "q_c_global_upper": qc_hi,
        "source_motion_inside_window_allowed": True,
        "translation_reference": "word_endpoint",
        "endpoint_referenced_observability": True,
        "endpoint_p_sign_similarity_applied": True,
        "forward_propagation_after_endpoint_reconstruction": False,
        "S_observation_spacing_uses_vector_PE": False,
        "S_observation_window_layout": "[0,g],[2g,3g],[4g,5g]",
        "S_observation_inverse_route": "OUTWARD_CLOSED_FORM_LAGRANGE_VANDERMONDE",
        "generic_interval_gauss_jordan_used_for_S_observation": False,
        "S_and_vector_PE_memories_overlap_at_endpoint": True,
        "covariance_memory_is_max_not_sum": True,
        "full_word_process_noise_dominator_retained": True,
        "qualification": QUALIFICATION,
    }


def install(base=BASE) -> None:
    """Install the endpoint theorem into the single canonical BASE build."""
    if base is not BASE:
        raise RuntimeError("endpoint covariance patch must target canonical SEA3 tube module")
    base._global_translation_upper = global_translation_upper
