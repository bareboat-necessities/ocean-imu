#!/usr/bin/env python3
"""Numerically stable backend for the canonical SEA3 moving-Riccati tube.

This is not a second theorem route.  It delegates the complete tube construction
and validation to ``ou3_sea3_riccati_tube`` and replaces one algebraically
identical numerical primitive: the scaled integrated-OU process covariance on
the implementation's small-x branch.

The original implementation first encloses Qbar(x) and then divides interval
entries by powers of the *same* interval x to form

    D^-1 Q D^-T,  D=diag(sigma*h, sigma*h^2, sigma*h^3, sigma).

At the long-tau edge x=h/tau ~= 4.17e-4, that introduces avoidable interval
dependency and can lose positive definiteness even after deep subdivision.
Here each small-x series is divided by x^k symbolically first, leaving a direct
polynomial with leading order x.  The polynomial is exactly the same truncated
Maclaurin branch already audited in the canonical producer; only its interval
evaluation order changes.
"""
from __future__ import annotations

import ou3_sea3_riccati_tube as BASE


_BASE_STEP_SCALED_Q = BASE.step_scaled_q


def _p(x, terms):
    return BASE.poly(x, tuple(terms))


def _small_scaled_q_factored(x):
    """Exact symbolic x-power cancellation for BASE's small-x Qbar series."""
    qvv = _p(x, ((1,2/3),(2,-1/2),(3,7/30),(4,-1/12),(5,31/1260),(6,-1/160),(7,127/90720)))
    qvp = _p(x, ((1,1/4),(2,-1/6),(3,5/72),(4,-1/45),(5,17/2880),(6,-41/30240)))
    qvS = _p(x, ((1,1/15),(2,-1/24),(3,41/2520),(4,-7/1440),(5,109/90720)))
    qva = _p(x, ((1,1),(2,-1),(3,7/12),(4,-1/4),(5,31/360),(6,-1/40),(7,127/20160),(8,-17/12096)))

    qpp = _p(x, ((1,1/10),(2,-1/18),(3,5/252),(4,-1/180),(5,17/12960)))
    qpS = _p(x, ((1,1/36),(2,-1/72),(3,13/2880),(4,-1/864)))
    qpa = _p(x, ((1,1/3),(2,-1/3),(3,11/60),(4,-13/180),(5,19/840),(6,-1/168),(7,247/181440)))

    qSS = _p(x, ((1,1/126),(2,-1/288),(3,13/12960)))
    qSa = _p(x, ((1,1/12),(2,-1/12),(3,2/45),(4,-1/60),(5,11/2240),(6,-73/60480)))

    qaa = _p(x, ((1,2),(2,-2),(3,4/3),(4,-2/3),(5,4/15),(6,-4/45),(7,8/315),(8,-2/315),(9,4/2835)))

    return [
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ]


def step_scaled_q(x):
    if x.hi < BASE.BRANCH_X:
        return _small_scaled_q_factored(x)
    return _BASE_STEP_SCALED_Q(x)


# Patch only the numerical primitive used by BASE.build/split_x_cell.  All
# source contracts, cell construction, covariance ceilings, Joseph comparison,
# H/A aggregation, validation, and CLI semantics remain BASE's implementation.
BASE.step_scaled_q = step_scaled_q

SCHEMA = BASE.SCHEMA
QUALIFICATION = BASE.QUALIFICATION
USEFUL_GATE = BASE.USEFUL_GATE
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN


def build(domain_path=DEFAULT_DOMAIN):
    return BASE.build(domain_path)


def validate(payload):
    return BASE.validate(payload)


def main():
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
