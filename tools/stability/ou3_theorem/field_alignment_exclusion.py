"""Exact scalar checks for the conditional OU-III field-axis exclusion lemma.

This module is NOT an end-to-end stability certificate. All inputs describe
one retained, accepted-update execution. It does not assert that the execution
enters/stays in the region or that its learned magnetic reference is accurate.
Source audit: PR #637, 1bbd3b3bc58879d7f95863cef5726ee10c80616f.
"""
from fractions import Fraction
from math import factorial
from typing import Sequence

F = Fraction


def _nonnegative(value: Fraction, name: str) -> Fraction:
    result = F(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def trapezoidal_weights(gaps: Sequence[Fraction]) -> tuple[Fraction, ...]:
    """Weights on the actual sample endpoints; not a runtime resampler."""
    steps = tuple(F(gap) for gap in gaps)
    if not steps or any(gap <= 0 for gap in steps):
        raise ValueError("at least one strictly positive accepted-update gap is required")
    length = sum(steps, F(0))
    return (steps[0] / (2 * length),) + tuple(
        (left + right) / (2 * length)
        for left, right in zip(steps[:-1], steps[1:])
    ) + (steps[-1] / (2 * length),)


def physical_mean_charge(*, duration: Fraction, velocity_bound: Fraction,
                         jerk_bound: Fraction, effective_gap: Fraction) -> Fraction:
    """2 V/L + J sum(gap^2)/(4 L), with an upper bound for effective_gap."""
    length = F(duration)
    if length <= 0:
        raise ValueError("duration must be strictly positive")
    velocity = _nonnegative(velocity_bound, "velocity_bound")
    jerk = _nonnegative(jerk_bound, "jerk_bound")
    gap = _nonnegative(effective_gap, "effective_gap")
    return 2 * velocity / length + jerk * gap / 4


def geometric_margin(*, transverse_gravity_lower: Fraction, duration: Fraction,
                     velocity_bound: Fraction, jerk_bound: Fraction,
                     effective_gap: Fraction, acceleration_error_factor: Fraction,
                     retained_radius: Fraction, force_defect: Fraction) -> Fraction:
    """Conditional lower bound for weighted nominal transverse-force norm.

    force_defect is mandatory: missing reference/calibration bounds must not be
    silently interpreted as zero. For a varying reference it must include the
    full reference-direction charge derived in the accompanying proof note.
    """
    gravity = _nonnegative(transverse_gravity_lower, "transverse_gravity_lower")
    factor = _nonnegative(acceleration_error_factor, "acceleration_error_factor")
    radius = _nonnegative(retained_radius, "retained_radius")
    defect = _nonnegative(force_defect, "force_defect")
    return gravity - factor * radius - defect - physical_mean_charge(
        duration=duration, velocity_bound=velocity_bound, jerk_bound=jerk_bound,
        effective_gap=effective_gap,
    )


def small_branch_process_bound() -> Fraction:
    """Uniform real-arithmetic Q_aa/[sigma^2(1-phi^2)] upper bound.

    Absolute coefficients below are the literal source polynomials, with
    sigma^2 * h^p * x factored out, x=h/tau. Use h<=.006 and x<.01.
    The first 3x3 regularizer is included before the additional S row.
    Successful spectral clipping cannot increase the operator norm.
    Floating-point sanitization/eigensolver fallbacks are NOT certified here.
    """
    coefficients = {
        "aa": (0, [F(2), F(2), F(4, 3), F(2, 3), F(4, 15), F(4, 45), F(8, 315), F(2, 315), F(4, 2835)]),
        "va": (1, [F(1), F(1), F(7, 12), F(1, 4), F(31, 360), F(1, 40), F(127, 20160), F(17, 12096)]),
        "pa": (2, [F(1, 3), F(1, 3), F(11, 60), F(13, 180), F(19, 840), F(1, 168), F(247, 181440)]),
        "vv": (2, [F(2, 3), F(1, 2), F(7, 30), F(1, 12), F(31, 1260), F(1, 160), F(127, 90720)]),
        "vp": (3, [F(1, 4), F(1, 6), F(5, 72), F(1, 45), F(17, 2880), F(41, 30240)]),
        "pp": (4, [F(1, 10), F(1, 18), F(5, 252), F(1, 180), F(17, 12960)]),
        "vS": (4, [F(1, 15), F(1, 24), F(41, 2520), F(7, 1440), F(109, 90720)]),
        "pS": (5, [F(1, 36), F(1, 72), F(13, 2880), F(1, 864)]),
        "SS": (6, [F(1, 126), F(1, 288), F(13, 12960)]),
        "Sa": (3, [F(1, 12), F(1, 12), F(2, 45), F(1, 60), F(11, 2240), F(73, 60480)]),
    }
    xmax, hmax = F(1, 100), F(3, 500)
    bounds = {
        name: hmax**power * sum((c * xmax**i for i, c in enumerate(values)), F(0))
        / (2 * (1 - xmax))
        for name, (power, values) in coefficients.items()
    }
    marginal = max(
        bounds["vv"] + bounds["vp"] + bounds["va"],
        bounds["vp"] + bounds["pp"] + bounds["pa"],
        bounds["va"] + bounds["pa"] + bounds["aa"],
    )
    return marginal + bounds["vS"] + bounds["pS"] + bounds["Sa"] + bounds["SS"]


def sin_ten_degrees_interval() -> tuple[Fraction, Fraction]:
    """Rational enclosure using Machin's identity and alternating series."""
    def atan_interval(x: Fraction) -> tuple[Fraction, Fraction]:
        terms = 18
        total = sum(((-1)**k * x**(2*k + 1) / (2*k + 1)
                     for k in range(terms)), F(0))
        next_sum = total + (-1)**terms * x**(2*terms + 1) / (2*terms + 1)
        return min(total, next_sum), max(total, next_sum)

    a_lo, a_hi = atan_interval(F(1, 5))
    b_lo, b_hi = atan_interval(F(1, 239))
    x_lo, x_hi = (16*a_lo - 4*b_hi) / 18, (16*a_hi - 4*b_lo) / 18
    # Six terms end negative and underbound sin; seven end positive and overbound.
    lo = sum((((x_lo if k % 2 == 0 else x_hi)**(2*k + 1))
              * (-1)**k / factorial(2*k + 1) for k in range(6)), F(0))
    hi = sum((((x_hi if k % 2 == 0 else x_lo)**(2*k + 1))
              * (-1)**k / factorial(2*k + 1) for k in range(7)), F(0))
    return lo, hi
