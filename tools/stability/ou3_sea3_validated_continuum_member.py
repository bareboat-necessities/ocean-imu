#!/usr/bin/env python3
"""Validated analytic member of the complete SEA3 continuum hard-driver ball.

This is source machinery for the mandatory complete-word feasibility experiment,
not a second source family or a P4 certificate.  A single continuum coefficient
field is supported on the continuous band [0.295, 0.305] Hz and is carried
without reseeding through startup and the 3 s word.  It cancels the admitted
JONSWAP/vertical-acceleration spectral factor on that band, so the physical
vertical acceleration is the exact continuum integral

  a_z(t)=C int[f1,f2] cos(2*pi*f*(t-tc)) df.

The integral is evaluated with outward Decimal arithmetic and a validated sine
Taylor enclosure.  No quadrature nodes, finite harmonics, replay samples, or
independent per-sample amplitudes occur anywhere in the construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
import math

PREC = 90
PI_LO = Decimal("3.14159265358979323846264338327950288419716939937510")
PI_HI = Decimal("3.14159265358979323846264338327950288419716939937511")
F1 = Decimal("0.295")
F2 = Decimal("0.305")
C = Decimal("6")
PHASE_CENTER_S = Decimal("30")
DT = Decimal("0.005")
PREHISTORY_S = Decimal("60")
WINDOW_SAMPLES = 601


def _lo(fn) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PREC
        ctx.rounding = ROUND_FLOOR
        return +fn()


def _hi(fn) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PREC
        ctx.rounding = ROUND_CEILING
        return +fn()


@dataclass(frozen=True)
class DInterval:
    lo: Decimal
    hi: Decimal

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed Decimal interval")

    @staticmethod
    def point(x: Decimal | int | str) -> "DInterval":
        d = Decimal(x)
        return DInterval(d, d)

    def __add__(self, other: "DInterval") -> "DInterval":
        return DInterval(_lo(lambda: self.lo + other.lo), _hi(lambda: self.hi + other.hi))

    def __sub__(self, other: "DInterval") -> "DInterval":
        return DInterval(_lo(lambda: self.lo - other.hi), _hi(lambda: self.hi - other.lo))

    def __neg__(self) -> "DInterval":
        return DInterval(-self.hi, -self.lo)

    def __mul__(self, other: "DInterval") -> "DInterval":
        pairs = [(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        lows = [_lo(lambda a=a, b=b: a * b) for a, b in pairs]
        highs = [_hi(lambda a=a, b=b: a * b) for a, b in pairs]
        return DInterval(min(lows), max(highs))

    def reciprocal(self) -> "DInterval":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        lows = [_lo(lambda x=x: Decimal(1) / x) for x in (self.lo, self.hi)]
        highs = [_hi(lambda x=x: Decimal(1) / x) for x in (self.lo, self.hi)]
        return DInterval(min(lows), max(highs))

    def __truediv__(self, other: "DInterval") -> "DInterval":
        return self * other.reciprocal()

    def float_bounds(self) -> list[float]:
        # Decimal->binary64 is nearest; widen one ulp on each side.
        return [
            math.nextafter(float(self.lo), -math.inf),
            math.nextafter(float(self.hi), math.inf),
        ]


def _pi() -> DInterval:
    return DInterval(PI_LO, PI_HI)


def _mul_pi(value: Decimal, factor: int = 1) -> DInterval:
    vals = [Decimal(factor) * PI_LO * value, Decimal(factor) * PI_HI * value]
    return DInterval(_lo(lambda: min(vals)), _hi(lambda: max(vals)))


def _sin_small(x: DInterval) -> DInterval:
    half_pi_hi = PI_HI / Decimal(2)
    if x.lo < -half_pi_hi or x.hi > half_pi_hi:
        raise ValueError("Taylor sine input outside [-pi/2,pi/2]")
    total = DInterval.point(0)
    x2 = x * x
    power = x
    fact = Decimal(1)
    sign = 1
    # Terms x^(2n+1)/(2n+1)! through degree 61.
    for n in range(31):
        if n:
            fact *= Decimal(2 * n) * Decimal(2 * n + 1)
            power = power * x2
        term = power / DInterval.point(fact)
        total = total + (term if sign > 0 else -term)
        sign *= -1
    with localcontext() as ctx:
        ctx.prec = PREC
        ctx.rounding = ROUND_CEILING
        rem = +(half_pi_hi ** 63) / Decimal(math.factorial(63))
    return DInterval(_lo(lambda: total.lo - rem), _hi(lambda: total.hi + rem))


def sin_interval(x: DInterval) -> DInterval:
    """Outward sine enclosure for the narrow argument boxes used by this source."""
    two_pi = DInterval(_lo(lambda: Decimal(2) * PI_LO), _hi(lambda: Decimal(2) * PI_HI))
    pi_mid = (PI_LO + PI_HI) / Decimal(2)
    x_mid = (x.lo + x.hi) / Decimal(2)
    n = int((x_mid / (Decimal(2) * pi_mid)).to_integral_value(rounding=ROUND_HALF_EVEN))
    r = x - (DInterval.point(n) * two_pi)
    # Verify/recover one neighboring period if midpoint selection was on a boundary.
    if r.lo > PI_HI:
        r = r - two_pi
    elif r.hi < -PI_HI:
        r = r + two_pi
    if r.lo < -PI_HI or r.hi > PI_HI:
        raise RuntimeError("sine range reduction failed")

    half_lo = PI_LO / Decimal(2)
    half_hi = PI_HI / Decimal(2)
    if r.lo >= half_hi:
        r = _pi() - r
    elif r.hi <= -half_hi:
        r = -(_pi() + r)
    elif r.lo < -half_lo or r.hi > half_lo:
        # Only possible when the tiny pi enclosure straddles a sine extremum.
        return DInterval(Decimal(-1), Decimal(1))
    return _sin_small(r)


def acceleration_interval(t_s: Decimal | str | float) -> DInterval:
    """Outward enclosure of the same continuum member at absolute source time t."""
    t = Decimal(str(t_s)) if not isinstance(t_s, Decimal) else t_s
    tau = t - PHASE_CENTER_S
    if tau == 0:
        return DInterval.point(C * (F2 - F1))
    arg2 = _mul_pi(F2 * tau, 2)
    arg1 = _mul_pi(F1 * tau, 2)
    numerator = sin_interval(arg2) - sin_interval(arg1)
    denominator = _mul_pi(tau, 2)
    return DInterval.point(C) * numerator / denominator


def driver_norm_certificate() -> dict:
    """Analytic lower-kernel proof that the coefficient field is in the unit ball."""
    H = Decimal("1.5")
    Tp = Decimal("6")
    gamma = Decimal("3.3")
    fp = Decimal(1) / Tp
    m0 = H * H / Decimal(16)
    # raw JONSWAP integral <= gamma/(5 fp^4), hence normalized scale >= below.
    scale_lo = m0 * Decimal(5) * fp**4 / gamma
    x_hi = Decimal(5) / Decimal(4) * (fp / F1) ** 4
    exp_lo = Decimal(1) - x_hi  # exp(-x) >= 1-x
    if exp_lo <= 0:
        raise RuntimeError("driver support does not admit positive exponential lower bound")
    shape_lo = F2**Decimal(-5) * exp_lo  # gamma^peak >= 1
    spectrum_lo = scale_lo * shape_lo
    two_pi_f_lo = Decimal(2) * PI_LO * F1
    kernel2_lo = spectrum_lo * two_pi_f_lo**4  # vertical RAO is exactly one on band
    norm2_hi = C * C * (F2 - F1) / kernel2_lo
    return {
        "normalization_scale_lower": str(scale_lo),
        "jonswap_exponential_lower": str(exp_lo),
        "spectrum_density_lower_on_driver_band": str(spectrum_lo),
        "acceleration_kernel_squared_lower_on_driver_band": str(kernel2_lo),
        "driver_norm_squared_upper": str(norm2_hi),
        "driver_norm_strictly_below_one": norm2_hi < Decimal(1),
    }


def self_check() -> dict:
    cert = driver_norm_certificate()
    if not cert["driver_norm_strictly_below_one"]:
        raise RuntimeError("continuum member escaped complete-SEA3 hard driver ball")
    max_width = Decimal(0)
    max_abs = Decimal(0)
    # Validate the complete 60 s prehistory and the 601-sample word using the
    # same analytic formula, without storing 12,601 samples here.
    total = int(PREHISTORY_S / DT) + WINDOW_SAMPLES
    for k in range(total):
        v = acceleration_interval(Decimal(k) * DT)
        max_width = max(max_width, v.hi - v.lo)
        max_abs = max(max_abs, abs(v.lo), abs(v.hi))
    return {
        "driver_norm_certificate": cert,
        "prehistory_samples_checked": int(PREHISTORY_S / DT),
        "word_samples_checked": WINDOW_SAMPLES,
        "max_interval_width": str(max_width),
        "max_abs_acceleration_mps2": str(max_abs),
        "quadrature_used": False,
        "finite_harmonic_source_used": False,
        "independent_sample_boxes_used": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2, sort_keys=True))
