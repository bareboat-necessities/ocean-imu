"""Exact rational certificate for the literal 0x5f375a86 / one-Newton step.

No trajectories or mantissa enumeration. On each of three exponent strata,
s*y_seed(s)^2 is a cubic with one possible interior maximum. Its extrema are
rational. IEEE round-to-nearest and gradual underflow are explicit premises.
This closes a LOCAL norm bound, not the Mahony pitch/integral invariant tube.
"""
from __future__ import annotations
from fractions import Fraction as F

MAGIC = 0x5F375A86
M = 2**23
U = F(1, 2**24)


def cubic_range(lo, hi, a, b):
    """Exact range of s*(a-b*s)^2 on a positive-seed interval."""
    assert 0 < lo <= hi and a - b*hi > 0
    points = [lo, hi]
    critical = a / (3*b)
    if lo <= critical <= hi:
        points.append(critical)
    values = [s*(a-b*s)**2 for s in points]
    return min(values), max(values)


def seed_strata():
    # Reducing the exponent by an even integer is exact power-of-four scaling.
    # floor(m/2) = m/2 - epsilon, epsilon in {0,1/2}; bound both parities.
    d = MAGIC - 190*M
    return (
        (F(1), F(2)-F(1, M), F(MAGIC, 2*M)-94, F(1, 4), F(1, 4*M)),
        (F(2), 2*(1+F(2*d+1, M)), F(MAGIC, 2*M)-F(377, 4), F(1, 8), F(1, 4*M)),
        (2*(1+F(2*d+2, M)), F(4)-F(2, M), F(MAGIC, 4*M)-F(375, 8), F(1, 16), F(1, 8*M)),
    )


def certificate():
    zlo, zhi = F(93, 100), F(108, 100)
    strata = []
    for lo, hi, a, b, eps in seed_strata():
        lower = cubic_range(lo, hi, a, b)[0]
        upper = cubic_range(lo, hi, a+eps, b)[1]
        assert zlo <= lower <= upper <= zhi
        strata.append({"input_interval": [str(lo), str(hi)],
                       "seed": {"a": str(a), "b": str(b), "epsilon_max": str(eps)},
                       "s_seed_squared_bounds": [str(lower), str(upper)]})
    # x*0.5 exact in the stated input domain. The two products before the
    # subtraction supply t; the subtraction and final multiplication supply v.
    # x*y_out^2 = z*(3-t*z)^2/4 * v, with t=(1+d1)(1+d2),
    # v=(1+d3)^2(1+d4)^2. Keep the repeated z linked.
    tlo, thi = (1-U)**2, (1+U)**2
    lower = min(z*(3-thi*z)**2/4 for z in (zlo, zhi))*(1-U)**4
    # For t fixed the only interior extremum is a maximum at z=1/t.
    upper = (1+U)**4/tlo
    # Rounded dot product: 4 products, 3 additions, path length <=4.
    # Standard relative gamma_4 plus <=8*min_subnormal absolute covers every
    # gradual-underflow rounding. pre-norm^2>=1/4 gives <=32*tiny relative;
    # 64*tiny is a rational conservative bound. Final component products give
    # squared relative (1+-u)^2 plus <=16*tiny (the norm is <2).
    tiny = F(1, 2**149)
    gamma = 4*U/(1-4*U)
    dot_rel = gamma + 64*tiny
    raw_lower = lower/(1+dot_rel)*(1-U)**2 - 16*tiny
    raw_upper = upper/(1-dot_rel)*(1+U)**2 + 16*tiny
    assert raw_lower > F(995, 1000) and raw_upper < F(1001, 1000)
    # The same exponent reduction and dot-product bound apply through 128.
    # x/2 and all seed/Newton intermediates remain normal. This includes the
    # planar accelerometer (norm squared near 96), not only quaternions.
    scalar_lower, scalar_upper = lower/(1+dot_rel), upper/(1-dot_rel)
    assert scalar_lower > F(995, 1000) and scalar_upper < F(1001, 1000)
    return {
        "qualification": "OU3_MAHONY_RAW_NORMALIZATION_RATIONAL_V1",
        "result_type": "PROVED — analytical",
        "method": "exact rational extrema of three cubic strata; no sampled mantissas",
        "source": "src/ahrs/Mahony_AHRS.h: invSqrt<float> and update normalization",
        "magic_hex": hex(MAGIC),
        "premises": ["IEEE binary32, round-to-nearest, gradual underflow, no unsafe reassociation",
                     "represented pre-normalization quaternion has exact squared norm in [1/4,4]",
                     "literal bit seed followed by exactly one Newton step"],
        "seed_strata": strata,
        "linked_seed_squared_enclosure": [str(zlo), str(zhi)],
        "newton_x_y_squared_bounds": [str(lower), str(upper)],
        "raw_output_squared_norm_bounds": [str(raw_lower), str(raw_upper)],
        "certified_coarse_raw_squared_norm_interval": ["199/200", "1001/1000"],
        "extended_input_squared_norm_domain": ["1/4", "128"],
        "extended_domain_argument": "even exponent shifts leave s*y_seed^2 unchanged; all intermediate products remain normal; identical dot/rounding bound",
        "pre_component_scalar_squared_norm_bounds": [str(scalar_lower), str(scalar_upper)],
        "structures_preserved": ["literal magic constant", "one Newton step", "raw quaternion", "rounding charges"],
        "relaxations_introduced": ["seed parity error enclosed by its interval", "rounding error outer bounds"],
        "all_time_pre_normalization_domain_verified": False,
        "pitch_integral_tube_verified": False,
        "all_time_magnetic_service_verified": False,
        "theorem_closed": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
