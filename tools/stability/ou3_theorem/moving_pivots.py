"""Conditional six-column bridge retaining chronological transport defects.

Exact algebra subordinate to the historical-reader covariance comparison.
Supplied scalar bounds are premises; this module does not obtain them from
physical span or convert finite diagnostics into all-history certificates.
"""
from fractions import Fraction as F
from operator import index


def gyro_transport_prefixes(operations):
    """op=(h, ||R|| ceiling, ||R-I|| ceiling, ||D-hI|| ceiling).

    h=0 represents a reset and requires the fourth entry zero. Observations
    themselves are identity transitions in the latent historical array.
    h and D must be expressed in the same fixed gyro coordinate scaling.
    """
    a, b, t = F(1), F(0), F(0)
    rows = []
    for h, r, u, v in operations:
        h, r, u, v = map(F, (h, r, u, v))
        if min(h, r, u, v) < 0 or r == 0 or (h == 0 and v != 0):
            raise ValueError("invalid chronological transport bounds")
        a, b, t = r * a, r * b + u * t + v, t + h
        rows.append({"attitude_norm_upper": a, "gyro_identity_defect_upper": b,
                     "elapsed_scaled_time": t, "gyro_singular_lower": max(F(0), t - b)})
    return rows


def two_group_six_column_floor(attitude_row_floor, gyro_floor,
                                attitude_transport_upper, row_defect_upper):
    c, b, a, epsilon = map(F, (attitude_row_floor, gyro_floor,
                               attitude_transport_upper, row_defect_upper))
    if min(c, b) <= 0 or min(a, epsilon) < 0:
        raise ValueError("positive row/gyro floors and nonnegative ceilings required")
    inverse_upper = 1 + (a + 1) / b
    margin = c / inverse_upper - epsilon
    return {"block_inverse_upper": inverse_upper, "signed_margin": margin,
            "six_column_singular_lower": max(F(0), margin),
            "six_columns_certified_conditionally": margin > 0}


def six_pivot_floor_valid(singular_floor, row_count, proposed_pivot_floor):
    s, p = map(F, (singular_floor, proposed_pivot_floor))
    m = index(row_count)
    if s <= 0 or p <= 0 or m < 6:
        raise ValueError("positive floors and at least six actual rows required")
    return m * p * p <= s * s


def certificate():
    from .matrix_certificates import add, identity, ldlt, matmul, transpose
    # Supplied normalized exact example, NOT an observed shipping word.
    # C0=C1=I, A=I, B=I/2, E diagonal 1/100. Include a reset in the
    # separate scalar recurrence to audit chronological defect accounting.
    ops = [(F(1, 4), F(1), F(0), F(0)),
           (F(0), F(11, 10), F(1, 10), F(0)),
           (F(1, 4), F(1), F(0), F(0))]
    transport = gyro_transport_prefixes(ops)[-1]
    bound = two_group_six_column_floor(1, F(1, 2), 1, F(1, 100))
    s = bound["six_column_singular_lower"]
    o = [[F(0) for _ in range(6)] for _ in range(6)]
    for i in range(3):
        o[i][i] = F(101, 100)
        o[i + 3][i] = F(1)
        o[i + 3][i + 3] = F(51, 100)
    gram = matmul(transpose(o), o)
    _, pivots = ldlt(add(gram, identity(6), -s*s))
    assert all(x > 0 for x in pivots)
    p = s / 3  # 6/9 <= 1, no floating square root in the certificate.
    assert six_pivot_floor_valid(s, 6, p)
    failed = two_group_six_column_floor(1, F(1, 2), 1, F(1, 4))
    assert not failed["six_columns_certified_conditionally"]
    # Exact quiet nominal subcase: C=[skew(g ez);skew(B ex)], with
    # g>=9, B>=9; two coincident acc/mag groups eight qualified predictions apart, zero
    # corrections/resets. C'C >=81 I, so this is an analytic bound,
    # not an eigenvalue sampled from a Gramian or a physical bias theorem.
    quiet = two_group_six_column_floor(9, F(4, 125), 1, 0)
    quiet_s = quiet["six_column_singular_lower"]
    quiet_p = quiet_s / 4
    assert quiet_s == F(18, 127) and six_pivot_floor_valid(quiet_s, 12, quiet_p)
    return {
        "qualification": "OU3_CHRONOLOGICAL_SIX_PIVOT_IMPLICATION_V1",
        "quiet_nominal_two_group_singular_floor": str(quiet_s),
        "quiet_nominal_all_six_pivot_floor": str(quiet_p),
        "quiet_nominal_scope": "zero-residual regular nominal source history, g>=9 and B>=9, actually applied acc/mag groups eight qualified predictions apart, raw AG scaling",
        "quiet_nominal_floor_identifies_physical_attitude_and_bias": False,
        "conditional_chronological_gyro_floor": True,
        "conditional_two_group_six_column_floor": True,
        "singular_floor_to_all_six_greedy_pivots": True,
        "supplied_example_singular_floor": str(s),
        "supplied_example_all_six_pivot_floor": str(p),
        "supplied_example_exact_gram_residual_pivots": [str(x) for x in pivots],
        "supplied_reset_sequence_gyro_floor": str(transport["gyro_singular_lower"]),
        "supplied_failed_defect_margin": str(failed["signed_margin"]),
        "actual_chronological_resets_may_be_dropped": False,
        "physical_span_supplies_uniform_row_and_defect_budget": False,
        "two_temporal_margins_imply_six_pivots": False,
        "uniform_historical_AG_readout_action": False,
        "rho0_certified": False,
    }
