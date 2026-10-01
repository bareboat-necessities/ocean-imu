"""Exact same-word dissipation/supply algebra; never a shipping certificate.

All operands must come from ONE frozen literal history. b is the signed,
transported endpoint defect, not a freely chosen innovation control. See
 docs/ou3-linked-finite-supply.md. Algebra verification does not establish
source-uniform loss, reachable defect bounds, release, or prefix retention.
"""
from __future__ import annotations

from fractions import Fraction as F
from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, ldlt, matmul, transpose


def matrix(value):
    """Convert a nonempty rectangular matrix to exact rational entries."""
    result = [[F(x) for x in row] for row in value]
    if not result or not result[0] or any(len(r) != len(result[0]) for r in result):
        raise ValueError("nonempty rectangular matrix required")
    return result


def _spd(value):
    result = matrix(value)
    if len(result) != len(result[0]) or result != transpose(result):
        raise ValueError("symmetric square matrix required")
    ldlt(result)
    return result


def quad(value, precision):
    return matmul(matmul(transpose(value), precision), value)[0][0]


def compose_defects(events, dimension):
    """Return exact M,b for e_next=A_i e+d_i in chronological order.

    No word-root reset or source independence is assumed. d_i includes
    nonlinear, physical, projection and arithmetic defects as applicable.
    """
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
        raise ValueError("positive integer dimension required")
    m = identity(dimension)
    b = [[F(0)] for _ in range(dimension)]
    for event in events:
        a, d = matrix(event['A']), matrix(event['d'])
        if len(a) != dimension or len(a[0]) != dimension or len(d) != dimension or len(d[0]) != 1:
            raise ValueError("operation dimension mismatch")
        m, b = matmul(a, m), add(matmul(a, b), d)
    return m, b


def linked_supply(root_precision, end_precision, transport, endpoint_defect, gamma):
    """Complete the square without separating forcing from loss directions.

    G=J0-M'JN M-gamma J0 must be SPD. With z=M'JN b,
    chi=b'JN b+z'G^-1 z is the SMALLEST fixed-b constant in
      (M e+b)'JN(M e+b) <= (1-gamma)e'J0 e+chi
    for unrestricted e. A physical theorem still needs a bound on chi over
    the actual linked reachable pairs (M,b); this routine proves no such bound.
    """
    j0, jn = _spd(root_precision), _spd(end_precision)
    m, b, gamma = matrix(transport), matrix(endpoint_defect), F(gamma)
    if not 0 < gamma < 1:
        raise ValueError("gamma must lie strictly between zero and one")
    if len(m) != len(jn) or len(m[0]) != len(j0) or len(b) != len(jn) or len(b[0]) != 1:
        raise ValueError("word dimension mismatch")
    loss = add(j0, matmul(matmul(transpose(m), jn), m), F(-1))
    g = add(loss, j0, -gamma)
    _spd(g)  # Fail closed: never replace an indefinite/singular G by clipping.
    z = matmul(matmul(transpose(m), jn), b)
    center = matmul(inverse(g), z)
    chi = quad(b, jn) + matmul(transpose(z), center)[0][0]
    return {'loss': loss, 'G': g, 'center': center, 'chi': chi, 'gamma': gamma,
            'algebra_verified': True, 'source_uniform_verified': False,
            'theorem_closed': False}


def completed_square_check(root_precision, end_precision, transport, endpoint_defect,
                           gamma, root_error):
    cert = linked_supply(root_precision, end_precision, transport, endpoint_defect, gamma)
    e, b, m = matrix(root_error), matrix(endpoint_defect), matrix(transport)
    if len(e) != len(root_precision) or len(e[0]) != 1:
        raise ValueError("root-error dimension mismatch")
    terminal = add(matmul(m, e), b)
    displaced = add(e, cert['center'], F(-1))
    lhs = quad(terminal, matrix(end_precision))
    rhs = ((1-cert['gamma'])*quad(e, matrix(root_precision)) + cert['chi']
           - quad(displaced, cert['G']))
    if lhs != rhs:
        raise ArithmeticError("linked completed-square identity failed")
    return {**cert, 'identity_verified': True, 'terminal_storage': lhs}


def retained_entry_budget(gamma, supply_upper, radius_squared):
    """Budget arithmetic only; inputs must already be uniform proved bounds.

    Equality at the boundary proves a root self-map, not finite entry into
    that boundary: x_next=(1-gamma)x+gamma*C can approach C from above forever.
    Every-prefix bounds are a separate obligation.
    """
    gamma, supply, radius = map(F, (gamma, supply_upper, radius_squared))
    if not 0 < gamma < 1 or supply < 0 or radius <= 0:
        raise ValueError("positive radius, nonnegative supply, gamma in (0,1) required")
    return {'root_self_map_budget': supply <= gamma*radius,
            'strict_entry_budget': supply < gamma*radius,
            'limiting_storage_bound': supply/gamma,
            'prefix_retention_verified': False, 'source_uniform_verified': False}


def imu_supply_outer(operators, steps_s, limits):
    """Conditional SF5 sensor support on ONE already-frozen literal word.

    operators contains accel_slow, accel_fast, gyro_slow, gyro_fast matrices
    (one output-by-3 coefficient per complete delivered cell). Include actual
    prediction/correction transport and bias-model mismatch in these matrices.
    A finite outer support is NOT the linked chi supremum: loss and forcing
    still share the same physical/covariance/tuner history. Do not use this
    norm relaxation to revive the failed independent-extrema proof route.
    """
    from .imu_temporal import (OpenTemporalQualification, slow_weighted_outer,
                               fast_weighted_outer)
    import numpy as np
    keys = {'accel_slow', 'accel_fast', 'gyro_slow', 'gyro_fast'}
    if set(operators) != keys:
        raise ValueError('both slow and fast transported operators required for BOTH sensors')
    if not limits.temporal_parameters_present:
        raise OpenTemporalQualification('both delivered-stream temporal profiles remain required')
    dt = np.asarray(steps_s, dtype=float)
    if dt.ndim != 1 or len(dt) == 0 or not np.isfinite(dt).all() or np.any(dt <= 0):
        raise ValueError('positive finite complete-cell durations required')
    shapes = {np.asarray(operators[k]).shape for k in keys}
    if (len(shapes) != 1 or len(next(iter(shapes))) != 3
            or next(iter(shapes))[0] != len(dt) or next(iter(shapes))[2] != 3):
        raise ValueError('all four operators must share one word and output coordinates')
    terms = {}
    for name,bs,ds,bf,fw in (
            ('accel',limits.B_a_s_mps2,limits.D_a_s_mps3,limits.B_a_f_mps2,limits.accel_fast_window),
            ('gyro',limits.B_g_s_rad_s,limits.D_g_s_rad_s2,limits.B_g_f_rad_s,limits.gyro_fast_window)):
        terms[name+'_slow'] = slow_weighted_outer(operators[name+'_slow'],dt[:-1],bs,ds)
        terms[name+'_fast'] = fast_weighted_outer(operators[name+'_fast'],dt,bf,fw)
    return {'terms':terms,'conditional_support_outer':sum(terms.values()),
            'source_uniform_verified':False,'linked_chi_supremum_verified':False,
            'physical_device_qualified':False,'decomposition_reselected':False}
