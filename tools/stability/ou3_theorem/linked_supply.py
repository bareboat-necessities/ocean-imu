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
