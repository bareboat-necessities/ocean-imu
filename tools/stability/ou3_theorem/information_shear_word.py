"""Exact information-shear word action, without a uniform-gap assertion.

Analytical identities are proved in app:information-shear-word. Rational
operands test identities, not reachable cells. Prefix maps must be derivatives
of one complete causal shipping history, including auxiliary-state feedback.
"""
from __future__ import annotations
from fractions import Fraction as F

from .lin_path_certificate import inverse
from .matrix_certificates import add, is_psd, ldlt, transpose
from .planar_complete_word_storage import product


def zeros(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def trace(a):
    return sum((row[i] for i, row in enumerate(a)), F(0))


def symmetric_basis(n):
    """Unscaled symmetric coordinates; metric retains both off-diagonal terms."""
    result = []
    for i in range(n):
        for j in range(i, n):
            a = zeros(n, n)
            a[i][j] = a[j][i] = F(1)
            result.append(a)
    return result


def joint_metric(P, weight):
    """eta' P^-1 eta + weight tr(P^-1 dP P^-1 dP)."""
    ldlt(P)
    weight = F(weight)
    if weight <= 0:
        raise ValueError('positive covariance weight required')
    n, J = len(P), inverse(P)
    basis = symmetric_basis(n)
    out = zeros(n+len(basis), n+len(basis))
    for i in range(n):
        out[i][:n] = J[i]
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            out[n+i][n+j] = weight*trace(product(J, a, J, b))
    return out


def joint_congruence(B):
    """Base tangent map (eta,dP) -> (B eta,B dP B')."""
    n = len(B)
    if any(len(row) != n for row in B):
        raise ValueError('square base map required')
    basis = symmetric_basis(n)
    out = zeros(n+len(basis), n+len(basis))
    for i in range(n):
        out[i][:n] = B[i]
    indices = [(i, j) for i in range(n) for j in range(i, n)]
    for col, a in enumerate(basis):
        image = product(B, a, transpose(B))
        for row, (i, j) in enumerate(indices):
            out[n+row][n+col] = image[i][j]
    return out


def linked_action_balance(base_maps, metrics, prefix_maps):
    """Exact root quadratic identity with state-dependent ports retained.

    v_i=T_i a, p_i=(T_(i+1)-L_i T_i)a. The ONE root a includes inherited
    auxiliary/source coordinates or their justified lift. T_i may be
    rectangular. No Markov closure on mean/P alone is assumed. Zero-action
    kernel is intersection ker(Q_i T_i), Q_i=J_i-L_i'J_(i+1)L_i >= 0.
    Signed port work is NOT an exogenous forcing bound.
    """
    if not base_maps or len(metrics) != len(base_maps)+1 or len(prefix_maps) != len(metrics):
        raise ValueError('nonempty paired maps and boundary metrics required')
    for metric in metrics:
        ldlt(metric)
    r = len(prefix_maps[0][0])
    action, work = zeros(r, r), zeros(r, r)
    local_losses, ports = [], []
    for i, L in enumerate(base_maps):
        T, Tnext, Jnext = prefix_maps[i], prefix_maps[i+1], metrics[i+1]
        Q = add(metrics[i], product(transpose(L), Jnext, L), -1)
        if not is_psd(Q):
            raise ValueError('D_SUFFICIENT_BOUND_FAILURE: selected base action is not PSD')
        LT = product(L, T)
        E = add(Tnext, LT, -1)
        cross = product(transpose(LT), Jnext, E)
        action = add(action, product(transpose(T), Q, T))
        work = add(work, add(add(cross, transpose(cross)), product(transpose(E), Jnext, E)))
        local_losses.append(Q)
        ports.append(E)
    initial = product(transpose(prefix_maps[0]), metrics[0], prefix_maps[0])
    final = product(transpose(prefix_maps[-1]), metrics[-1], prefix_maps[-1])
    gap = add(action, work, -1)
    if gap != add(initial, final, -1):
        raise ArithmeticError('linked word identity failed')
    return {'action': action, 'signed_port_work': work, 'gap': gap,
            'initial': initial, 'terminal': final, 'local_losses': local_losses,
            'root_port_maps': ports}


def precision_projector(J, basis):
    """Fixed-rank basis only. Rank loss fails at the actual Gram inverse."""
    ldlt(J)
    gram = product(transpose(basis), J, basis)
    ldlt(gram)
    return product(basis, inverse(gram), transpose(basis), J)


def projector_differential(J, basis, dJ, dbasis):
    """Exact dPi, retaining precision and physical-basis transport."""
    gram = product(transpose(basis), J, basis)
    ldlt(gram)
    U = inverse(gram)
    dgram = add(add(product(transpose(dbasis), J, basis),
                     product(transpose(basis), dJ, basis)),
                 product(transpose(basis), J, dbasis))
    out = add(product(dbasis, U, transpose(basis), J),
              product(basis, U, transpose(dbasis), J))
    out = add(out, product(basis, U, transpose(basis), dJ))
    return add(out, product(basis, U, dgram, U, transpose(basis), J), -1)


def congruence_shear(P, G, e, mismatch, dP, dG, de, dmismatch):
    """P+=GPG', e+=Ge+t; t is the LITERAL chart/injection discrepancy.

    G is shipping covariance reset, not an asserted exact finite SO(3) map.
    No bound on t or dt is inferred. G must be invertible on this branch.
    """
    C = product(G, P, transpose(G))
    ldlt(P)
    ldlt(C)
    dC = add(add(product(G, dP, transpose(G)), product(dG, P, transpose(G))),
             product(G, P, transpose(dG)))
    ep = add(product(G, e), mismatch)
    dep = add(add(product(G, de), product(dG, e)), dmismatch)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(dep, product(dC, inverse(C), ep), -1)
    linked = add(product(G, eta), product(G, P, transpose(dG), inverse(transpose(G)), inverse(P), e), -1)
    linked = add(add(linked, dmismatch), product(dC, inverse(C), mismatch), -1)
    assert direct == linked
    return {'eta_plus': direct, 'covariance_plus': C, 'covariance_tangent_plus': dC,
            'base': product(G, eta), 'port': add(linked, product(G, eta), -1)}


def covariance_increment_shear(P, increment, e, dP, dincrement, de):
    """Default AW: actual PSD increment, including its face/target derivative."""
    if not is_psd(increment):
        raise ValueError('actual PSD increment required')
    C = add(P, increment)
    ldlt(P)
    ldlt(C)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(de, product(add(dP, dincrement), inverse(C), e), -1)
    linked = add(eta, product(dP, add(inverse(P), inverse(C), -1), e))
    linked = add(linked, product(dincrement, inverse(C), e), -1)
    assert linked == direct
    return {'eta_plus': direct, 'base': eta, 'port': add(linked, eta, -1)}


def mean_shift_shear(P, e, shift, dP, de, dshift):
    """Bias projection with unchanged P: no covariance-metric contraction assumed."""
    ldlt(P)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(add(de, dshift), product(dP, inverse(P), add(e, shift)), -1)
    linked = add(add(eta, dshift), product(dP, inverse(P), shift), -1)
    assert direct == linked
    return {'eta_plus': direct, 'port': add(linked, eta, -1)}


def certificate():
    return {
        'qualification': 'OU3_INFORMATION_SHEAR_WORD_V1',
        'result_type': 'PROVED — analytical balance and kernel characterization; uniform coercivity OPEN',
        'proof_appendix_label': 'app:information-shear-word',
        'governing_strategy': 'symmetry/kernel -> justified quotient/tube -> linked dissipation -> forcing -> minimal activation domain',
        'general_theorem_assumes_existing_magnetic_service': True,
        'general_theorem_waits_for_planar_admission': False,
        'storage': 'W=eta^T J eta+lambda tr(J dP J dP), eta=de-dP J e, lambda>0',
        'word_balance': 'W_N-W_0=-sum v_i^T Q_i v_i+sum[2(L_i v_i)^T J_(i+1) p_i+p_i^T J_(i+1) p_i]',
        'port_scope': 'literal linked state-dependent/source ports, not independent disturbances; inherited auxiliary tangents remain causal',
        'covariance_zero_loss': '(I-U^T U) Z=0, Z=P^-1/2 dP P^-1/2',
        'optimal_correction_zero_loss': 'H eta=0 and H dP=0; reached held BA uses active block and retained held ports',
        'word_action_kernel': 'intersection_i ker(Q_i T_i) for ACTUAL complete causal prefix tangents T_i',
        'quotient_balance': 'Wperp_N-Wperp_0=-A_W+F_W+gauge_energy_0-gauge_energy_N',
        'moving_projector_differential_retained': True,
        'physical_compatibility_equals_homogeneous_kernel': False,
        'source_uniform_coercivity_requirement': 'B^T[G_action-G_port-G_gauge,0+G_gauge,N]B >= c B^T Jperp,0 B, c>0, on a justified transverse lift and every admitted activated word',
        'relative_Schur_obligation': 'for G-cJ=[[A,B],[B^T,N]], N>0 and A-B N^-1 B^T>=0 on justified service-coordinate strata; actual cross precision retained',
        'uniform_c': None,
        'minimal_obligations_note': 'docs/ou3-information-shear-strategy.md',
        'full_literal_derivative_domain_certified': False,
        'nonlinear_remainder_and_precision_transfer_certified': False,
        'structures_preserved': ['complete 21-state mean and covariance', 'literal Joseph/masks and effective held noise',
            'prediction mismatch and dF/dQ', 'dH/dR/dr', 'injection/reset and bias projection',
            'default AW pending-target/face differential', 'inherited reference/gates/clocks',
            'same-history gauge/source and moving projectors'],
        'relaxations_introduced': ['regular real-operation tangent identities; nonsmooth/arithmetic closure remains separate',
            'dissipative base plus exact signed ports is an identity, not a surrogate execution'],
        'radius_solve_performed': False,
        'all_time_planar_magnetic_service_verified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
