"""Linked complete-word storage algebra, not a claimed contraction factor.

All columns in an input forcing matrix refer to ONE persistent root/history
coordinate. They are added before norms. Actual joint quotient/covariance
Jacobians, moving charts, reset/AW faces and gauge forcing are required inputs.
No frozen-gain factor word or arbitrary SPD covariance certifies the caller.
"""
from __future__ import annotations

from fractions import Fraction as F

from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, matmul, transpose


def product(*values):
    out = values[0]
    for value in values[1:]:
        out = matmul(out, value)
    return out


def compose_linked_word(factors, forcing_maps):
    """v_(i+1)=A_i v_i+B_i w, with the SAME w at every operation.

    w includes the persistent gauge/root/source coordinates, after their
    causal propagation has been included in B_i; it is not fresh noise.
    """
    if not factors or len(factors) != len(forcing_maps):
        raise ValueError('nonempty paired operation and forcing maps required')
    n, m = len(factors[0]), len(forcing_maps[0][0])
    M, N = identity(n), [[F(0) for _ in range(m)] for _ in range(n)]
    for A, B in zip(factors, forcing_maps):
        M, N = product(A, M), add(product(A, N), B)
    return M, N


def storage_identity(factors, forcing_maps, root_metric, terminal_metric):
    """W_N-W_0=-v'Qv+2v'Xw+w'Zw, an exact complete-word identity."""
    ldlt(root_metric)
    ldlt(terminal_metric)
    M, N = compose_linked_word(factors, forcing_maps)
    return {'M': M, 'N': N,
            'Q': add(root_metric, product(transpose(M), terminal_metric, M), -1),
            'X': product(transpose(M), terminal_metric, N),
            'Z': product(transpose(N), terminal_metric, N)}


def telescoped_gap(factors, metrics):
    """Sum T_i' (J_i-A_i'J_(i+1)A_i) T_i, preserving all cross terms."""
    if len(metrics) != len(factors)+1 or not factors:
        raise ValueError('one metric per operation boundary required')
    n = len(factors[0])
    T = identity(n)
    gap = [[F(0) for _ in range(n)] for _ in range(n)]
    for i, A in enumerate(factors):
        local = add(metrics[i], product(transpose(A), metrics[i+1], A), -1)
        gap = add(gap, product(transpose(T), local, T))
        T = product(A, T)
    return gap


def storage_inequality(identity_record, root_metric, epsilon):
    """Exact matrix sufficient inequality, NOT source-domain certification.

    The caller must prove the same inequality over the entire inherited
    causal domain. Evaluating it at a point word does not do that.
    """
    epsilon = F(epsilon)
    if not 0 < epsilon <= 1:
        raise ValueError('strict positive storage margin <=1 required')
    if not is_psd(add(identity_record['Q'], root_metric, -epsilon)):
        raise ValueError('D_SUFFICIENT_BOUND_FAILURE: complete-word gap lacks the claimed margin')
    X, Z = identity_record['X'], identity_record['Z']
    supply = add(Z, product(transpose(X), inverse(root_metric), X), 2/epsilon)
    return {'rho': 1-epsilon/2, 'supply_matrix': supply}


def optimal_correction_storage(P, H, R, e, gain_residual_force):
    """Exact regular unmasked correction identity, with dK*r retained.

    Applies ONLY where gain and innovation use the same unmasked H. On the
    held-BA branch use its literal covariance and the general word identity.
    e and force are columns; force=dK*r plus any separately identified chart
    terms. The identity never bounds the force by an independent point max.
    """
    S = add(product(H, P, transpose(H)), R)
    K = product(P, transpose(H), inverse(S))
    C = add(P, product(K, H, P), -1)
    J, Jc = inverse(P), inverse(C)
    A = add(identity(len(P)), product(K, H), -1)
    out = add(product(A, e), gain_residual_force)
    loss = product(transpose(e), transpose(H), inverse(S), H, e)[0][0]
    cross = 2*product(transpose(e), J, gain_residual_force)[0][0]
    forcing = product(transpose(gain_residual_force), Jc, gain_residual_force)[0][0]
    initial = product(transpose(e), J, e)[0][0]
    final = product(transpose(out), Jc, out)[0][0]
    assert final == initial-loss+cross+forcing
    return {'initial': initial, 'final': final, 'loss': loss,
            'cross': cross, 'forcing': forcing}


def covariance_residual_charge(P, H, R, dP, residual):
    """Exact evaluation of the proved f_P energy <= (NIS/2)*L_P identity.

    L_P is the SAME operation's homogeneous covariance Fisher-storage loss.
    The universal proof uses D+E=I, with D,E PSD, in the appendix. This
    evaluation is useful for algebra regression, not future-domain sampling.
    """
    ldlt(P)
    ldlt(R)
    if dP != transpose(dP):
        raise ValueError('symmetric covariance tangent required')
    S = add(product(H, P, transpose(H)), R)
    Sinv = inverse(S)
    K = product(P, transpose(H), Sinv)
    A = add(identity(len(P)), product(K, H), -1)
    C = product(A, P)
    J, Jc = inverse(P), inverse(C)
    dC = product(A, dP, transpose(A))
    def tr(matrix):
        return sum((row[i] for i, row in enumerate(matrix)), F(0))
    loss = tr(product(J, dP, J, dP))-tr(product(Jc, dC, Jc, dC))
    force = product(A, dP, transpose(H), Sinv, residual)
    charge = product(transpose(force), Jc, force)[0][0]
    nis = product(transpose(residual), Sinv, residual)[0][0]
    assert 0 <= charge <= nis*loss/2
    return {'covariance_storage_loss': loss, 'gain_residual_energy': charge,
            'NIS': nis, 'linked_charge_upper': nis*loss/2}


def certificate():
    return {
        'qualification': 'OU3_PLANAR_COMPLETE_WORD_STORAGE_IDENTITY_V1',
        'result_type': 'PROVED — analytical identity; strict complete-word inequality OPEN',
        'controlling_inequality': 'W_N <= rho W_0 + c_d ||d||^2, rho<1, with every-prefix retention',
        'word': 'M=A_(N-1)...A_0; N=sum ordered_suffix_i B_i; one shared persistent root/history input',
        'exact_storage_identity': 'W_N-W_0=-v^T Q v+2 v^T X w+w^T Z w',
        'Q': 'J_0-M^T J_N M', 'X': 'M^T J_N N', 'Z': 'N^T J_N N',
        'exact_telescoping': 'Q=sum T_i^T (J_i-A_i^T J_(i+1) A_i) T_i',
        'regular_correction_identity': 'eplus^T C^-1 eplus=e^T P^-1 e-(H e)^T S^-1 H e+2 e^T P^-1 f+f^T C^-1 f, f=dK*r',
        'covariance_homogeneous_metric': '||C^-1/2 A dP A^T C^-1/2||_F<=||P^-1/2 dP P^-1/2||_F on regular optimal correction',
        'linked_gain_residual_bound': '||A dP H^T S^-1 r||_(C^-1)^2 <= (r^T S^-1 r)/2 * [||dP||_P^2-||A dP A^T||_C^2]',
        'linked_gain_residual_constant': '1/2',
        'constant_scope': 'unmasked optimal correction with actual same-operation P,H,R,r; dH/dR and held-BA terms retained separately',
        'derived_partial_joint_inequality': 'Wplus <= (1+NIS/lambda)*(Vmean-Lmean)+lambda*Vcov-((lambda-NIS)/2)*Lcov, lambda>0, for only the covariance-induced gain residual port',
        'linked_terms_required': ['dK*r', 'same-model dr and dH', 'dR and dF/dQ',
            'held-BA distinct gain/innovation rows', 'actual injection/reset differential',
            'AW divided differences and active-face inclusion', 'BA/BG projection faces',
            'moving P-normalized quotient chart', 'C_Q alpha and physical curvature',
            'literal nonlinear and arithmetic remainders in persistent history forcing',
            'inherited generator/reference/gates and every-prefix switching'],
        'strict_margin_requirement': 'Q >= epsilon J_0 uniformly on inherited same-history domain, epsilon>0; then bound the linked forcing quadratic',
        'derived_storage_inequality': 'W_N <= (1-epsilon/2) W_0 + w^T [Z+(2/epsilon) X^T J_0^-1 X] w',
        'uniform_epsilon': None,
        'failed_closure': {'classification': 'D_SUFFICIENT_BOUND_FAILURE',
            'exact_uncertified_quantity': 'the complete-word cross term from dK*r, mean-dependent rows, reset/AW/quotient derivatives in Q',
            'reason': 'available inherited domain does not bound nominal aw, actual innovations, precision, or their linked word action',
            'not_a_shipping_counterexample': True},
        'frontend_gain_is_MEKF_gain': False,
        'strict_complete_word_storage_verified': False,
        'radius_solve_performed': False,
        'structures_preserved': ['full joint quotient/covariance derivative', 'persistent gauge/source forcing',
            'actual ordered word and endpoint precision metrics', 'all linked cross terms'],
        'relaxations_introduced': ['fixed differentiable operation branches for local identities; face crossings need finite inclusions'],
        'all_time_magnetic_service_verified': False, 'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
