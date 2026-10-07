"""Exact innovation/storage algebra for the inherited-domain proof.

No replay, finite difference, arbitrary covariance domain or contraction fit is
used. The additive correction is an intermediate operation: literal attitude
injection, reset, projection, AW and moving quotient maps must follow it.
"""
from __future__ import annotations

from fractions import Fraction as F
import struct

from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, transpose
from .planar_linked_riccati_mean import (
    gain_differential, optimal_covariance_differential, product,
)


def energy(x, metric):
    return product(transpose(x), metric, x)[0][0]


def innovation_variational_identity(P, H, R, residual, comparison):
    """NIS=min_z z'Jz+(r-Hz)'R^-1(r-Hz), with exact square gap.

    H is the INNOVATION row. In held BA, the algebraic minimizer is generally
    not the applied masked gain. This identity does not replace that gain.
    """
    ldlt(P)
    ldlt(R)
    J, Ri = inverse(P), inverse(R)
    S = add(product(H, P, transpose(H)), R)
    minimizer = product(P, transpose(H), inverse(S), residual)
    remainder = add(residual, product(H, comparison), -1)
    precision = add(J, product(transpose(H), Ri, H))
    nis = energy(residual, inverse(S))
    candidate = energy(comparison, J)+energy(remainder, Ri)
    gap = energy(add(comparison, minimizer, -1), precision)
    assert candidate-nis == gap
    return {'NIS': nis, 'comparison_upper': candidate, 'square_gap': gap,
            'minimizer': minimizer, 'remainder': remainder}


def physical_acceleration_remainder_form(dimension, theta_indices, aw_indices,
                                         gravity_plus_physical_acceleration, split):
    """B such that ||r-H e-noise|| <= e'B e for the default CoG model.

    e uses a LEFT rotation vector and physical-minus-nominal AW/BA. Exact
    rotation Taylor remainder gives (g+A)|theta|^2/2+|theta| |delta_aw|.
    Young's inequality with split>0 bounds the mixed term. A is a PHYSICAL
    acceleration cap, not a substitution for the unknown nominal AW state.
    """
    split = F(split)
    g_a = F(gravity_plus_physical_acceleration)
    if split <= 0 or g_a < 0:
        raise ValueError('positive split and nonnegative physical force cap required')
    if len(theta_indices) != 3 or len(aw_indices) != 3:
        raise ValueError('complete attitude and AW triples required')
    if len(set(theta_indices+aw_indices)) != 6 or any(
            i < 0 or i >= dimension for i in theta_indices+aw_indices):
        raise ValueError('disjoint in-range attitude and AW indices required')
    B = [[F(0) for _ in range(dimension)] for _ in range(dimension)]
    for i in theta_indices:
        B[i][i] = (g_a+split)/2
    for i in aw_indices:
        B[i][i] = 1/(2*split)
    return B


def innovation_domain_bound(P, remainder_form, kappa, V_cap, noise_cap, R_floor):
    """Exact sufficient NIS bound V+(kappa V+noise)^2/R_floor.

    Unlike a triangle bound, the variational identity has no cross term
    between state precision energy and the residual remainder. The caller
    must establish V_cap, R_floor and kappa*P^-1 >= B on EVERY inherited
    prefix. A rational operand check is an identity/domain regression only.
    No supplied cap is claimed reachable or invariant by this function.
    """
    kappa, V_cap, noise_cap, R_floor = map(F, (kappa, V_cap, noise_cap, R_floor))
    ldlt(P)
    if min(kappa, V_cap, noise_cap) < 0 or R_floor <= 0:
        raise ValueError('nonnegative caps and strictly positive noise floor required')
    if not is_psd(add([[kappa*x for x in row] for row in inverse(P)], remainder_form, -1)):
        raise ValueError('D_SUFFICIENT_BOUND_FAILURE: linked remainder domination not established')
    return V_cap+(kappa*V_cap+noise_cap)**2/R_floor


def information_shear_correction(P, H, R, e, residual, dP, de, dH, dR, dr):
    """Exact differential in eta=de-dP P^-1 e; comparison is not reseeded.

    Keeps dH, dR and the literal residual differential through dr+H*de.
    e is an additive pre-injection comparison coordinate, not a claim that
    physical attitude errors add on SO(3). Chart/injection terms are separate.
    """
    ldlt(P)
    ldlt(R)
    S = add(product(H, P, transpose(H)), R)
    K = product(P, transpose(H), inverse(S))
    A = add(identity(len(P)), product(K, H), -1)
    C = product(A, P)
    dK = gain_differential(P, H, K, inverse(S), dP, dH, dR)
    dC = optimal_covariance_differential(P, H, K, dP, dH, dR)
    eplus = add(e, product(K, residual))
    deplus = add(add(de, product(dK, residual)), product(K, dr))
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(deplus, product(dC, inverse(C), eplus), -1)
    a = add(product(H, e), residual)
    mismatch = add(dr, product(H, de))
    linked = product(A, eta)
    linked = add(linked, product(C, transpose(dH), inverse(R), a))
    linked = add(linked, product(K, dH, e))
    linked = add(linked, product(K, dR, inverse(R), a), -1)
    linked = add(linked, product(K, mismatch))
    assert linked == direct
    return {'eta': eta, 'eta_plus': direct, 'homogeneous': product(A, eta),
            'linked_row_noise_chart_force': add(linked, product(A, eta), -1),
            'posterior_covariance': C, 'posterior_covariance_tangent': dC}


def shear_coercivity_factor(comparison_energy_cap, covariance_weight):
    """C^-1 Wbase <= Wshear <= C Wbase, with actual same-P comparison cap."""
    E, weight = F(comparison_energy_cap), F(covariance_weight)
    if E < 0 or weight <= 0:
        raise ValueError('nonnegative comparison energy and positive covariance weight required')
    return max(F(2), 1+2*E/weight)


def information_shear_prediction(P, Fmap, Q, e, source, dP, dF, dQ, de, dsource):
    """Exact prediction shear, including physical OU mismatch and dF/dQ."""
    Pnext = add(product(Fmap, P, transpose(Fmap)), Q)
    ldlt(P)
    ldlt(Pnext)
    enext = add(product(Fmap, e), source)
    dPnext = add(add(add(product(dF, P, transpose(Fmap)),
                         product(Fmap, dP, transpose(Fmap))),
                     product(Fmap, P, transpose(dF))), dQ)
    denext = add(add(product(dF, e), product(Fmap, de)), dsource)
    eta = add(de, product(dP, inverse(P), e), -1)
    gnext = product(inverse(Pnext), enext)
    direct = add(denext, product(dPnext, gnext), -1)
    linked = add(product(Fmap, eta), product(Fmap, dP,
        add(product(inverse(P), e), product(transpose(Fmap), gnext), -1)))
    linked = add(linked, product(dF, add(e, product(P, transpose(Fmap), gnext), -1)))
    linked = add(linked, product(Fmap, P, transpose(dF), gnext), -1)
    linked = add(add(linked, dsource), product(dQ, gnext), -1)
    assert linked == direct
    return {'eta_plus': direct, 'homogeneous': product(Fmap, eta),
            'linked_process_force': add(linked, product(Fmap, eta), -1)}


def linked_prediction_charge(P, Fmap, Q, e, dP):
    """Fixed F,Q, zero source ONLY: shear port energy <= L_mean*L_P/2.

    The actual source/dF/dQ terms above remain in the complete word. This
    identity cannot impose the estimator OU law on physical acceleration.
    """
    ldlt(P)
    if not is_psd(Q) or dP != transpose(dP):
        raise ValueError('PSD process covariance and symmetric tangent required')
    Pnext = add(product(Fmap, P, transpose(Fmap)), Q)
    ldlt(Pnext)
    J, Jnext = inverse(P), inverse(Pnext)
    propagated = product(Fmap, e)
    dPnext = product(Fmap, dP, transpose(Fmap))
    force = product(Fmap, dP, add(product(J, e), product(transpose(Fmap), Jnext, propagated), -1))
    mean_loss = energy(e, J)-energy(propagated, Jnext)
    def tr(matrix):
        return sum((row[i] for i, row in enumerate(matrix)), F(0))
    covariance_loss = tr(product(J, dP, J, dP))-tr(product(Jnext, dPnext, Jnext, dPnext))
    charge = energy(force, Jnext)
    assert 0 <= charge <= mean_loss*covariance_loss/2
    return {'mean_process_loss': mean_loss, 'covariance_process_loss': covariance_loss,
            'sheared_covariance_port_energy': charge,
            'linked_charge_upper': mean_loss*covariance_loss/2}


def held_ba_reachable_correction(P_active, P_ba, H_active, R):
    """Active 18-state block of the literal held-BA Joseph correction.

    Only on the reachable P_active,BA=0 manifold. All three BA means and
    covariance entries remain inherited; effective R=R+P_ba is linked and
    dR_eff=dR+dP_ba. Not a reduction of the full 21-state theorem target.
    """
    ldlt(P_active)
    ldlt(P_ba)
    ldlt(R)
    effective_R = add(R, P_ba)
    S = add(product(H_active, P_active, transpose(H_active)), effective_R)
    B = product(P_active, transpose(H_active))
    K = product(B, inverse(S))
    C = add(add(add(P_active, product(K, transpose(B)), -1),
                product(B, transpose(K)), -1), product(K, S, transpose(K)))
    assert C == add(P_active, product(K, transpose(B)), -1)
    return {'effective_R': effective_R, 'innovation': S, 'active_K': K,
            'active_C': C, 'held_P_ba': P_ba}


def inherited_handoff_bounds():
    """Phase-uniform real-operation bound at the literal FIRST Live boundary.

    Qualified implication: the default planar private-angle/reference cone
    must hold and handoff must occur. This neither proves finite acquisition
    nor transfers the bound through subsequent predictions/corrections or
    binary32 covariance/reference arithmetic. Physical S starts at zero at
    construction and is carried, never zeroed at handoff.
    """
    def represented(value):
        return F(struct.unpack('f', struct.pack('f', value))[0])
    g, c = F(196133, 20000), F(39999, 40001)
    angle = F(11, 5000)
    frequency_upper = F(11, 35)  # pi/10 < 11/35
    amplitude = F(1, 50)
    acceleration = amplitude*frequency_upper**2
    velocity = amplitude*frequency_upper
    primitive = F(2, 15)  # 2*amplitude/(pi/10) < 2/15, with S(0)=0
    terms = {
        'attitude': angle**2/represented(.035)**2,
        'gyro_bias': F(0),
        'velocity': velocity**2,
        'position': amplitude**2/F(400),
        'physical_S': primitive**2/F(2500),
        'world_acceleration': acceleration**2/represented(.05)**2,
        'slow_BA_central_fibre': (g*(1-c))**2/represented(.004)**2,
    }
    V = sum(terms.values(), F(0))
    gravity_mismatch = abs(represented(9.80665)-g)
    acc_remainder = ((g+acceleration)*angle**2/2+angle*acceleration
                     +(1+angle)*gravity_mismatch)
    ref_difference = 75*angle*(1+angle**2/8)
    mag_remainder = 75*angle**2/2+(1+angle)*ref_difference
    nis_acc = V+acc_remainder**2/represented(.2)**2
    nis_mag = V+mag_remainder**2/represented(.8)**2
    assert V < F(21, 1000)
    assert nis_acc < F(21, 1000)
    assert nis_mag < F(13, 200)
    return {'comparison': 'beta=0 central physical fibre; no physical/estimator restart',
            'boundary': 'after first goLive/enterLive and before next estimator operation',
            'arithmetic': 'real operation lift with represented configuration scalars; reference averaging real',
            'premises': ['default planar source and no external state/covariance writes',
                'private pitch error and reference cone with a=11/5000 at handoff',
                'zero planar yaw and default S_factor=1', 'physical S(0)=0 at construction'],
            'precision_energy_terms': {k: str(v) for k, v in terms.items()},
            'precision_energy_upper_exact': str(V),
            'precision_energy_strict_upper': '21/1000',
            'represented_gravity_mismatch': str(gravity_mismatch),
            'accelerometer_NIS_upper_exact': str(nis_acc),
            'accelerometer_NIS_strict_upper': '21/1000',
            'magnetic_NIS_upper_exact': str(nis_mag),
            'magnetic_NIS_strict_upper': '13/200',
            'NIS_scope': 'hypothetical immediate boundary correction using boundary operands, not the following sample',
            'finite_handoff_time_verified': False,
            'post_handoff_prefix_invariance_verified': False,
            'binary32_reference_covariance_transfer_verified': False}


def certificate():
    return {
        'qualification': 'OU3_PLANAR_INNOVATION_STORAGE_V1',
        'result_type': 'PROVED — analytical operation identities and conditional domain bounds; inherited domain OPEN',
        'controlling_inequality': 'positive complete-word gap before any radius or all-future service floor',
        'NIS_identity': 'NIS=min_z z^T P^-1 z+(r-Hz)^T R^-1(r-Hz)',
        'NIS_square_gap': '(z-Kr)^T (P^-1+H^T R^-1 H) (z-Kr)',
        'NIS_domain_bound': 'NIS<=V+(kappa V+d)^2/R_floor when remainder<=e^T B e+d and B<=kappa P^-1',
        'acceleration_remainder': '||(exp([theta])-I-[theta]) Rhat(a_true-g)+[theta] Rhat delta_aw|| <= (g+Aphysical)|theta|^2/2+|theta||delta_aw|',
        'held_BA_reachable_manifold': 'P_active,BA=0; P_BA constant while held; default held BA mean remains zero',
        'held_BA_effective_noise': 'R_eff=R_acc+P_BA; dR_eff=dR_acc+dP_BA; held BA mean error remains a residual/source port',
        'held_BA_scope': 'regular default CoG branch from construction or explicit hold decoupling, before release, without external covariance writes',
        'inherited_first_handoff': inherited_handoff_bounds(),
        'information_shear': 'eta=de-dP P^-1 e=P d(P^-1 e)',
        'exact_sheared_correction': 'eta+=A eta+C dH^T R^-1 a+K dH e-K dR R^-1 a+K(dr+H de), a=H e+r',
        'fixed_row_covariance_residual_cancellation': 'dH=dR=0 and dr=-H de imply eta+=A eta for arbitrary same-operation residual',
        'coercivity_factor': 'max(2,1+2 E/lambda) when e^T P^-1 e<=E and lambda>0',
        'exact_sheared_prediction': 'F eta+F dP(J e-F^T Jnext enext)+dF(e-P F^T Jnext enext)-F P dF^T Jnext enext+dsource-dQ Jnext enext',
        'linked_process_port_bound': 'fixed F,Q and zero source: ||f_P||_(Pnext^-1)^2<=L_mean L_P/2; actual physical source/dF/dQ retained separately',
        'comparison_reseed_allowed': False,
        'physical_attitude_additive_assumption': False,
        'remaining_ports': ['mean-dependent dH', 'same-history dR and dP_BA in H18',
            'prediction/source reference transport', 'literal attitude injection and reset',
            'BA/BG projection', 'AW faces and pending targets', 'moving quotient projector',
            'physical gauge and nonlinear fibre curvature', 'arithmetic'],
        'structures_preserved': ['actual innovation rows, P, R and residual',
            'held masked gain and literal Joseph', 'all 21 coordinates with invertible tangent shear',
            'same-history comparison and residual dependence', 'physical compatibility gauge'],
        'relaxations_introduced': ['regular real-operation branch',
            'Young quadratic bound on the exact rotation remainder; no inherited cap supplied'],
        'source_uniform_nominal_domain_verified': False,
        'uniform_NIS_cap': None, 'uniform_complete_word_epsilon': None,
        'radius_solve_performed': False, 'all_time_magnetic_service_verified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
