"""Sharp principal AW covariance ceiling under the isotropic stationary sync.

Proof: docs/ou3-world-frame-rows.md, section 7. Scope is the real-arithmetic
default profile of docs/ou3-nuisance-upper-proof.md: S_factor=1, additive
pending sync, tau<=12 s, Sigma_aw=sigma^2 I<=16 I, no frame/relock
reconfiguration. The ceiling enters Corollary A through
|e_aw|^2<=lambda_max(P_aw) V; it is not propagated into the nuisance upper
comparison. Non-promoting.
"""
from fractions import Fraction as F

from .lin_path_certificate import small_x_source_defect
from .matrix_certificates import add, identity, is_psd, matmul, transpose
from .world_frame import decimal_lower, nominal_attitude_column_floor, quaternion_rotation

SIGMA2_MAX = F(16)             # MAX_SIGMA_A=4 m/s^2, S_factor=1
INITIAL_VARIANCE = F(121, 25)  # construction seed 2.2^2 I
TAU_MAX = F(12)
INHERITED_STD = F(156)


def ceiling(eps, sigma2_max=SIGMA2_MAX, initial=INITIAL_VARIANCE):
    """Invariant ceiling of lambda_max(P_aw) over every regular operation."""
    if eps < 0 or sigma2_max < 0 or initial < 0:
        raise ValueError('nonnegative defect, variance ceiling and seed required')
    return max(F(initial), (1+F(eps))*F(sigma2_max))


def prediction_upper(m, phi2, sigma2, eps):
    """lambda_max after P_aw <- phi^2 P_aw+Q_aa with Q_aa<=(1+eps)(1-phi^2)sigma^2 I."""
    if not 0 <= phi2 <= 1 or m < 0 or sigma2 < 0 or eps < 0:
        raise ValueError('phi^2 in [0,1] and nonnegative m, sigma^2, eps required')
    return phi2*m+(1+eps)*(1-phi2)*sigma2


def diagonal(values):
    return [[F(values[i]) if i == j else F(0) for j in range(len(values))] for i in range(len(values))]


def outer(u, v, weight=F(1)):
    return [[weight*a*b for b in v] for a in u]


def isotropic_sync_witness():
    """P+Pi_+(sigma^2 I-P) is the spectral max of P and sigma^2 I (exact)."""
    q = quaternion_rotation((5, 1, -2, 3))
    beta, sigma2 = [F(1, 2), F(3), F(20)], F(4)
    p = matmul(q, matmul(diagonal(beta), transpose(q)))
    # sigma^2 I-P=Q diag(sigma^2-beta) Q', so its positive part is exact.
    delta = matmul(q, matmul(diagonal([max(F(0), sigma2-b) for b in beta]), transpose(q)))
    after = add(p, delta)
    expected = matmul(q, matmul(diagonal([max(b, sigma2) for b in beta]), transpose(q)))
    assert after == expected
    assert is_psd(add(after, identity(3), -sigma2))
    assert is_psd(add(identity(3), after, -F(1, 20))) and is_psd(add(p, identity(3), -F(1, 2)))
    return {'lambda_max_before': str(max(beta)), 'sigma2': str(sigma2),
            'lambda_max_after': str(max(max(beta), sigma2)), 'post_sync_floor_sigma2': True}


def anisotropic_sync_counterexample():
    """S_factor=2: the pending sync can raise lambda_max above both operands.

    Sigma=diag(16,16,4); P=Sigma+q1 q1'-2 q2 q2' in the x-z plane has
    0<=P<=16 I, yet P+Pi_+(Sigma-P)=Sigma+q1 q1' has (1,1) entry 16+9/25.
    """
    sigma = diagonal([16, 16, 4])
    q1, q2 = [F(3, 5), F(0), F(4, 5)], [F(-4, 5), F(0), F(3, 5)]
    x = add(outer(q1, q1), outer(q2, q2), F(-2))
    p = add(sigma, x)
    after = add(sigma, outer(q1, q1))  # X has eigenpairs (1,q1),(-2,q2),(0,e_y)
    assert is_psd(p) and is_psd(add(identity(3), p, -F(1, 16)))
    assert after[0][0] > 16
    return {'S_factor': '2', 'sigma_diagonal': ['16', '16', '4'],
            'lambda_max_operands_upper': '16', 'after_x_diagonal': str(after[0][0]),
            'isotropy_required': True}


def storage_radius(threshold, eps, sigma2):
    """Rational r<=threshold/sqrt((1+eps)sigma^2), using sqrt(1+eps)<=1+eps/2."""
    root_sigma = {F(16): F(4), F(1): F(1)}.get(F(sigma2))
    if root_sigma is None:
        raise ValueError('exact square root of sigma^2 required')
    return threshold/(root_sigma*(1+eps/2))


def certificate():
    eps = small_x_source_defect()[0]
    m = ceiling(eps)
    # Linear in phi^2: equality at every phi^2 in [0,1] from both endpoints.
    for phi2 in (F(0), F(1, 2), F(1)):
        assert prediction_upper(m, phi2, SIGMA2_MAX, eps) <= m
    # Excess over a level s>=(1+eps)sigma^2 contracts by phi^2 per prediction.
    level, start, phi2 = (1+eps)*F(1, 25), F(3), F(99, 100)
    assert max(F(0), prediction_upper(start, phi2, F(1, 25), eps)-level) <= phi2*(start-level)
    threshold = nominal_attitude_column_floor(16, 0, F(1, 5))['aw_error_threshold']
    clamp_radius = storage_radius(threshold, eps, SIGMA2_MAX)
    unit_radius = storage_radius(threshold, eps, 1)
    return {
        'qualification': 'OU3_AW_COVARIANCE_CEILING_V1',
        'scope': 'real-arithmetic default regular profile: S_factor=1, additive pending sync, tau<=12 s, Sigma_aw=sigma^2 I<=16 I, no frame/relock reconfiguration',
        'process_covariance_relative_defect': str(eps),
        'aw_variance_ceiling': str(m),
        'aw_standard_deviation_ceiling_upper_mps2': str(4*(1+eps/2)),
        'inherited_aw_standard_deviation_ceiling_mps2': str(INHERITED_STD),
        'standard_deviation_improvement_factor_lower': str(decimal_lower(INHERITED_STD/(4*(1+eps/2)))),
        'prediction_step_invariant': True,
        'corrections_decrease_aw_marginal': True,
        'resets_and_projections_leave_aw_block': True,
        'isotropic_sync_witness': isotropic_sync_witness(),
        'anisotropic_sync_counterexample': anisotropic_sync_counterexample(),
        'excess_contraction_rate_per_s': str(2/TAU_MAX),
        'post_sync_floor_makes_ceiling_tight': True,
        'covariance_ceiling_below_stationary_variance_possible': False,
        'corollary_A_threshold_16s_mps2': str(threshold),
        'corollary_A_storage_radius_lower_at_sigma_clamp': str(decimal_lower(clamp_radius)),
        'corollary_A_storage_radius_lower_at_unit_sigma': str(decimal_lower(unit_radius)),
        'inherited_storage_radius': str(threshold/INHERITED_STD),
        'propagated_to_nuisance_upper_comparison': False,
        'uniform_AW_tracking_bound': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
