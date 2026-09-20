"""Coupled covariance reduction and actual-gain finite-error word bounds.

Enters sqrt(V_end) <= sqrt(1-delta)*sqrt(V_root) + sum(operation supplies).
The nuisance constants are source bounds. The missing six-column corrected
loss is NOT supplied by this module; a conditional implication is not a
source-uniform contraction certificate. See docs/ou3-corrected-word-proof.md.
"""
from fractions import Fraction as F

from .lin_matrix_certificate import SCALES, action_matrix
from .lin_path_certificate import inverse, rational_record, small_x_source_defect
from .matrix_certificates import add, encoded, identity, is_psd, ldlt
from .nuisance_upper_certificate import bounds
from .root_covariance_certificate import process_floors
from .word_energy import full_loss_margin, word_identity


def nuisance_root_bounds():
    """Embedded nuisance floor at the root immediately BEFORE prediction.

    Each 4x4 LIN matrix is tensored with I3. An acc correction has nuisance
    H=[0,0,0,R,I], hence H'Racc^-1 H <= 2/sigma_acc,min^2 on AW and BA.
    At most one acc and one S correction occur per prediction. Magnetic
    corrections and attitude resets preserve an embedded nuisance floor.
    """
    a = action_matrix()
    scales = tuple(map(F, SCALES))
    precision = [[2*a[i][j]/(scales[i]*scales[j])
                  for j in range(4)] for i in range(4)]
    precision[2][2] += 1/F('.075')**2
    precision[3][3] += 2/F('.05')**2
    _, ba = process_floors()
    ba_precision = 2/ba + 2/F('.05')**2
    upper = bounds()[-1]
    # Bound lambda_max(U^1/2 L^-1 U^1/2) by its 4x4 trace.
    # The three identical axes do not introduce a factor of three.
    alpha = 1/max(sum(upper[i]*precision[i][i] for i in range(4)),
                  upper[4]*ba_precision)
    lower = inverse(precision)
    test = [[lower[i][j] - (alpha*upper[i] if i == j else 0)
             for j in range(4)] for i in range(4)]
    if not 0 < alpha < 1 or not is_psd(test) or 1/ba_precision < alpha*upper[4]:
        raise ArithmeticError('embedded nuisance lower/upper comparison failed')
    return lower, 1/ba_precision, upper, alpha


def prediction_floor():
    """Positive raw-coordinate Q floor; existence bound, not a useful radius.

    The integrated LIN Gram is full rank. Natural step scales are
    (h,h^2,h^3,1), not three driving-noise columns at a single instant.
    """
    eps, _, b0, _ = small_x_source_defect()
    b0_inverse_norm = max(sum(abs(v) for v in row) for row in inverse(b0))
    hmin, hmax, tau_min, tau_max = F('.004'), F('.006'), F('.02'), F(12)
    lin = ((1-eps)*F('.05')**2*(hmin/tau_max)*hmin**6
           / ((1+hmax/tau_min)**2*b0_inverse_norm))
    ag, ba = process_floors()
    q = min(ag*F('.02')**2, ba, lin)
    if q <= 0:
        raise ArithmeticError('nonpositive full process floor')
    return q


def conditional_scalar_margin(mu, alpha, nuisance_ceiling, q, f_norm):
    """Consequence ONLY if D_AG,AG >= mu I at this pre-prediction root.

    P >= E_n L E_n', P_nn <= U <= c I, L >= alpha U imply
    P <= (1/(alpha*mu)+c) I. The first prediction gives the displayed loss.
    Prefer the matrix version in the proof when certifying an actual margin.
    """
    mu, alpha, c, q, f = map(F, (mu, alpha, nuisance_ceiling, q, f_norm))
    if min(mu, alpha, c, q, f) <= 0 or alpha >= 1:
        raise ValueError('positive premises and 0 < alpha < 1 required')
    cap = 1/(alpha*mu)+c
    return q/(f*f*cap+q)


def reset_remainder_bound(injection_norm, error_difference_norm, chart_radius):
    """Exact exponential reset vs shipping G=I+[injection]x/2.

    Bounds Log(Exp(theta) Exp(-d))-G(d)(theta-d); |theta-d| <= r < 2.
    The real small-angle polynomial defect is charged separately below.
    Float32 trigonometric/normalization errors are not included.
    """
    d, v, r = map(F, (injection_norm, error_difference_norm, chart_radius))
    if min(d, v, r) < 0 or v > r or r >= 2:
        raise ValueError('nonnegative norms and |theta-d| <= r < 2 required')
    return d*d*v/6 + (F(1, 2)/(2-r)+F(1, 4))*v*v


def polynomial_injection_defect(injection_norm):
    """Angular defect of normalized literal degree-four quaternion branch."""
    d = F(injection_norm)
    if d < 0:
        raise ValueError('nonnegative injection norm required')
    if d >= F('.01'):
        return F(0)  # exact-real sin/cos branch; arithmetic remains separate
    return 4*(d**6/46080+d**7/645120)


def shipping_reset_remainder_bound(injection_norm, error_difference_norm, chart_radius):
    d, v, r = map(F, (injection_norm, error_difference_norm, chart_radius))
    base = reset_remainder_bound(d, v, r)
    defect = polynomial_injection_defect(d)
    if r+defect >= 2:
        raise ValueError('implemented injection leaves the certified log chart')
    return base+2*defect/(2-r-defect)


def projection_storage_guard():
    """Source BA marginal turns the existing Euclidean guard into a V guard.

    At the literal pre-projection boundary, |e_ba| <= sqrt(P_ba,ba V)
    <= sqrt(V)/40. This bound retains every cross covariance and uses the
    inherited raw BA marginal, not the five-block Cauchy comparison. The
    guard must hold at every such prefix; it is not an invariance assertion.
    """
    import json
    from pathlib import Path
    constants = json.loads(Path(__file__).with_name('constants.json').read_text())
    # The contract stores physical and estimator bias bounds separately.
    physical = F(str(constants['imu_bias']['B_a_mps2']))
    radius = F('0.4')
    gap = radius-physical
    if gap <= 0:
        raise ValueError('projection interior is empty')
    return {'physical_bias_bound': str(physical), 'literal_projection_radius': str(radius),
            'BA_marginal_variance_ceiling': '1/1600',
            'sqrt_V_strict_ceiling': str(40*gap),
            'certified_test_radius': '6', 'remaining_bias_distance_at_test_radius': str(gap-F(6, 40)),
            'pre_projection_prefix_required': True,
            'projection_defect_under_guard': '0',
            'prefix_retention_proved': False, 'explicit_retained_region_proved': False}


def coupled_example():
    """Exact adversarial cross-covariance check, not a shipping history."""
    p = [[F(100), F(0), F(10)], [F(0), F(2), F(0)], [F(10), F(0), F(3)]]
    f = identity(3)
    f[0][1] = F(1, 200)
    q = [[x/100 for x in row] for row in identity(3)]
    word = word_identity(p, [
        {'kind': 'prediction', 'F': f, 'Q': q},
        {'kind': 'correction', 'H': [[1, 0, 1], [0, 1, 0]], 'R': identity(2)},
    ])
    mu, alpha = F(1, 125), F(1, 3)
    ldlt(add([r[:2] for r in word['loss'][:2]], identity(2), -mu))
    embedded = [[F(0) for _ in range(3)] for _ in range(3)]
    embedded[2][2] = 1
    ldlt(add(p, embedded, F(-1)))
    delta = conditional_scalar_margin(mu, alpha, 3, F('.01'), F('1.005'))
    if not full_loss_margin(word, delta):
        raise ArithmeticError('conditional full-loss consequence failed')
    return {'root_covariance': encoded(p), 'AG_loss_floor': str(mu),
            'nuisance_ratio': str(alpha), 'strict_full_loss_margin': str(delta),
            'exact_full_loss_check': True, 'shipping_history': False}


def certificate():
    lower, ba, upper, alpha = nuisance_root_bounds()
    return {
        'qualification': 'OU3_CORRECTED_WORD_REDUCTION_V1',
        'verified': True,
        'scope': 'regular default real-arithmetic A21 after 17 s; roots immediately before prediction; inherited construction nuisance bounds',
        'LIN_embedded_lower_before_prediction': encoded(lower),
        'BA_embedded_lower_before_prediction': str(ba),
        'nuisance_upper_diagonal': list(map(str, upper)),
        'embedded_nuisance_ratio': rational_record(alpha),
        'full_process_scalar_floor': rational_record(prediction_floor()),
        'prediction_operator_norm_ceiling': '1.012',
        'matrix_implication': 'D_AG,AG >= J > 0 implies P <= diag((1+eta)/alpha J^-1, (1+1/eta) U), eta > 0',
        'remaining_uniform_premise': 'six AG columns of actual complete corrected word loss, in fixed raw coordinates, have a common SPD lower bound J',
        'exact_coupled_example': coupled_example(),
        'actual_gain_finite_error_word_composition_proved': True,
        'joint_prediction_measurement_input_action_proved': True,
        'reset_comparison': 'exact finite-angle log reset, with normalized source small-angle polynomial defect',
        'projection_storage_guard': projection_storage_guard(),
        'independent_filter_trajectories_required': False,
        'six_column_source_uniform_loss_verified': False,
        'full_21_covariance_upper_verified': False,
        'source_uniform_A21_linear_dissipativity': False,
        'explicit_nonlinear_retained_region_verified': False,
        'float32_transfer_verified': False,
        'theorem_closed': False,
    }
