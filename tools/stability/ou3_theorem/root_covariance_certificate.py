"""Joint post-prediction covariance bound from two full Loewner comparisons.

This supplies coercivity for the complete-word energy identity; it is not an
information/loss certificate. The actual carried P is never restarted.
If P >= X and P >= Y, then P >= (X+Y)/2, even when X/Y are singular.
Combine the corrected 16-s LIN path factor with fresh one-step AG/BA process
noise at the *same* post-prediction root, before its first correction.
"""
from fractions import Fraction as F
from .lin_matrix_certificate import action_matrix
from .lin_path_certificate import inverse, rational_record
from .matrix_certificates import encoded, ldlt


def process_floors():
    hmin,hmax=F('.004'),F('.006')
    sigma_g,qbg,sbg=F('.00135'),F('1e-10'),F('.02')
    # Shipping isotropic Qtt >= sigma_g^2*h*I, Qbb=qbg*h*I.
    # Simpson B Qbg B' is PSD. The exact integral of B has norm <=h^2/2;
    # h^2 also covers the literal ||omega||<1e-7 polynomial branch.
    # Scale theta by 1 and bg by .02, then use a two-block norm bound.
    theta=sigma_g**2*hmin
    bias=qbg*hmin/sbg**2
    cross=qbg*hmax**2/sbg
    ag=min(theta,bias)-cross
    # q_ba = drive*tau/2*(1-exp(-2h/tau)); exp(-x)<=1/(1+x).
    drive,tau=F('2.5e-7'),F(5000)
    ba=drive*hmin/(1+2*hmin/tau)
    if min(ag,ba)<=0:
        raise ArithmeticError('fresh process bounds are not positive')
    return ag,ba


def certificate():
    a=action_matrix(); ag,ba=process_floors()
    c=[[v/2 for v in row] for row in inverse(a)]
    l,d=ldlt(c)
    return {
        'qualification':'OU3_JOINT_ROOT_COVARIANCE_CONVEX_FACTORS_V1',
        'verified':True,
        'scope':'real-arithmetic regular A21 post-prediction roots after a complete 16-s regular A21 LIN window; default bound shipping profile',
        'root_location':'after process prediction/PSD inflation, before the first current-sample correction',
        'comparisons':[
            'P >= diag(Q_AG,0,Q_BA) from the current prediction',
            'P >= E_LIN A_inverse E_LIN_transpose from the preceding path action',
        ],
        'combination_weight':'1/2',
        'joint_bound':'P >= diag(q_AG I_6, A_inverse tensor I_3, q_BA I_3)/2',
        'AG_coordinate_scales':['1','0.02'],
        'LIN_coordinate_scales':['2.4','18','132','4'],
        'BA_coordinate_scale':'1',
        'AG_covariance_diagonal_lower':rational_record(ag/2),
        'LIN_covariance_lower_matrix':encoded(c),
        'LIN_covariance_ldlt_unit_lower':encoded(l),
        'LIN_covariance_ldlt_diagonal':[str(x) for x in d],
        'BA_covariance_diagonal_lower':rational_record(ba/2),
        'root_metric_gamma':'1',
        'cross_covariance_role':'retained implicitly in both full Loewner comparisons; no cross blocks discarded',
        'source_constants':{'gyro_noise_density':'0.00135','gyro_bias_density':'1e-10',
                            'accel_bias_density':'2.5e-7','accel_bias_tau':'5000'},
        'recurring_post_prediction_lower_bound_verified':True,
        'uniform_covariance_upper_bound_verified':False,
        'every_prefix_retention_verified':False,
        'constructive_full_A21_mu_rho_enclosure':False,
        'float32_covariance_factor_verified':False,
        'theorem_closed':False,
    }
