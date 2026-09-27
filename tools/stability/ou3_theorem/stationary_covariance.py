"""Historical covariance upper on the exact quiet nominal source subcase.

This is the SAME historical reader, with two actual acc/mag groups. It
cancels arbitrary AG root covariance and all its nuisance cross terms.
It does not identify physical attitude/BA or prove stationary nonlinear ISS.
See docs/ou3-stationary-detectability.md for the complete scope and proof.
"""
from fractions import Fraction as F

from .nuisance_upper_certificate import bounds as nuisance_bounds


def quiet_covariance_bounds(*, gap_min=F(4, 125), gap_max=F(6, 125),
                            prefix_max=F(6, 125), gravity_min=F(9), field_min=F(9),
                            aw_std=F(156), ba_std=F(1, 40),
                            acc_std=F(1, 3), mag_std=F(1),
                            gyro_density=F(1, 500000), bias_density=F(1, 10**9)):
    """PSD comparisons in raw (theta,bg), then (theta,bg,v,p,S,aw,ba).

    Density arguments are covariance densities, NOT standard deviations.
    AW/BA ceilings concern the unconditioned auxiliary history after its
    actual root. They are not independently selected nominal mean boxes.
    All measurement/process factors and all within-root correlations remain.
    """
    lo, hi, delay, g, field, aw, ba, acc, mag, qg, qb = map(F, (
        gap_min, gap_max, prefix_max, gravity_min, field_min, aw_std, ba_std,
        acc_std, mag_std, gyro_density, bias_density))
    if min(lo, hi, g, field, acc, mag, qg, qb) <= 0 or min(delay, aw, ba) < 0 or lo > hi:
        raise ValueError('positive ordered gap/noise/geometry bounds required')
    # Each normalized attitude observation is theta+eta. AW and BA can
    # correlate arbitrarily; the two epochs need not be independent either.
    eta = [(aw + ba + acc)**2 / g**2] * 2 + [mag**2 / field**2]
    theta = [2 * e for e in eta]
    bg = [2 * (4 * e / lo**2 + qg / lo + qb * hi / 3) for e in eta]
    # Dropping intervening optimal corrections is conservative; quiet resets
    # are identity. Bound both principal blocks, then retain cross terms
    # through a full Loewner comparison with the factor two.
    theta_prefix = [2 * (a + delay**2 * b + qg * delay + qb * delay**3 / 3)
                    for a, b in zip(theta, bg)]
    bg_prefix = [2 * (b + qb * delay) for b in bg]
    nuisance = [x for value in nuisance_bounds()[4] for x in [value] * 3]
    return {'normalized_observation_error_upper': eta,
            'historical_AG_action_upper_diagonal': theta + bg,
            'every_prefix_AG_upper_diagonal': theta_prefix + bg_prefix,
            'full_21_upper_diagonal': [2 * x for x in theta_prefix + bg_prefix + nuisance]}


def certificate():
    from .corrected_word import prediction_floor
    bounds = quiet_covariance_bounds()
    # Existing positive full process comparison and bounded quiet F give
    # existence of epsilon>0 with Q>=epsilon F C_quiet F'. Do not turn
    # this coarse existence check into a useful numerical rho or radius.
    process = prediction_floor()
    return {
        'qualification': 'OU3_QUIET_NOMINAL_HISTORICAL_COVARIANCE_V1',
        'scope': 'real-arithmetic regular A21, zero corrected rate/AW/injections, identity attitude, fixed B e_x; actual coincident acc/mag groups eight qualified predictions apart; inherited AW/BA bounds and 17-second nuisance comparison',
        'profile': {'gravity_min': '9', 'field_min': '9', 'acc_std_max': '1/3',
                    'mag_std_max': '1', 'gyro_covariance_density_max': '1/500000',
                    'gyro_bias_covariance_density_max': '1/1000000000',
                    'sample_gap_min_s': '.004', 'sample_gap_max_s': '.006'},
        'process_lower_profile': 'existing corrected_word.prediction_floor default theorem profile; positive process densities retained',
        **{k: [str(x) for x in v] for k, v in bounds.items()},
        'all_AG_root_and_cross_covariance_cancelled': True,
        'uniform_quiet_nominal_historical_action': True,
        'uniform_quiet_nominal_full_covariance_upper': True,
        'quiet_nominal_every_operation_covariance_upper': True,
        'quiet_nominal_homogeneous_linear_loss_exists': process > 0,
        'numerical_useful_rho_certified': False,
        'stationary_compatible_class_practical_stability': False,
        'source_uniform_moving_A21_covariance_upper': False,
        'nonlinear_every_prefix_retention': False,
        'float32_totality': False,
        'theorem_closed': False,
    }
