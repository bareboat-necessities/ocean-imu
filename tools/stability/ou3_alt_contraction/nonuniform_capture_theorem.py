"""Exact scalar comparison identities for the conditional ALT tail theorem.

These verify finite algebra, not any shipping infinite-history premise. See
ou3-alt-nonuniform-capture-theorem.md for the proof and quantifier boundary.
"""
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import binary32_covariance_coercivity as COV


def exact(value):
    if isinstance(value, (float, bool)):
        raise TypeError('use exact rational inputs, not floating-point evidence')
    return F(value)


def finite_comparison(initial, factors, supplies):
    initial = exact(initial)
    factors, supplies = tuple(map(exact, factors)), tuple(map(exact, supplies))
    if initial < 0 or len(factors) != len(supplies):
        raise ValueError('nonnegative initial storage and paired sequences required')
    if any(x < 0 for x in (*factors, *supplies)):
        raise ValueError('nonnegative factors and supplies required')
    product, forced = F(1), F(0)
    for rho, supply in zip(factors, supplies):
        product *= rho
        forced = rho*forced + supply
    return {'initial_multiplier': product, 'forced_response': forced,
            'storage_upper_bound': product*initial + forced,
            'finite_prefix_only': True, 'infinite_tail_certified': False}


def fixed_history_bound(initial, rho, supply, steps):
    initial, rho, supply = map(exact, (initial, rho, supply))
    if initial < 0 or supply < 0 or not 0 <= rho < 1:
        raise ValueError('nonnegative data and a strict history-specific factor required')
    if type(steps) is not int or steps < 0:
        raise ValueError('nonnegative integer word index required')
    product = rho**steps
    return {'storage_upper_bound': product*initial + supply*(1-product)/(1-rho),
            'conditional_ultimate_storage': supply/(1-rho),
            'history_uniform_rate_claimed': False,
            'shipping_tail_premises_proved': False}


def normalized_supply_bound(initial, factors, level):
    initial, level = exact(initial), exact(level)
    factors = tuple(map(exact, factors))
    if initial < 0 or level < 0 or any(not 0 <= rho <= 1 for rho in factors):
        raise ValueError('nonnegative storage/level and factors in [0,1] required')
    report = finite_comparison(initial, factors, ((1-rho)*level for rho in factors))
    product = report['initial_multiplier']
    if report['storage_upper_bound'] != product*initial + (1-product)*level:
        raise ArithmeticError('telescoping identity failed')
    report['product_tends_to_zero_proved'] = False
    return report


def selection_status():
    return {
        'qualification': 'OU3_ALT_HISTORY_DEPENDENT_CAPTURE_NONUNIFORM_TAIL_V1',
        'target': 'conditional history-wise eventual practical boundedness',
        'capture_quantifier': 'for each admitted history/initial state, some finite capture time',
        'capture_is_an_explicit_assumption': True,
        'common_capture_deadline_required': False,
        'same_history_tail_contraction_and_coercivity_required': True,
        'abstract_comparison_lemma_documented': True,
        'finite_capture_proved_for_shipping_filter': False,
        'shipping_nonlinear_word_inequality_proved': False,
        'history_wise_tail_coercivity_proved': False,
        'finite_SPD_binary32_implies_tail_coercivity': COV.bounds()['conditional_format_coercivity_lemma_proved'],
        'finite_SPD_binary32_tail_membership_proved': False,
        'history_wise_prefix_retention_proved': False,
        'variable_factors_below_one_alone_sufficient': False,
        'uniform_ISS_claimed': False,
        'lyapunov_stability_from_startup_claimed': False,
        'source_uniform_rho_certified': False,
        'storage_search_allowed': False,
        'ALT_STARTUP_PASS': False, 'ALT_LIVE_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
