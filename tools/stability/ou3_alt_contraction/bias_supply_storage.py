"""Conditional full-storage input cross-term absorption, without metric fitting.

For x=Uu+Ss in SPD M0, ||Uu||_M0 <= ||x||_M0+||Ss||_M0.
A word norm estimate ||F||_M1 <= (a+eps)||Uu||_M0+b D and an
INDEPENDENT bound ||Ss||_M0 <= k D imply the coefficients below.
No native operator, remainder, projection or covariance preservation is proved
by evaluating these rational coefficients.
"""
from fractions import Fraction as F


def exact_nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError('exact rational coefficient required')
    value = F(value)
    if value < 0:
        raise ValueError('nonnegative coefficient required')
    return value


def coefficients(tangent_norm, remainder_norm, bias_embedding_gain,
                 other_supply_gain, young):
    a, eps, k, b, eta = map(exact_nonnegative, (
        tangent_norm, remainder_norm, bias_embedding_gain, other_supply_gain, young))
    if eta == 0:
        raise ValueError('strictly positive Young parameter required')
    alpha = a+eps
    return {'rho': (1+eta)*alpha*alpha,
            'C': (1+1/eta)*(alpha*k+b)**2,
            'native_premises_proved': False}


def bias_coordinate_squared_bound(estimate_radius, true_bias_radius):
    """||s||² <= (R+B)²+B² for s=(beta-ba_hat,beta).

    Requires an independently proved ||ba_hat||<=R and ||beta||<=B.
    A constant R supplies a practical floor; it is not zero-input ISS decay.
    """
    R, B = map(exact_nonnegative, (estimate_radius, true_bias_radius))
    return (R+B)**2+B*B
