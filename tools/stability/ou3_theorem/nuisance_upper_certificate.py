"""Rational constants for the recurring nuisance covariance upper comparison.

Proof: docs/ou3-nuisance-upper-proof.md. Enters covariance coercivity in
V=e'P^-1 e; this closes only the LIN/BA part of the upper comparison.
No contraction, AG bound or nonlinear retained radius is inferred.
"""
from fractions import Fraction as F
from .lin_path_certificate import small_x_source_defect


def interpolation_rows(a, b, d=F(0)):
    """Map three S observations at 0,a,a+b to (v,p,S) at a+b+d."""
    if min(a, b) <= 0 or d < 0:
        raise ValueError('positive observation gaps and nonnegative delay required')
    v = [2/(a*(a+b)), -2/(a*b), 2/(b*(a+b))]
    p = [b/(a*(a+b)), -(a+b)/(a*b), (a+2*b)/(b*(a+b))]
    s = [F(0), F(0), F(1)]
    return [v, [p[i]+d*v[i] for i in range(3)],
            [s[i]+d*p[i]+d*d*v[i]/2 for i in range(3)]]


def bounds():
    eps, _, _, _ = small_x_source_defect()
    hmin, hmax, tau, sig2 = F('.004'), F('.006'), F(12), F(16)
    gap, separation, horizon = F('.156'), F(8), F(17)
    # 1/(1-exp(-x)) <= 1+1/x. Keep the realized additive sync
    # increment fixed in the comparison, never recompute its positive part.
    aw = (1+eps)*sig2 + sig2*(1+tau/(2*hmin))
    fresh = (1+eps)*(2*sig2/F('.02'))*hmax*hmax*horizon/3
    assert aw <= 156**2 and fresh <= 1
    assert 2*(separation+gap)+gap < horizon
    wv = 4/separation**2
    wp0 = 4*(separation+gap)/separation**2
    wp = wp0+gap*wv
    ws = 1+gap*wp0+gap**2*wv/2
    ev = 156*horizon+1
    ep = 156*horizon**2/2+horizon
    es = 156*horizon**3/6+horizon**2/2
    sd = [ev+wv*(es+100), ep+wp*(es+100), es+ws*(es+100), F(156)]
    # Four LIN blocks and one BA block; Cauchy controls every cross block.
    diagonal = [5*x*x for x in sd]+[F(5,1600)]
    return aw, fresh, [wv, wp, ws], sd, diagonal


def certificate():
    aw, fresh, weights, sd, diagonal = bounds()
    return {
        'qualification': 'OU3_RECURRING_NUISANCE_UPPER_V1',
        'verified': True,
        'scope': 'real-arithmetic default regular A21; every operation after 17 s of regular S service; inherited construction AW/BA bounds; no frame/relock reconfiguration',
        'coordinate_order': ['v', 'p', 'S', 'a_w', 'b_a'],
        'covariance_comparison': 'P_nn <= diag(upper_diagonal tensor I_3), all nuisance cross blocks retained',
        'upper_diagonal': [str(x) for x in diagonal],
        'AW_unconditioned_variance_ceiling': str(aw),
        'fresh_neutral_noise_coefficient_squared': str(fresh),
        'interpolation_row_l1_ceilings': [str(x) for x in weights],
        'trial_error_standard_deviation_ceilings': [str(x) for x in sd]+['1/40'],
        'applied_S_gap_ceiling_s': '0.156',
        'window_s': '17',
        'nuisance_dimension': 15,
        'full_21_covariance_upper_verified': False,
        'source_uniform_A21_linear_dissipativity': False,
        'theorem_closed': False,
    }
