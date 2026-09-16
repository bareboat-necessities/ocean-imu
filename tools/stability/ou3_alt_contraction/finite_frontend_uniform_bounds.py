"""All-time finite supplies for the configured OU-III sigma frontend.

This induction supplies the LPF, stillness, band, variance and sigma operations
inside a complete startup/Live event. It does not prove startup capture or
MEKF totality. The input is the SAME Mahony vertical already bounded by the
physical MEMS contract; no reset, trace fit, or startup deadline is used.

Every bound is rational. Exp may have relative error 2^-20; sqrt may have
relative error 2^-20. Scalar operations are gradual-underflow RN32 with any
of the existing local no-reassociation contractions. Target object selection
and the observer seed/transparent-guard premises are composed separately.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as W
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as BC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as CFG
from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as STATS
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL_CFG
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as DEPLOY

QUALIFICATION = 'OU3_ALT_CONFIGURED_FRONTEND_UNIFORM_FINITE_SUPPLIES_V1'
ROOT = Path(__file__).resolve().parents[3]
SOURCE_HASHES = {
    'src/tuner/SeaStateFusionTunerCommon.h': '218c6db9a9305b1d27c320d6471d4e5d3324b65ddfd33cdcd73d6859c0cb91f3',
    'src/tuner/AdaptiveWaveBandPass.h': 'ddc8c59444008552de60f8ebe2298b0db08c8bc28aa9fff679ef1772a0ad5021',
    'src/tuner/SeaStateAutoTuner.h': '00eab99b6806323836277b7986047a51c7631945f0494b6d95c2999edd49172e',
    'src/tuner/SeaStateAdaptationLimits.h': '91129a09c6fe328380410fee927936eeab19d695428328d7825aed06ca082abb',
}
DT = W.DT
G = B.rn32(F(980665, 100000))
BENCH_SIGMA = B.rn32(F(12, 100))
VERTICAL = F(512)
LPF = F(1024)
LOWPASS = F(1024)
BAND = F(2048)
ENERGY = F(32768)
MEAN = F(4096)
SQUARE = F(8388608)
WEIGHT = F(2)
LOW_ALPHA = (F(1, 4096), F(1, 32))
HIGH_ALPHA = (F(1, 512), F(1, 4))
STATS_ALPHA = (F(1, 16384), F(1, 512))
LPF_ALPHA = (F(1, 8), F(1, 4))  # input weight, not the LPF's named decay

# Each band covariance monomial contains at most eight rounded operations,
# including rounded a10/b1 and left-associated c2. The longest sum adds three
# further roundings. 32 covers these and both permitted fused alternatives.
# An additive underflow error can be multiplied only by coefficients <=2 or
# by a covariance state <=8; summing at most32 such paths is dominated by
# 2^64*ETA. The polynomials are linear in each covariance predecessor. Stats
# have at most one state factor<=2^23 downstream of a coefficient operation,
# no repeated products of states, and fewer than32 sites: the same tail covers
# their additive errors. Signal/EMA expressions have shorter operation paths.
ROUND_FACTOR = (1 + W.U)**32
ABS_TAIL = 2**64 * W.ETA


def _round_graph(ideal_abs_sum):
    return ROUND_FACTOR * F(ideal_abs_sum) + ABS_TAIL


def default_band_config():
    return CFG.BandConfig(F(1, 2), F(4), F(1, 100), F(6), F(3, 100), F(6, 5))


def default_stats_config():
    return CFG.StatsConfig(F(4), F(3, 10), F(60), F(1, 20), F(5))


def _compiled_fields(cfg, names):
    return tuple(B.rn32(getattr(cfg, name)) for name in names)


def require_config(*, band_cfg, stats_cfg, still_cfg, cutoff_hz, dt,
                   bench_noise_sigma):
    """Bind the induction to actual compiled settings, not arbitrary setters."""
    band_names = ('low_ratio', 'high_ratio', 'min_hz', 'max_hz',
                  'tune_freq_floor', 'tune_freq_ceil', 'tune_freq_prior')
    stats_names = ('K_periods', 'tau_var_min', 'tau_var_max', 'f_min', 'f_max')
    still_names = ('gravity', 'energy_alpha', 'energy_thresh')
    for actual, expected, names in ((band_cfg, default_band_config(), band_names),
            (stats_cfg, default_stats_config(), stats_names),
            (still_cfg, STILL_CFG.Config(), still_names)):
        if _compiled_fields(actual, names) != _compiled_fields(expected, names):
            raise ValueError('frontend uniform supplies require configured shipping constants')
    if (F(cutoff_hz), F(dt), F(bench_noise_sigma)) != (F(6), DT, BENCH_SIGMA):
        raise ValueError('frontend cutoff/dt/bench-noise detached from configured induction')


def _alpha_interval(xlo, xhi):
    # Monotonicity covers every intermediate rounded argument, not samples.
    elo = W.exp_result_interval(-xhi, W.ERROR_PROFILE)[0]
    ehi = W.exp_result_interval(-xlo, W.ERROR_PROFILE)[1]
    assert F(1, 2) < elo <= ehi < 1
    # Sterbenz: 1-e and 1-(1-e) are exact for these binary32 e>=1/2.
    # The second subtraction recovers representable e exactly.
    return 1-ehi, 1-elo


@lru_cache(maxsize=1)
def _build():
    # WPE's source gate pins the complete deployment header and constants.
    W.build(4, W.ERROR_PROFILE)
    for path, expected in SOURCE_HASHES.items():
        if sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise ValueError('frontend source changed: ' + path)
    assert BC._source_shape_matches() and STATS._source_shape_matches()
    assert DEPLOY._source_shape_matches()
    input_proof = W.physical_input_certificate(W.ERROR_PROFILE)
    assert input_proof['vertical_abs_upper'] < VERTICAL
    from tools.stability.ou3_alt_contraction import target_wpe_libm as TARGET
    target_exp = TARGET.certificate()['exp']
    target_profile = TARGET.profile_correspondence()
    assert target_exp['argument_interval'][0] <= -60
    assert target_exp['argument_interval'][1] >= 0
    assert target_exp['total_relative_exp_error'] < W.EXP_SQRT_REL_ERROR

    f_lo, f_hi = DEPLOY.MIN_FREQ, DEPLOY.MAX_FREQ
    # Source corner clamps are inactive on these default extrema. All corner
    # operations are monotone, so the endpoint intervals cover lagged/current
    # frequency mismatch and every possible switching history.
    low = tuple(B.mul(B.rn32(F(1, 2)), f) for f in (f_lo, f_hi))
    high = tuple(B.mul(F(4), f) for f in (f_lo, f_hi))
    upper = min(F(6), B.div(BC.NYQUIST_FACTOR, DT))
    assert B.rn32(F(1, 100)) < low[0] <= low[1] < B.div(upper, BC.SPACING)
    assert high[1] < upper and all(B.mul(BC.SPACING, l) < h for l, h in zip(low, high))
    arg = lambda f: B.mul(B.mul(BC.TWO_PI, f), DT)
    actual_alpha = {
        'LPF': _alpha_interval(arg(F(6)), arg(F(6))),
        'band_low': _alpha_interval(arg(low[0]), arg(low[1])),
        'band_high': _alpha_interval(arg(high[0]), arg(high[1])),
        # The source clamps sea_time to[.5,6], hence 4*(2*sea)>=4;
        # clampDynamicEmaHorizonSec limits the actual horizon to35.
        'stats': _alpha_interval(B.div(DT, F(35)), B.div(DT, F(4))),
    }
    boxes = {'LPF': LPF_ALPHA, 'band_low': LOW_ALPHA,
             'band_high': HIGH_ALPHA, 'stats': STATS_ALPHA}
    for name, (lo, hi) in actual_alpha.items():
        assert boxes[name][0] <= lo <= hi <= boxes[name][1]

    al, ah, av = LOW_ALPHA[0], HIGH_ALPHA[0], STATS_ALPHA[0]
    lp = LPF_ALPHA[0]
    after = {
        'LPF': _round_graph((1-lp)*LPF + lp*VERTICAL),
        'band_lowpass': _round_graph((1-al)*LOWPASS + al*VERTICAL),
    }
    hp = W.rounded_abs(W.rounded_abs(VERTICAL + LOWPASS))  # q_low<=1
    assert hp < BAND
    after['band'] = _round_graph((1-ah)*BAND + ah*hp)

    # Same-history q_l=1-a_l, q_h=1-a_h is essential. These are upper bounds
    # on the absolute monomial sums, before the generic graph roundoff factor.
    # p00: 2q_l^2+a_l^2 = 2-4a_l+3a_l^2, decreasing on the chosen interval.
    # p01: q_l(2a_h*q_l+2q_h)+a_l*a_h*q_l
    #      = q_l(2-a_l*a_h) <=2q_l.
    # p11: 3a_h^2*q_l^2+4a_h*q_l*q_h+8q_h^2
    #      <=8-12a_h+7a_h^2, decreasing for a_h<=1/4.
    after['p00'] = _round_graph(2-4*al+3*al*al)
    after['p01_abs'] = _round_graph(2*(1-al))
    after['p11'] = _round_graph(8-12*ah+7*ah*ah)

    an = W.rounded_abs(LPF/G)
    inst = W.rounded_abs(an*an)
    ea = B.rn32(F(1, 20)); ed = B.sub(F(1), ea)
    after['energy'] = _round_graph(ed*ENERGY + ea*inst)
    x2 = W.rounded_abs(BAND*BAND)
    after['mean'] = _round_graph((1-av)*MEAN + av*BAND)
    after['square'] = _round_graph((1-av)*SQUARE + av*x2)
    after['weight'] = _round_graph((1-av)*WEIGHT + av)
    limits = {'LPF': LPF, 'band_lowpass': LOWPASS, 'band': BAND,
              'p00': F(2), 'p01_abs': F(2), 'p11': F(8), 'energy': ENERGY,
              'mean': MEAN, 'square': SQUARE, 'weight': WEIGHT}
    margins = {name: limits[name]-value for name, value in after.items()}
    assert all(margin > 0 for margin in margins.values())

    # Variance is evaluated only if BOTH actual weights exceed READY_WEIGHT.
    # Squaring the debiased mean is eager even if the result will be floored.
    mu = W.rounded_abs(MEAN/STATS.READY_WEIGHT)
    mu2 = W.rounded_abs(mu*mu)
    second = W.rounded_abs(SQUARE/STATS.READY_WEIGHT)
    variance_intermediate = W.rounded_abs(second+mu2)
    variance = W.rounded_abs(second)  # subtract nonnegative square, then floor
    noise_root = W.rounded_abs((1+W.EXP_SQRT_REL_ERROR)*F(3))  # sqrt8<3
    noise = W.rounded_abs(BENCH_SIGMA*noise_root)
    assert noise < F(1, 2)
    var_wave = W.rounded_abs(variance)
    # sqrt is monotone; a power-of-two rational enclosure avoids any diagnostic
    # floating arithmetic being mistaken for part of the proof.
    sigma_root_cap = F(2**22)
    assert sigma_root_cap**2 > var_wave
    sigma_unclamped = W.rounded_abs(DEPLOY.SIGMA_COEFF *
        W.rounded_abs((1+W.EXP_SQRT_REL_ERROR)*sigma_root_cap))
    assert max(DEPLOY.MAX_SIGMA, noise, B.rn32(F(1, 20))) == 4
    largest = max(variance_intermediate, sigma_unclamped, SQUARE, ENERGY)
    max_finite = F((2**24-1)*2**104)
    assert largest < max_finite
    return {
        'qualification': QUALIFICATION,
        'source_hashes': dict(SOURCE_HASHES),
        'input_qualification': input_proof['input_qualification'],
        'vertical_input_abs_upper': VERTICAL,
        'actual_alpha_intervals': actual_alpha,
        'inductive_state_bounds': limits,
        'inductive_after_upper': after,
        'strict_inductive_margins': margins,
        'highpass_abs_upper': hp,
        'stillness_instant_energy_upper': inst,
        'still_time_interval': (F(0), F(60)),
        'attenuation_interval': (F(0), F(1)),
        'variance_upper': variance,
        'band_noise_sigma_upper': noise,
        'sigma_target_upper': F(4),
        'largest_ordinary_intermediate_abs_upper': largest,
        'binary32_overflow_margin': max_finite-largest,
        'exp_argument_interval': (F(-60), F(0)),
        'source_uniform_configured_frontend_arithmetic_totality_closed': True,
        'all_finite_prefixes_from_literal_frontend_reset_covered': True,
        'startup_deadline_required': False,
        'same_history_band_noise_covariance_and_stats_input_retained': True,
        'ordinary_rounding_and_local_FMA_alternatives_included': True,
        'exp_relative_error_allowance': W.EXP_SQRT_REL_ERROR,
        'sqrt_relative_error_allowance': W.EXP_SQRT_REL_ERROR,
        'target_exp_all_configured_frontend_arguments_covered': True,
        'configured_frontend_exp_runtime_accepts_target_error_profile': True,
        'pinned_frontend_exp_approximation_correspondence_closed': target_profile[
            'pinned_WPE_exp_log_approximation_correspondence_closed'],
        'startup_Racc_dormant_guard_RAOne_stage_branch_is_identity_or_nominal_restore': True,
        'Racc_nominal_positive_finite_source_required': True,
        'startup_capture_and_complete_MEKF_totality_proved_here': False,
        'target_frontend_compiler_and_runtime_error_profile_attachment_closed': False,
        'storage_search_allowed': False,
    }


def build():
    return deepcopy(_build())


def validate(document):
    if document != build():
        raise ValueError('frontend finite-supply certificate changed')


def require_state(vertical, frontend):
    """Check retained states without replacing them by a zero/reset witness."""
    from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
    from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as FS
    if not isinstance(vertical, VS.State) or not isinstance(frontend, FS.State):
        raise TypeError('same-history vertical and frontend machine states required')
    if vertical.samples != frontend.samples:
        raise ValueError('vertical and sigma frontend sample ordinals detached')
    if vertical.lpf.cutoff_hz != 6 or abs(vertical.lpf.value) > LPF:
        raise ValueError('LPF state outside configured uniform supply')
    still = vertical.stillness
    if not 0 <= still.energy <= ENERGY or not 0 <= still.still_time <= 60:
        raise ValueError('stillness state outside uniform supply')
    band = frontend.band.machine
    if not (abs(band.lowpass_low) <= LOWPASS and abs(band.band) <= BAND
            and 0 <= band.p00 <= 2 and abs(band.p01) <= 2 and 0 <= band.p11 <= 8):
        raise ValueError('carried band signal/covariance outside uniform supply')
    stats = frontend.stats
    if not (abs(stats.mean_value) <= MEAN and 0 <= stats.sq_value <= SQUARE
            and 0 <= stats.mean_weight <= WEIGHT and 0 <= stats.sq_weight <= WEIGHT):
        raise ValueError('carried statistics moments outside uniform supply')


def require_event(vertical_result, frontend_result, *, sigma_target=None):
    """Attach the induction to one already source-executed complete event.

    The surrounding source join owns raw admission and configuration. This
    binder preserves all predecessor/successor objects and checks the SAME
    current vertical, current band output, noise covariance and variance.
    """
    from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
    from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as FS
    if not isinstance(vertical_result, VS.Result) or not isinstance(frontend_result, FS.Result):
        raise TypeError('source-executed vertical and frontend results required')
    proof = build()
    require_state(vertical_result.before, frontend_result.before)
    require_state(vertical_result.state, frontend_result.state)
    if abs(vertical_result.band_input) > VERTICAL:
        raise ValueError('frontend input detached from physical Mahony envelope')
    bc = frontend_result.band.coefficients
    if not bc.active or bc.dt != DT or not DEPLOY.MIN_FREQ <= bc.f_ref <= DEPLOY.MAX_FREQ:
        raise ValueError('frontend band coefficients detached from configured active domain')
    if frontend_result.band.envelope.x != vertical_result.band_input:
        raise ValueError('band did not consume SAME new machine Mahony vertical')
    for alpha, box in ((bc.alpha_low, LOW_ALPHA), (bc.alpha_high, HIGH_ALPHA),
                      (frontend_result.stats_envelope.coefficients.alpha, STATS_ALPHA)):
        if not box[0] <= alpha <= box[1]:
            raise ValueError('frontend coefficient outside all-input library enclosure')
    if not (0 <= frontend_result.accel_variance <= proof['variance_upper']
            and 0 <= frontend_result.band_noise_sigma <= proof['band_noise_sigma_upper']):
        raise ValueError('source-derived variance/noise outside uniform supply')
    if sigma_target is not None:
        FS.require_sigma(frontend_result, sigma_target)
        VS.require_sigma_stillness(vertical_result, sigma_target)
        if sigma_target.cfg.sigma_coeff != DEPLOY.SIGMA_COEFF or not sigma_target.cfg.clamp_enabled:
            raise ValueError('sigma target detached from configured clamp/coefficient')
        if sigma_target.cfg.max_sigma != 4 or not 0 <= sigma_target.sigma_target <= 4:
            raise ValueError('sigma target outside configured uniform supply')
    return proof
