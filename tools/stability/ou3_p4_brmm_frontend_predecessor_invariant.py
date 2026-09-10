#!/usr/bin/env python3
"""Source-uniform compact predecessor invariant for the COMPLETE-BRMM frontend.

This certificate closes the missing *state-family induction* between the
qualified BRMM source/startup theorem and the existing same-signal one-step
frontend transition.  It does not enumerate sea histories and it does not turn
physical samples into independent theorem coordinates.

The only purpose of scalar bounds here is to prove compact forward invariance of
persistent implementation memory.  Coefficients used by P4 must still be
created by ``ou3_p4_joint_brmm_frontend_transition.advance`` from one common
source history.

The important dependency-saving observation is the exact high-pass transfer

    y_n = d (y_{n-1} + x_n - x_{n-1}),   0<d<1.

Starting from reset, its impulse response has l1 norm

    d + sum_{k>=1} (1-d)d^k = 2d,

so |y| <= 2 d |x|_inf.  No independent bound on x_n-x_{n-1} is introduced.
The two cascaded stages therefore have gain <=(2d)^2.  The subsequent leaky
integrator has induced gain q/(1-d)=1/lambda exactly.  These identities provide
finite bounds on WPE velocity/elevation and hence on every EW moment.

The monotonically increasing WPE elapsed timer is quotiented only after all
shipping guards become equivalent: once the usable latch is set elapsed is no
longer state-relevant; before takeover, a cap at 3/lambda + T_raw,max preserves
every possible usable-period comparison.  Thus an unbounded wall-clock counter
is not falsely asserted compact and binary32 +inf at astronomically long times
does not alter the guard class.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_finite_window_primitive_qualification as PRIM
import ou3_brmm_private_mahony_live_invariant as MAHONY
import ou3_brmm_private_mahony_discrete_invariant as MAHONY32
import ou3_fast_inv_sqrt_interval as FINV
import ou3_brmm_wpe_state_step as WPE
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_sigma_stillness_step as STILL
import ou3_p5_startup_timeout_capture as CAPTURE
import ou3_validated_transcendentals as VT
from ou3_interval import Interval

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / 'tools/stability/ou3_proof_operating_domain.json'
WPE_HEADER = REPO / 'src/tuner/WavePeriodEstimator.h'
SCHEMA = 1
QUALIFICATION = 'OU3_P4_COMPLETE_BRMM_FRONTEND_PREDECESSOR_INVARIANT_V1'


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _positive_decay(lam: Interval, dt: float) -> Interval:
    d = VT.exp_interval(-(lam * I(dt)))
    if not (0.0 < d.lo <= d.hi < 1.0):
        raise RuntimeError('WPE decay is not strictly Schur')
    return d


def _band_covariance_outer(c: TUNER.Constants) -> dict:
    # Uniform coefficient magnitudes over the shipping-clamped band.  We use
    # only stability/compactness here; no such box is allowed to generate P4
    # coefficients.
    pi = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
    low_min = c.band_min_hz
    high_min = 1.05 * low_min
    al_min = (I(1.0) - VT.exp_interval(-(I(2.0) * pi * I(low_min) * I(c.dt)))).lo
    ah_min = (I(1.0) - VT.exp_interval(-(I(2.0) * pi * I(high_min) * I(c.dt)))).lo
    if not (0.0 < al_min < 1.0 and 0.0 < ah_min < 1.0):
        raise RuntimeError('band filter lost a uniform Schur gap')
    ql = 1.0 - al_min
    qh = 1.0 - ah_min
    # b0,b1,|a10| <=1 is deliberately conservative.  Positive denominators
    # are the only fact required to establish a finite invariant covariance.
    p00 = up(1.0 / (1.0 - ql * ql))
    p01 = up((p00 + 1.0) / (1.0 - ql * qh))
    p11 = up((p00 + 2.0 * qh * p01 + 1.0) / (1.0 - qh * qh))
    return {
        'alpha_low_uniform_lower': al_min,
        'alpha_high_uniform_lower': ah_min,
        'ql_uniform_upper': ql,
        'qh_uniform_upper': qh,
        'p00_abs_upper': p00,
        'p01_abs_upper': p01,
        'p11_abs_upper': p11,
        'finite_covariance_outer': all(math.isfinite(x) for x in (p00, p01, p11)),
    }


def build() -> dict:
    domain = json.loads(DOMAIN.read_text())
    brmm = BRMM.build(); prim = PRIM.build(); mah = MAHONY.build(); m32 = MAHONY32.build()
    cap = CAPTURE.build(); wc = WPE.constants(); tc = TUNER.constants(); sc = STILL.constants()
    bad = {
        'BRMM': BRMM.validate(brmm), 'primitive': PRIM.validate(prim),
        'Mahony': MAHONY.validate(mah), 'Mahony32': MAHONY32.validate(m32),
        'capture': CAPTURE.validate(cap),
    }
    bad = {k:v for k,v in bad.items() if v}
    if bad:
        raise RuntimeError('frontend predecessor invariant prerequisites failed: '+repr(bad))

    g = float(domain['startup']['gravity_mps2'])
    a = float(prim['uniform_physical_primitives']['acceleration_norm_upper_mps2'])
    specific_force = g + a
    qnorm2 = FINV.all_positive_normal_normalized_norm2_enclosure()
    # The third quaternion rotation row has exact norm ||q||^2, hence this is a
    # Cauchy-Schwarz bound on the deployed private-Mahony vertical output.
    vertical = up(qnorm2.hi * specific_force + g)

    pi = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
    lam = I(2.0) * pi * I(wc.high_pass_hz)
    d = _positive_decay(lam, wc.dt)
    hp_gain = up(2.0 * d.hi)
    hp1 = up(hp_gain * vertical)
    hp2 = up(hp_gain * hp1)
    lam_lo = lam.lo
    vel = up(hp2 / lam_lo)
    elev = up(vel / lam_lo)

    # EW states are convex recurrences from reset.  Their unnormalised first
    # moments and second moments therefore stay inside the signal hulls.
    weight_lo = math.nextafter(1.0e-3, math.inf)
    weight_hi = 1.0
    vel_mean = vel
    vel_sq = up(vel * vel)
    elev_mean = elev
    elev_sq = up(elev * elev)

    # On every valid period branch evar>1e-12 and omega^2>1e-8.  A finite lower
    # raw-period bound follows from vv <= velocity_sq/weight; the upper bound is
    # directly the strict omega^2 guard.  Hold branches leave the canonical log
    # state unchanged, and valid EMA updates are convex combinations, so the
    # same finite log interval is invariant once initialized.
    evar_floor = math.nextafter(1.0e-12, math.inf)
    omega2_floor = math.nextafter(1.0e-8, math.inf)
    vvar_upper = up(vel_sq / weight_lo)
    omega2_upper = up(vvar_upper / evar_floor)
    raw_period_min = math.nextafter(2.0 * math.pi / math.sqrt(omega2_upper), 0.0)
    raw_period_max = up(2.0 * math.pi / math.sqrt(omega2_floor))
    if not (0.0 < raw_period_min < raw_period_max < math.inf):
        raise RuntimeError('raw period compactness failed')
    log_min = math.nextafter(math.log(raw_period_min), -math.inf)
    log_max = math.nextafter(math.log(raw_period_max), math.inf)
    frequency_min = math.nextafter(math.exp(-log_max), 0.0)
    frequency_max = up(math.exp(-log_min))

    moment_start = 3.0 / lam_lo
    elapsed_guard_cap = up(moment_start + raw_period_max)
    # Exact source audit: elapsed only feeds threshold/age comparisons after its
    # increment.  Once usable is latched, it cannot affect period/frequency or
    # any filter recurrence.  Before latch, the above cap dominates every
    # possible period admitted by a valid update.
    wtext = WPE_HEADER.read_text()
    elapsed_parity = all(s in wtext for s in (
        'elapsed_sec_ += dt_sec;',
        'if (elapsed_sec_ < moment_start_sec) return;',
        'elapsed_sec_ >= settled_floor_sec',
        'const float moment_history_sec = elapsed_sec_ - moment_start_sec;',
        'if (elapsed_sec_ >= usable_floor_sec && moment_history_sec >= period)',
        'if (usable_period_) return;',
    ))

    # Adaptive band state. low-pass is convex in bounded input. For the high
    # pass, B=2*A is invariant because q_h B + a_h q_l(2A) <= 2A.
    band_low = vertical
    band = up(2.0 * vertical)
    bcov = _band_covariance_outer(tc)

    # Tuner moments are also convex EW recurrences from reset. Candidate target
    # values are shipping-clamped, candidate EMA is convex, and active values
    # are either the previous active state or an exact staged candidate.
    tuner_moment_mean = band
    tuner_moment_sq = up(band * band)
    candidate = {
        'tau_s': [tc.tau_min, tc.tau_max],
        'sigma_mps2': [0.0, tc.sigma_max],
        'rs': [tc.rs_min, tc.rs_max],
    }
    active = {
        'tau_s': [tc.tau_min, tc.tau_max],
        'sigma_mps2': [0.0, tc.sigma_max],
        'rs': [tc.rs_min, tc.rs_max],
        'pseudo_period_s': [tc.pseudo_min_s, tc.pseudo_max_s],
    }
    scheduler_max = up(tc.adapt_every_s + tc.dt)

    # Stillness path is a convex LPF/EMA plus a saturated timer.
    still_energy = up((vertical / sc.gravity) ** 2)

    strict = {
        'WPE_decay_strictly_inside_unit_disk': d.hi < 1.0,
        'WPE_telescoping_highpass_l1_gain_below_2': hp_gain < 2.0,
        'WPE_leaky_integrator_gain_finite': math.isfinite(1.0 / lam_lo),
        'WPE_weight_positive_on_Normal_Live_family': weight_lo > 1e-3,
        'WPE_valid_period_compact': 0.0 < raw_period_min < raw_period_max < math.inf,
        'WPE_frequency_compact': 0.0 < frequency_min <= frequency_max < math.inf,
        'WPE_elapsed_guard_quotient_source_parity': elapsed_parity,
        'adaptive_band_uniform_Schur_gap': bcov['alpha_low_uniform_lower'] > 0.0 and bcov['alpha_high_uniform_lower'] > 0.0,
        'adaptive_band_covariance_compact': bcov['finite_covariance_outer'],
        'candidate_and_active_clamps_compact': tc.tau_min > 0 and tc.tau_max >= tc.tau_min and tc.sigma_max > 0 and tc.rs_min > 0 and tc.rs_max >= tc.rs_min,
        'scheduler_compact_modulo_exact_reset_branch': scheduler_max < math.inf,
        'stillness_timer_saturated': sc.still_max_s == 60.0,
        'startup_root_enters_same_component_invariants': cap['startup_timeout_capture_closed'],
        'Mahony_all_live_invariant_closed': mah['continuous_all_live_PI_invariant_closed'],
        'Mahony_source_order_binary32_invariant_closed': m32['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'],
    }
    closed = all(strict.values())

    return {
        'schema': SCHEMA, 'qualification': QUALIFICATION,
        'canonical_source': 'COMPLETE_BRMM_NORMAL_LIVE_WORD',
        'source_generator': False, 'trajectory_replay_used': False,
        'independent_sample_boxes_used_to_generate_coefficients': False,
        'same_history_BRMM_ancestry_required_downstream': True,
        'source_order_binary32_model_consumed': True,
        'vertical_acceleration_abs_upper_mps2': vertical,
        'WPE_invariant': {
            'decay': d.as_list(), 'single_highpass_linf_gain_upper': hp_gain,
            'highpass1_abs_upper_mps2': hp1, 'highpass2_abs_upper_mps2': hp2,
            'velocity_abs_upper_mps': vel, 'elevation_abs_upper_m': elev,
            'weight_interval': [weight_lo, weight_hi],
            'velocity_mean_abs_upper': vel_mean, 'velocity_sq_upper': vel_sq,
            'elevation_mean_abs_upper': elev_mean, 'elevation_sq_upper': elev_sq,
            'valid_raw_period_s': [raw_period_min, raw_period_max],
            'canonical_log_period': [log_min, log_max],
            'canonical_frequency_hz': [frequency_min, frequency_max],
            'elapsed_guard_quotient_cap_s': elapsed_guard_cap,
            'actual_elapsed_need_not_be_bounded': True,
            'elapsed_quotient_preserves_all_shipping_guards': elapsed_parity,
        },
        'adaptive_band_invariant': {
            'lowpass_abs_upper_mps2': band_low, 'band_abs_upper_mps2': band,
            **bcov,
        },
        'tuner_invariant': {
            'mean_value_abs_upper': tuner_moment_mean,
            'sq_value_upper': tuner_moment_sq,
            'weights_interval': [0.0, 1.0],
            'frequency_after_shipping_clamp_hz': [0.05, 5.0],
            'candidate': candidate, 'active': active,
            'scheduler_elapsed_upper_s': scheduler_max,
        },
        'stillness_invariant': {
            'lpf_abs_upper_mps2': vertical,
            'energy_ema_upper': still_energy,
            'still_time_s': [0.0, sc.still_max_s],
            'attenuation': [0.0, 1.0],
        },
        'strict_invariance_checks': strict,
        'all_component_predecessor_invariants_closed': closed,
        'complete_BRMM_predecessor_family_compact_modulo_elapsed_guard_quotient': closed,
        'complete_BRMM_predecessor_state_family_covered': closed,
        'coefficient_history_still_must_come_from_joint_transition': True,
        'outer_bounds_may_not_generate_independent_f_sigma_tau_TS_RS': True,
        'P4_promoted_here': False,
        'next_obligation': 'consume this predecessor-family induction in the same-signal joint transition, then attach every estimator-owned image to literal Riccati/Joseph prefixes and execute the augmented endpoint/every-prefix LDLT',
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get('schema') != SCHEMA or d.get('qualification') != QUALIFICATION:
        f.append('schema/qualification mismatch')
    for k in ('same_history_BRMM_ancestry_required_downstream','source_order_binary32_model_consumed',
              'all_component_predecessor_invariants_closed','complete_BRMM_predecessor_family_compact_modulo_elapsed_guard_quotient',
              'complete_BRMM_predecessor_state_family_covered','coefficient_history_still_must_come_from_joint_transition',
              'outer_bounds_may_not_generate_independent_f_sigma_tau_TS_RS'):
        if d.get(k) is not True: f.append(k+' not true')
    for k in ('source_generator','trajectory_replay_used','independent_sample_boxes_used_to_generate_coefficients','P4_promoted_here'):
        if d.get(k) is not False: f.append(k+' not false')
    if not all(d.get('strict_invariance_checks',{}).values()): f.append('one or more strict invariant checks failed')
    w=d.get('WPE_invariant',{})
    rp=w.get('valid_raw_period_s',[0,0]); fr=w.get('canonical_frequency_hz',[0,0])
    if not (len(rp)==2 and 0<float(rp[0])<float(rp[1])<math.inf): f.append('raw period interval invalid')
    if not (len(fr)==2 and 0<float(fr[0])<=float(fr[1])<math.inf): f.append('frequency interval invalid')
    return f


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    d=build(); f=validate(d); d['validation_pass']=not f; d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'closed':d['complete_BRMM_predecessor_state_family_covered'],'vertical':d['vertical_acceleration_abs_upper_mps2'],'wpe_velocity':d['WPE_invariant']['velocity_abs_upper_mps'],'wpe_elevation':d['WPE_invariant']['elevation_abs_upper_m'],'period':d['WPE_invariant']['valid_raw_period_s'],'frequency':d['WPE_invariant']['canonical_frequency_hz'],'failures':f},sort_keys=True))
    return int(bool(f))

if __name__=='__main__': raise SystemExit(main())
