"""Inductive deployed-tau domain through bounded startup plus the 600-step word.

Live does not in general begin at the construction seed 1.1f: once Cold ends,
the measurement-only tuner may adapt during startup.  This theorem therefore
starts at the actual construction seed and covers the *entire* possible bounded
startup followed by the canonical 600 Live IMU transitions.

The wrapper is designed for 200 Hz (5 ms), while the bootstrap timeout is 150 s.
Using 30,000 possible pre-Live updates is conservative because Cold returns
before adaptation for its first 5 s.  Add the 600 Live updates for 30,600 total.

The source-locked target remains in [0.4,12] (the true lower rail at f=1.2 is
about 0.4167).  The tightened binary32 theorem gives at most B=2^-20 s local
deployed-minus-convex-shadow discrepancy for either compiler contraction shape.
Hence after k updates

    tau_k >= min(1.1,0.4) - k B,
    tau_k <= max(1.1,12) + k B.

At k=30,600 the envelope is still about [0.3708,12.0292], strictly positive and
strictly below the 13 s premise of the local theorem.  This closes the bounded
startup + one-word scalar domain.  It does not prove indefinite word tiling or
machine lifetime.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_tau_roundoff_bound as ROUND

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
PRELIVE_MAX_UPDATES=30000
WORD_STEPS=600
TOTAL_UPDATES=PRELIVE_MAX_UPDATES+WORD_STEPS
INITIAL=B32.rn32(F(11,10)); TARGET_LOWER=F(2,5); TARGET_UPPER=F(12); DOMAIN_ABS_MAX=F(13)


def lower_after(k:int):
    if not isinstance(k,int) or not 0<=k<=TOTAL_UPDATES: raise ValueError('k outside bounded startup plus canonical word')
    return min(INITIAL,TARGET_LOWER)-k*ROUND.UNIFORM_RESIDUAL_MAX


def upper_after(k:int):
    if not isinstance(k,int) or not 0<=k<=TOTAL_UPDATES: raise ValueError('k outside bounded startup plus canonical word')
    return max(INITIAL,TARGET_UPPER)+k*ROUND.UNIFORM_RESIDUAL_MAX


def assert_induction_closes():
    hi=TARGET.evaluate(TARGET.CEIL)
    if hi.exact_target<TARGET_LOWER or hi.binary32_target<TARGET_LOWER:
        raise AssertionError('shipping maximum-frequency tau target below 0.4 s rail')
    if TARGET.TAU_MAX!=TARGET_UPPER: raise AssertionError('shipping upper tau target detached from 12 s rail')
    if not lower_after(TOTAL_UPDATES)>0: raise AssertionError('startup+word lower tau induction lost positivity')
    if not upper_after(TOTAL_UPDATES)<DOMAIN_ABS_MAX: raise AssertionError('startup+word upper tau induction escaped local roundoff domain')
    return True


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(n in s for n in (
      'float tau_applied   = 1.1f;',
      'constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;',
      'float proxy_startup_timeout_sec = 150.0f;',
      'constexpr float MAX_TUNE_FREQ_HZ = 1.2f;',
      'constexpr float MAX_TAU_S   = 12.0f;',
      'float tau_coeff_    = 1.0f;'))


def readiness():
    closed=bool(_source_shape_matches() and TARGET.readiness()['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'] and assert_induction_closes())
    return {
      'conservative_prelive_tau_update_cap':PRELIVE_MAX_UPDATES,
      'canonical_word_tau_updates':WORD_STEPS,
      'shipping_tau_initial_timeout_rate_and_target_domain_source_shape_matches':_source_shape_matches(),
      'tau_target_lower_rail_from_shipping_frequency_ceiling':True,
      'both_FMA_and_separate_roundoff_shapes_covered':True,
      'shipping_tau_predecessor_domain_inductively_closed_through_startup_and_600_word':closed,
      'source_uniform_tau_roundoff_supply_bound_closed_through_startup_and_600_word':closed,
      'successive_words_tau_domain_tiled_indefinitely':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
