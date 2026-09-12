"""Inductive deployed-tau domain for the canonical 600-transition Live word.

The local tau roundoff theorem needs |tau_applied| <= 13 s.  This module proves
that premise for the whole finite word without traces or sampling.

Shipping starts at 1.1f.  With f_tune clamped to <=1.2f and tau_coeff=1.0f, the
exact machine-operand target is at least 0.5/1.2 > 0.4 s and at most 12 s.  The
canonical exact shadow update is a convex combination of the current deployed
tau and that target because 0 <= 1-e_f <= 1.  The deployed successor differs
from that convex shadow by at most B=2^-14 s for either compiler contraction
shape.

Thus after k <= 600 updates

    tau_k >= min(tau_0, 0.4) - k B,
    tau_k <= max(tau_0, 12) + k B.

At k=600 these are >0.36 s and <12.037 s, respectively, so every predecessor
remains positive and strictly below the 13 s premise needed by the local bound.
This is a finite-word induction only; it does not solve indefinite machine
lifetime or successive-word tiling.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_tau_roundoff_bound as ROUND

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
WORD_STEPS=600
INITIAL=B32.rn32(F(11,10))
TARGET_LOWER=F(2,5)
TARGET_UPPER=F(12)
DOMAIN_ABS_MAX=F(13)


def lower_after(k:int):
    if not isinstance(k,int) or not 0<=k<=WORD_STEPS: raise ValueError('k outside canonical 600-step word')
    return min(INITIAL,TARGET_LOWER)-k*ROUND.UNIFORM_RESIDUAL_MAX


def upper_after(k:int):
    if not isinstance(k,int) or not 0<=k<=WORD_STEPS: raise ValueError('k outside canonical 600-step word')
    return max(INITIAL,TARGET_UPPER)+k*ROUND.UNIFORM_RESIDUAL_MAX


def assert_induction_closes():
    # Source-target theorem at the largest allowed frequency establishes the
    # lower target rail from actual compiled binary32 operands.
    hi=TARGET.evaluate(TARGET.CEIL)
    if hi.exact_target < TARGET_LOWER or hi.binary32_target < TARGET_LOWER:
        raise AssertionError('shipping maximum-frequency tau target below 0.4 s rail')
    if TARGET.TAU_MAX != TARGET_UPPER:
        raise AssertionError('shipping upper tau target detached from 12 s rail')
    if not lower_after(WORD_STEPS)>0:
        raise AssertionError('600-step lower tau induction lost positivity')
    if not upper_after(WORD_STEPS)<DOMAIN_ABS_MAX:
        raise AssertionError('600-step upper tau induction escaped local roundoff domain')
    return True


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(n in s for n in (
      'float tau_applied   = 1.1f;',
      'constexpr float MAX_TUNE_FREQ_HZ = 1.2f;',
      'constexpr float MAX_TAU_S   = 12.0f;',
      'float tau_coeff_    = 1.0f;'))


def readiness():
    closed=bool(_source_shape_matches() and TARGET.readiness()['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'] and assert_induction_closes())
    return {
      'canonical_word_tau_updates':WORD_STEPS,
      'shipping_tau_initial_and_target_domain_source_shape_matches':_source_shape_matches(),
      'tau_target_lower_rail_from_shipping_frequency_ceiling':True,
      'both_FMA_and_separate_roundoff_shapes_covered':True,
      'shipping_tau_predecessor_domain_inductively_closed_for_600_step_word':closed,
      'source_uniform_tau_roundoff_supply_bound_closed_for_600_step_word':closed,
      'successive_words_tau_domain_tiled_indefinitely':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
