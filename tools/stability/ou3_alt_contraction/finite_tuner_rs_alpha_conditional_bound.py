"""Conditional source-uniform bound for the SpectralMSE R_S alpha supply.

This lemma deliberately does NOT prove target ``std::exp`` correspondence.  It
answers the next question once a machine exp result has been qualified against
the rigorous same-argument real enclosure: how large can

    alpha_binary32 - alpha_exact

be over the entire deployed R_S smoothing-horizon domain?

For deployed slew_log=0, safe tau is in [0.5,6] and mult=1.5, so the exact
horizon is [0.75,9] s and the final [0.05,35] clamp is inactive for the canonical
5 ms sample.  Binary32 multiplication and division perturb x=dt/horizon by at
most two RNE relative errors.  The rigorous second-order exp enclosure has width
x^2/2; at the worst x this dominates the float arithmetic by orders of
magnitude.  We certify the simple rational envelope

    |alpha_binary32 - alpha_exact| <= 1/40000 = 2.5e-5

conditional on the machine exp witness satisfying the same-argument enclosure.
The bound is intentionally a little wider than the derived value so it remains
simple to propagate through the finite master word.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as J

DT=B.rn32(F(1,200))
U=F(1,1<<24)                 # binary32 RNE relative unit roundoff
SAFE_TAU_MIN=F(1,2)
MULT=F(3,2)
EXACT_HORIZON_MIN=MULT*SAFE_TAU_MIN  # 0.75 s
EXACT_X_MAX=DT/EXACT_HORIZON_MIN
# One rounded multiply in horizon and one rounded divide in x.  For positive
# normal values, each RNE step is contained in y*(1 +/- U).  Combining them:
X_RATIO_HI=(1+U)/(1-U)
MACHINE_X_MAX=EXACT_X_MAX*X_RATIO_HI
X_DISCREPANCY_MAX=EXACT_X_MAX*(X_RATIO_HI-1)
EXP_INTERVAL_WIDTH_MAX=max(EXACT_X_MAX**2,MACHINE_X_MAX**2)/2
# alpha is below 0.007 throughout this domain; one full 2^-30 allowance is far
# wider than its actual half-ulp and avoids depending on a particular exponent
# boundary in this lemma.
ALPHA_FINAL_ROUND_MAX=F(1,1<<30)
DERIVED_ABS_BOUND=X_DISCREPANCY_MAX+EXP_INTERVAL_WIDTH_MAX+ALPHA_FINAL_ROUND_MAX
CERTIFIED_ABS_BOUND=F(1,40000)
if DERIVED_ABS_BOUND > CERTIFIED_ABS_BOUND:
    raise AssertionError('simple RS alpha supply bound no longer covers deployed domain')


def validate(join:J.Join):
    """Return a qualified same-source alpha join under the uniform envelope."""
    if not isinstance(join,J.Join): raise TypeError('same-source RS alpha join required')
    if join.dt!=DT: raise ValueError('conditional uniform bound covers canonical compiled 5 ms dt only')
    # The local join itself proves the machine exp witness lies in the rigorous
    # same-rounded-argument enclosure.  Platform correspondence that supplies
    # such a witness for every execution remains external/open.
    lo,hi=join.machine_minus_exact_alpha_lo,join.machine_minus_exact_alpha_hi
    if lo < -CERTIFIED_ABS_BOUND or hi > CERTIFIED_ABS_BOUND:
        raise ValueError('local alpha supply escaped certified uniform envelope')
    return join


def readiness():
    return {
      'canonical_dt_binary32_fixed':True,
      'deployed_slew_zero_horizon_domain_reduced_to_safe_tau_interval':True,
      'binary32_horizon_multiply_and_divide_roundoff_uniformly_bounded':True,
      'second_order_exp_enclosure_width_uniformly_bounded':True,
      'conditional_abs_alpha_supply_le_2p5e_minus5':True,
      'derived_bound_strictly_inside_certified_simple_bound':DERIVED_ABS_BOUND < CERTIFIED_ABS_BOUND,
      'target_libm_exp_correspondence_closed':False,
      'unconditional_machine_execution_alpha_bound_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
