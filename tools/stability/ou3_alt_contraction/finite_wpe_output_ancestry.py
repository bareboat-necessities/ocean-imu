"""Bind canonical WPE output witnesses to the SAME log-period state.

``finite_wpe_runtime`` intentionally leaves exp/log/sqrt deployment arithmetic
open, but its legacy current/post output arguments can otherwise be supplied as
an unrelated reciprocal pair.  This wrapper removes that history-splicing
freedom: every visible period/frequency pair carries the exact log-period state
from which shipping evaluates ``exp(log_period)`` and ``exp(-log_period)``.

This is an ancestry theorem only.  It does NOT claim that the supplied rational
period/frequency are the target libm's binary32 exp results.  That numerical
correspondence remains an explicit deployment blocker.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
QUALIFICATION='OU3_ALT_WPE_CANONICAL_OUTPUT_ANCESTRY_V1'


@dataclass(frozen=True)
class CanonicalOutput:
    log_period:F
    period:F
    frequency:F
    def __post_init__(self):
        lp,p,f=map(P.rational,(self.log_period,self.period,self.frequency))
        if p<=0 or f<=0 or p*f!=1:
            raise ValueError('exact-real WPE shadow output must be one positive reciprocal pair')
        object.__setattr__(self,'log_period',lp)
        object.__setattr__(self,'period',p)
        object.__setattr__(self,'frequency',f)


def update(state:W.WPEState,cfg:W.WPEConfig,*,current:CanonicalOutput|None=None,
           post:CanonicalOutput|None=None,**kwargs):
    """Run one WPE edge with current/post outputs owned by its log state."""
    if not isinstance(state,W.WPEState) or not isinstance(cfg,W.WPEConfig):
        raise TypeError('WPE state/config required')
    if state.log_period is None:
        if current is not None: raise ValueError('no current canonical output before log-period state exists')
        cp=cf=None
    else:
        if not isinstance(current,CanonicalOutput):
            raise TypeError('existing log-period state requires source-bound current output')
        if current.log_period != state.log_period:
            raise ValueError('current canonical output detached from carried log-period state')
        cp,cf=current.period,current.frequency
    if post is None:
        post_shadow=None
    else:
        if not isinstance(post,CanonicalOutput): raise TypeError('source-bound post output required')
        post_shadow=W.CanonicalOutputWitness(post.period,post.frequency)
    out=W.update(state,cfg,current_period=cp,current_frequency=cf,
                 post_output=post_shadow,**kwargs)
    if out.produced_period:
        if post is None: raise AssertionError('produced WPE period without source-bound post output')
        if out.state.log_period != post.log_period:
            raise ValueError('post canonical output detached from computed log-period successor')
        if out.period != post.period or out.frequency != post.frequency:
            raise AssertionError('lower WPE output changed source-bound canonical pair')
    elif post is not None:
        raise ValueError('non-producing WPE branch consumes no post canonical output')
    return out


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(n in s for n in (
      'return std::isfinite(log_period_sec_) ? std::exp(log_period_sec_) : NAN;',
      'return std::isfinite(log_period_sec_) ? std::exp(-log_period_sec_) : NAN;'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_period_frequency_share_one_log_state_source_shape':_source_shape_matches(),
      'current_output_bound_to_carried_log_period_state':True,
      'post_output_bound_to_computed_log_period_successor':True,
      'cross_history_reciprocal_output_splicing_forbidden':True,
      'period_exp_binary32_libm_correspondence_closed':False,
      'frequency_exp_binary32_libm_correspondence_closed':False,
      'upstream_log_sqrt_binary32_correspondence_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
