"""Exact binary32 SeaStateAutoTuner frequency storage for ALT.

Shipping ``SeaStateAutoTuner::update`` does not smooth frequency.  For every
accepted update it copies the supplied float through the literal frequency
clamp and stores that float in ``frequency_hz``; ``getFrequencyHz`` returns the
stored value unchanged.  This module closes that storage/consumer edge without
claiming that the upstream WavePeriodEstimator has yet been reduced to a
binary32 deployment graph.

The theorem-facing use is therefore:

    WPE binary32 output (still open) -> store() -> stored frequency
      -> getFrequencyHz() identity -> tau target / EMA graph.

NaN/Inf are excluded here by requiring a finite normal binary32 input; shipping
rejects non-finite inputs before the assignment.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/SeaStateAutoTuner.h'
QUALIFICATION='OU3_ALT_TUNER_FREQUENCY_BINARY32_STORE_V1'


@dataclass(frozen=True)
class StoredFrequency:
    input_hz:F
    min_hz:F
    max_hz:F
    stored_hz:F
    def __post_init__(self):
        vals=tuple(F(x) for x in (self.input_hz,self.min_hz,self.max_hz,self.stored_hz))
        if not all(B.is_binary32(x) for x in vals):
            raise ValueError('tuner frequency store operands must be actual binary32 values')
        i,lo,hi,out=vals
        if i<=0 or lo<=0 or hi<lo:
            raise ValueError('positive finite tuner frequency domain required')
        expected=max(lo,min(hi,i))
        if out!=expected:
            raise ValueError('stored tuner frequency detached from shipping clamp')
        object.__setattr__(self,'input_hz',i); object.__setattr__(self,'min_hz',lo)
        object.__setattr__(self,'max_hz',hi); object.__setattr__(self,'stored_hz',out)


def store(input_hz, min_hz, max_hz):
    """One accepted ``SeaStateAutoTuner::update`` frequency assignment."""
    i,lo,hi=map(F,(input_hz,min_hz,max_hz))
    if not all(B.is_binary32(x) for x in (i,lo,hi)):
        raise ValueError('accepted shipping frequency operands must already be binary32')
    return StoredFrequency(i,lo,hi,max(lo,min(hi,i)))


def get_frequency_hz(state:StoredFrequency):
    if not isinstance(state,StoredFrequency): raise TypeError('StoredFrequency required')
    return state.stored_hz


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'float f_eff = f_input_hz;',
      'f_eff = std::max(f_min_hz, std::min(f_max_hz, f_eff));',
      'frequency_hz = f_eff;',
      'inline float getFrequencyHz() const { return frequency_hz; }',
    )
    return all(n in s for n in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_frequency_store_source_shape_matches':_source_shape_matches(),
      'accepted_input_must_be_actual_binary32':True,
      'frequency_clamp_and_store_exact_binary32':True,
      'getFrequencyHz_is_identity_on_stored_binary32':True,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
