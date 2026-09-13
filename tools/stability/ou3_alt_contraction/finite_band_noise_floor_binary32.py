"""Binary32 ``band_noise_floor_sigma_()`` from the persistent machine band state.

Before the adaptive band is ready shipping returns the configured bench-noise
sigma unchanged.  Once ready it reads the actual stored ``p11`` gain, computes
``sqrt(gain)``, and multiplies by the same bench sigma.  The sqrt result is tied
to the exact root through its binary32 RNE cell; platform sqrt correspondence
remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as L
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as RNE
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_BAND_NOISE_FLOOR_BINARY32_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Result:
    band:L.State
    bench_sigma:F
    sqrt_gain:F|None
    noise_sigma:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.band,L.State): raise TypeError('persistent machine band ledger required')
        object.__setattr__(self,'bench_sigma',_q(self.bench_sigma,'bench sigma'))
        object.__setattr__(self,'noise_sigma',_q(self.noise_sigma,'band noise sigma'))
        if self.sqrt_gain is not None: object.__setattr__(self,'sqrt_gain',_q(self.sqrt_gain,'sqrt gain'))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong band-noise-floor qualification')
        if self.bench_sigma<0 or self.noise_sigma<0: raise ValueError('nonnegative band-noise quantities required')
        if not self.band.machine.ready:
            if self.sqrt_gain is not None or self.noise_sigma!=self.bench_sigma:
                raise ValueError('unready band must return bench sigma exactly and consume no sqrt')
        else:
            if self.sqrt_gain is None: raise ValueError('ready band requires sqrt(gain) witness')
            if self.noise_sigma!=B.mul(self.bench_sigma,self.sqrt_gain):
                raise ValueError('band noise sigma detached from source-order bench*sqrt multiply')


def evaluate(band:L.State,*,bench_sigma,sqrt_gain=None,bits=96):
    if not isinstance(band,L.State): raise TypeError('persistent machine band State required')
    bench=_q(bench_sigma,'bench sigma')
    if bench<0: raise ValueError('bench sigma must be nonnegative')
    if not band.machine.ready:
        if sqrt_gain is not None: raise ValueError('unready machine band consumes no sqrt(gain) witness')
        return Result(band,bench,None,bench)
    if sqrt_gain is None: raise ValueError('ready machine band requires sqrt(gain) witness')
    sg=_q(sqrt_gain,'sqrt gain')
    gain=max(B.rn32(0),band.machine.p11)
    lo,hi=ROOT.sqrt_enclosure(gain,bits)
    if not RNE._interval_hits_rne_cell(lo,hi,sg):
        raise ValueError('band-noise sqrt witness detached from SAME stored p11 RNE cell')
    return Result(band,bench,sg,B.mul(bench,sg))


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=('if (!sigma_wave_band_.isReady()) {','return acc_noise_floor_sigma_;',
      'const float gain = sigma_wave_band_.whiteNoiseVarianceGain();',
      'return acc_noise_floor_sigma_ * std::sqrt(gain);')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_band_noise_floor_source_shape_matches':_source_shape_matches(),
      'unready_machine_band_returns_bench_sigma_by_identity':True,
      'ready_machine_band_consumes_same_stored_p11_gain':True,
      'sqrt_gain_related_to_exact_root_by_binary32_RNE_cell':True,
      'bench_times_sqrt_gain_binary32_multiply_materialized':True,
      'target_sqrt_libm_correspondence_closed':False,
      'machine_band_execution_membership_closed':False,
      'common_TuneState_boundary_can_consume_this_noise_floor':True,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
