"""Finite AdaptiveWaveBandPass + SeaStateAutoTuner variance recurrence.

Shipping timing is deliberately asymmetric inside one IMU sample.  The sigma
band/tuner runs BEFORE ``wave_period_.update``.  Therefore its external tuning
frequency is a read-only view of the WPE state carried into the sample: the
fixed prior until that pre-sample WPE state has latched usable, then the
canonical frequency represented by that state.  Only later does the current
vertical sample advance WPE for the next IMU sample.

The adaptive band's corner uses the previous SeaStateAutoTuner frequency when
that state is ready, otherwise the pre-update WPE/prior external frequency.
The band-noise floor follows shipping literally: before the adaptive band has
completed a valid step it is the raw bench noise sigma; only a ready band uses
bench_sigma*sqrt(p11).  Transcendental exp/sqrt binary32 ancestry remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE


def R(x): return P.rational(x)
def clamp(x,lo,hi): return min(max(x,lo),hi)

EMA_SCALE_MIN,EMA_SCALE_MAX=F(1,2),F(6)
EMA_HORIZON_MIN,EMA_HORIZON_MAX=F(1,20),F(35)


@dataclass(frozen=True)
class BandConfig:
    low_ratio:F; high_ratio:F; min_hz:F; max_hz:F
    tune_freq_floor:F; tune_freq_ceil:F; tune_freq_prior:F=F(1,5)
    def __post_init__(self):
        for n in ('low_ratio','high_ratio','min_hz','max_hz','tune_freq_floor','tune_freq_ceil','tune_freq_prior'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.low_ratio<=0 or self.high_ratio<=self.low_ratio or self.min_hz<=0 or self.max_hz<=self.min_hz:
            raise ValueError('invalid adaptive band configuration')
        if self.tune_freq_floor<=0 or self.tune_freq_ceil<self.tune_freq_floor or self.tune_freq_prior<=0:
            raise ValueError('invalid tune frequency bounds/prior')


@dataclass(frozen=True)
class BandState:
    lowpass_low:F=F(0); band:F=F(0)
    p00:F=F(0); p01:F=F(0); p11:F=F(0)
    low_hz:F|None=None; high_hz:F|None=None; ready:bool=False
    def __post_init__(self):
        for n in ('lowpass_low','band','p00','p01','p11'): object.__setattr__(self,n,R(getattr(self,n)))
        if self.low_hz is not None: object.__setattr__(self,'low_hz',R(self.low_hz))
        if self.high_hz is not None: object.__setattr__(self,'high_hz',R(self.high_hz))
        if not isinstance(self.ready,bool): raise TypeError('band ready flag must be literal')
        if self.p00<0 or self.p11<0: raise ValueError('band covariance diagonal must be nonnegative')


@dataclass(frozen=True)
class BandDecayWitness:
    q_low:F; q_high:F
    def __post_init__(self):
        for n in ('q_low','q_high'): object.__setattr__(self,n,R(getattr(self,n)))
        if not 0<self.q_low<=1 or not 0<self.q_high<=1: raise ValueError('band exp decays must lie in (0,1]')


@dataclass(frozen=True)
class StatsConfig:
    K_periods:F; tau_var_min:F; tau_var_max:F; f_min:F; f_max:F
    def __post_init__(self):
        for n in ('K_periods','tau_var_min','tau_var_max','f_min','f_max'): object.__setattr__(self,n,R(getattr(self,n)))
        if self.K_periods<=0 or self.tau_var_min<=0 or self.tau_var_max<self.tau_var_min or self.f_min<=0 or self.f_max<self.f_min:
            raise ValueError('invalid tuner-statistics configuration')


@dataclass(frozen=True)
class StatsState:
    frequency:F|None=None; tau_var:F=F(0)
    mean_value:F=F(0); mean_weight:F=F(0)
    sq_value:F=F(0); sq_weight:F=F(0)
    def __post_init__(self):
        if self.frequency is not None: object.__setattr__(self,'frequency',R(self.frequency))
        for n in ('tau_var','mean_value','mean_weight','sq_value','sq_weight'): object.__setattr__(self,n,R(getattr(self,n)))
        if self.mean_weight<0 or self.sq_weight<0: raise ValueError('negative EMA weight')
    @property
    def freq_ready(self): return self.frequency is not None and self.frequency>0
    @property
    def var_ready(self): return self.mean_weight>F(1,10**6) and self.sq_weight>F(1,10**6)


@dataclass(frozen=True)
class VarianceDecayWitness:
    decay:F
    def __post_init__(self):
        d=R(self.decay)
        if not 0<d<=1: raise ValueError('variance EMA decay must lie in (0,1]')
        object.__setattr__(self,'decay',d)


@dataclass(frozen=True)
class NoiseSqrtWitness:
    sqrt_gain:F
    def __post_init__(self):
        g=R(self.sqrt_gain)
        if g<0: raise ValueError('noise sqrt gain must be nonnegative')
        object.__setattr__(self,'sqrt_gain',g)


@dataclass(frozen=True)
class WPEFrequencyView:
    """Read-only tuner-visible output of the WPE state at sample entry."""
    state: WPE.WPEState
    period:F|None=None
    frequency:F|None=None
    def __post_init__(self):
        if not isinstance(self.state,WPE.WPEState): raise TypeError('carried WPE state required')
        p=None if self.period is None else R(self.period)
        f=None if self.frequency is None else R(self.frequency)
        if self.state.usable_period:
            if self.state.log_period is None or p is None or f is None or p<=0 or f<=0 or p*f!=1:
                raise ValueError('usable pre-update WPE state requires one canonical period/frequency pair')
        else:
            # Shipping ignores getFrequencyHz() until hasUsablePeriod(); do not
            # let an unnecessary frequency witness influence the tuner branch.
            if p is not None or f is not None:
                raise ValueError('unusable pre-update WPE tuner view consumes no canonical output witness')
        object.__setattr__(self,'period',p); object.__setattr__(self,'frequency',f)


def wpe_frequency_view(state:WPE.WPEState,*,period=None,frequency=None):
    return WPEFrequencyView(state,period,frequency)


@dataclass(frozen=True)
class FrontendResult:
    band_state:BandState
    stats_state:StatsState
    band_reference_frequency:F
    external_tuner_frequency:F
    current_tuner_frequency:F
    filtered_accel:F
    variance_ready:bool
    accel_variance:F
    band_noise_sigma:F


def band_step(s:BandState,cfg:BandConfig,*,x,dt,f_ref,decay:BandDecayWitness):
    x,dt,f_ref=map(R,(x,dt,f_ref))
    if dt<=0 or f_ref<=0: raise ValueError('positive band dt/frequency required')
    upper=min(cfg.max_hz,F(45,100)/dt)
    if upper<=cfg.min_hz: return s
    low=max(cfg.min_hz,cfg.low_ratio*f_ref); low=min(low,upper/F(105,100))
    high=min(upper,cfg.high_ratio*f_ref); high=max(high,low*F(105,100)); high=min(high,upper)
    if high<=low: return s
    ql,qh=decay.q_low,decay.q_high; al,ah=1-ql,1-qh
    low_new=ql*s.lowpass_low+al*x
    hp=ql*(x-s.lowpass_low)
    band_new=qh*s.band+ah*hp
    a00=ql; a10=-ah*ql; a11=qh; b0=al; b1=ah*ql
    p00=max(F(0),a00*a00*s.p00+b0*b0)
    p01=a00*(a10*s.p00+a11*s.p01)+b0*b1
    p11=max(F(0),a10*a10*s.p00+2*a10*a11*s.p01+a11*a11*s.p11+b1*b1)
    return BandState(low_new,band_new,p00,p01,p11,low,high,True)


def stats_step(s:StatsState,cfg:StatsConfig,*,dt,accel,external_frequency,decay:VarianceDecayWitness):
    dt,a,f=map(R,(dt,accel,external_frequency))
    if dt<=0 or f<=0: raise ValueError('positive tuner-stat dt/frequency required')
    fe=clamp(f,cfg.f_min,cfg.f_max)
    sea=clamp(F(1,2)/fe,EMA_SCALE_MIN,EMA_SCALE_MAX)
    Teff=2*sea
    requested=clamp(cfg.K_periods*Teff,cfg.tau_var_min,cfg.tau_var_max)
    lo=min(dt,EMA_HORIZON_MAX) if dt>EMA_HORIZON_MIN else EMA_HORIZON_MIN
    tau=clamp(requested,lo,EMA_HORIZON_MAX)
    d=decay.decay; alpha=1-d
    mv=d*s.mean_value+alpha*a; mw=d*s.mean_weight+alpha
    sv=d*s.sq_value+alpha*a*a; sw=d*s.sq_weight+alpha
    return StatsState(fe,tau,mv,mw,sv,sw)


def variance(s:StatsState):
    if not s.var_ready: return F(0)
    mu=s.mean_value/s.mean_weight
    return max(F(0),s.sq_value/s.sq_weight-mu*mu)


def tuner_frequency(view:WPEFrequencyView,cfg:BandConfig):
    if not isinstance(view,WPEFrequencyView): raise TypeError('pre-update WPE tuner view required')
    return view.frequency if view.state.usable_period else cfg.tune_freq_prior


def frontend_step_from_view(band:BandState,stats:StatsState,view:WPEFrequencyView,*,vertical_accel,dt,
                            band_cfg:BandConfig,stats_cfg:StatsConfig,band_decay:BandDecayWitness,
                            variance_decay:VarianceDecayWitness,bench_noise_sigma,
                            noise_sqrt:NoiseSqrtWitness|None=None):
    """Literal shipping tuner step using the WPE state carried into this sample."""
    external_f=tuner_frequency(view,band_cfg)
    fref=stats.frequency if stats.freq_ready else external_f
    fref=clamp(fref,band_cfg.tune_freq_floor,band_cfg.tune_freq_ceil)
    bnext=band_step(band,band_cfg,x=vertical_accel,dt=dt,f_ref=fref,decay=band_decay)
    snext=stats_step(stats,stats_cfg,dt=dt,accel=bnext.band,external_frequency=external_f,decay=variance_decay)
    bench=R(bench_noise_sigma)
    if bench<0: raise ValueError('bench noise sigma must be nonnegative')
    if not bnext.ready:
        if noise_sqrt is not None: raise ValueError('unready adaptive band consumes no sqrt(gain) witness')
        noise=bench
    else:
        if noise_sqrt is None: raise ValueError('ready adaptive band requires sqrt(gain) witness')
        if noise_sqrt.sqrt_gain*noise_sqrt.sqrt_gain != max(F(0),bnext.p11):
            raise ValueError('band-noise sqrt detached from same covariance recurrence')
        noise=bench*noise_sqrt.sqrt_gain
    return FrontendResult(bnext,snext,fref,external_f,snext.frequency,bnext.band,snext.var_ready,variance(snext),noise)


def frontend_step(band:BandState,stats:StatsState,wpe:WPE.UpdateResult,**kwargs):
    """Legacy conditional helper using a post-update WPE result.

    It remains for isolated algebra tests.  Shipping temporal composition must
    use ``frontend_step_from_view`` before the current WPE update.
    """
    if not isinstance(wpe,WPE.UpdateResult): raise TypeError('finite WPE result required')
    if wpe.state.usable_period:
        view=WPEFrequencyView(wpe.state,wpe.period,wpe.frequency)
    else:
        view=WPEFrequencyView(wpe.state)
    return frontend_step_from_view(band,stats,view,**kwargs)


def readiness():
    return {
      'adaptive_band_state_and_covariance_materialized':True,
      'adaptive_band_identity_branches_materialized':True,
      'unready_band_noise_floor_branch_materialized':True,
      'pre_update_WPE_usable_gate_and_tune_frequency_prior_materialized':True,
      'tuner_before_current_WPE_update_materialized':True,
      'previous_tuner_frequency_drives_band_corner':True,
      'external_prior_or_preupdate_WPE_frequency_drives_variance_horizon':True,
      'debiased_first_second_moments_materialized':True,
      'band_noise_variance_gain_materialized':True,
      'band_exp_and_noise_sqrt_binary32_ancestry_attached':False,
      'vertical_accel_frontend_same_history_attached':False,
      'stillness_state_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
