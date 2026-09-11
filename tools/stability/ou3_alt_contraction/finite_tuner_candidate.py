"""Finite default SeaStateAutoTuner candidate/smoothing relation for ALT.

The shipping-level entry consumes the finite adaptive-band/statistics successor:
its bounded tuner frequency, acceleration variance readiness/value and propagated
band-noise sigma all come from one same-history frontend state. It then applies
the deployed SpectralMSE/default-slew-zero target and EMA logic before the staged
commit. Stillness and transcendental binary32 ancestry remain open obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState, clamp

EMA_SCALE_MIN, EMA_SCALE_MAX = F(1,2), F(6)
EMA_HORIZON_MIN, EMA_HORIZON_MAX = F(1,20), F(35)
VAR_WAVE_MIN = F(1,10**6)
SIGMA_AB_MIN = F(1,10**6)


@dataclass(frozen=True)
class CandidateConfig:
    min_freq: F; max_freq: F
    tau_coeff: F; sigma_coeff: F
    min_tau: F; max_tau: F; max_sigma: F
    pseudo_tau_ratio: F; pseudo_min: F; pseudo_max: F
    min_RS: F; max_RS: F
    rs_mse_coeff: F; accel_noise_density: F; qeff_pow: F
    adapt_tau_sec: F; adapt_tau_sea_periods: F
    adapt_RS_mult: F; adapt_RS_slew_log: F
    adapt_every_sec: F
    clamp_enabled: bool = True
    def __post_init__(self):
        names=('min_freq','max_freq','tau_coeff','sigma_coeff','min_tau','max_tau','max_sigma','pseudo_tau_ratio','pseudo_min','pseudo_max','min_RS','max_RS','rs_mse_coeff','accel_noise_density','qeff_pow','adapt_tau_sec','adapt_tau_sea_periods','adapt_RS_mult','adapt_RS_slew_log','adapt_every_sec')
        vals=[P.rational(getattr(self,n)) for n in names]
        for n,v in zip(names,vals): object.__setattr__(self,n,v)
        if not isinstance(self.clamp_enabled,bool): raise TypeError('literal clamp branch required')
        if self.min_freq<=0 or self.max_freq<self.min_freq or self.tau_coeff<=0 or self.sigma_coeff<=0: raise ValueError('invalid frequency/tuner coefficients')
        if self.min_tau<=0 or self.max_tau<self.min_tau or self.max_sigma<0: raise ValueError('invalid tau/sigma clamps')
        if self.pseudo_tau_ratio<=0 or self.pseudo_min<=0 or self.pseudo_max<self.pseudo_min: raise ValueError('invalid pseudo cadence')
        if self.max_RS<self.min_RS or self.rs_mse_coeff<=0 or self.accel_noise_density<=0 or self.qeff_pow<=0: raise ValueError('invalid SpectralMSE configuration')
        if self.qeff_pow**14 != 2*self.accel_noise_density: raise ValueError('qeff_pow detached from same accel-noise density')
        if self.adapt_tau_sec<=0 or self.adapt_tau_sea_periods<0 or self.adapt_RS_mult<=0 or self.adapt_RS_slew_log!=0 or self.adapt_every_sec<0:
            raise ValueError('this finite lemma covers deployed slew_log=0 and valid horizons')


@dataclass(frozen=True)
class WaveBandSample:
    frequency_hz: F
    variance_ready: bool
    accel_variance: F
    band_noise_sigma: F
    still: bool
    still_time: F
    still_attenuation: F
    sigma_wave_sqrt: F
    def __post_init__(self):
        for n in ('frequency_hz','accel_variance','band_noise_sigma','still_time','still_attenuation','sigma_wave_sqrt'):
            object.__setattr__(self,n,P.rational(getattr(self,n)))
        if not isinstance(self.variance_ready,bool) or not isinstance(self.still,bool): raise TypeError('literal variance/stillness branches required')
        if self.accel_variance<0 or self.band_noise_sigma<0 or self.still_time<0 or not 0<=self.still_attenuation<=1 or self.sigma_wave_sqrt<0: raise ValueError('invalid wave-band sample witness')
        if not self.still and self.still_attenuation != 1: raise ValueError('non-still branch must not attenuate variance')


@dataclass(frozen=True)
class TargetState:
    frequency: F
    variance_wave: F
    tau_target: F
    sigma_target: F


@dataclass(frozen=True)
class SpectralWitness:
    u_pow_6_7: F
    sqrt_TS: F
    def __post_init__(self):
        for n in ('u_pow_6_7','sqrt_TS'): object.__setattr__(self,n,P.rational(getattr(self,n)))
        if self.u_pow_6_7<0 or self.sqrt_TS<=0: raise ValueError('positive spectral witnesses required')


@dataclass(frozen=True)
class EmaWitness:
    decay_tau_sigma: F
    decay_RS: F
    def __post_init__(self):
        for n in ('decay_tau_sigma','decay_RS'): object.__setattr__(self,n,P.rational(getattr(self,n)))
        if not 0<self.decay_tau_sigma<=1 or not 0<self.decay_RS<=1: raise ValueError('EMA decay witnesses must lie in (0,1]')


@dataclass(frozen=True)
class CandidateResult:
    frequency: F
    variance_wave: F
    tau_target: F
    sigma_target: F
    RS_target: F
    tune_next: TuneState
    pending_after: bool
    last_adapt_time_after: F
    adapt_horizon: F
    RS_horizon: F


def _horizon_lo(dt):
    dt=P.rational(dt); return min(dt,EMA_HORIZON_MAX) if dt>EMA_HORIZON_MIN else EMA_HORIZON_MIN

def _clamp_horizon(x,dt): return clamp(P.rational(x),_horizon_lo(dt),EMA_HORIZON_MAX)


def targets(sample:WaveBandSample,cfg:CandidateConfig):
    f=clamp(sample.frequency_hz,cfg.min_freq,cfg.max_freq)
    noise_var=sample.band_noise_sigma**2
    total=max(F(0),sample.accel_variance) if sample.variance_ready else noise_var
    wave=max(F(0),total-noise_var)
    if sample.still: wave*=sample.still_attenuation
    wave=max(wave,VAR_WAVE_MIN)
    if sample.sigma_wave_sqrt**2 != wave: raise ValueError('sigma-wave sqrt witness detached from SAME variance branch')
    tau_raw=cfg.tau_coeff*F(1,2)/f
    if cfg.clamp_enabled:
        tau_t=clamp(tau_raw,cfg.min_tau,cfg.max_tau)
        sigma_t=min(sample.sigma_wave_sqrt*cfg.sigma_coeff,cfg.max_sigma)
    else:
        tau_t=tau_raw; sigma_t=sample.sigma_wave_sqrt*cfg.sigma_coeff
    return TargetState(f,wave,tau_t,sigma_t)


def spectral_RS(cfg:CandidateConfig,target:TargetState,w:SpectralWitness):
    TS=clamp(cfg.pseudo_tau_ratio*target.tau_target,cfg.pseudo_min,cfg.pseudo_max)
    if w.sqrt_TS*w.sqrt_TS != TS: raise ValueError('sqrt(T_S) witness detached from SAME tau-derived cadence')
    sigma_aB=max(target.sigma_target/cfg.sigma_coeff,SIGMA_AB_MIN)
    u=sigma_aB*target.tau_target**4
    if w.u_pow_6_7**7 != u**6: raise ValueError('6/7 power witness detached from SAME tau/sigma target')
    raw=cfg.rs_mse_coeff*cfg.qeff_pow*w.u_pow_6_7/w.sqrt_TS
    return clamp(raw,cfg.min_RS,cfg.max_RS) if cfg.clamp_enabled else raw


def step(previous:TuneState,sample:WaveBandSample,cfg:CandidateConfig,*,dt,time,last_adapt_time,
         spectral:SpectralWitness,ema:EmaWitness):
    dt,time,last_adapt_time=map(P.rational,(dt,time,last_adapt_time))
    if dt<=0 or time<last_adapt_time: raise ValueError('positive dt and monotone tuner time required')
    target=targets(sample,cfg); rs_t=spectral_RS(cfg,target,spectral)
    sea_time=F(1,2)/target.frequency
    if cfg.adapt_tau_sea_periods>0:
        safe=clamp(sea_time,EMA_SCALE_MIN,EMA_SCALE_MAX)
        adapt_h=_clamp_horizon(cfg.adapt_tau_sea_periods*safe,dt)
    else: adapt_h=cfg.adapt_tau_sec
    a=1-ema.decay_tau_sigma
    tune_tau=previous.tau_applied+a*(target.tau_target-previous.tau_applied)
    tune_sigma=previous.sigma_applied+a*(target.sigma_target-previous.sigma_applied)
    safe_tau=clamp(target.tau_target,EMA_SCALE_MIN,EMA_SCALE_MAX)
    RS_h=_clamp_horizon(cfg.adapt_RS_mult*safe_tau,dt)
    aRS=1-ema.decay_RS
    tune_RS=previous.RS_applied+aRS*(rs_t-previous.RS_applied)
    nxt=TuneState(tune_tau,tune_sigma,tune_RS)
    fire=(time-last_adapt_time)>cfg.adapt_every_sec
    return CandidateResult(target.frequency,target.variance_wave,target.tau_target,target.sigma_target,rs_t,nxt,fire,time if fire else last_adapt_time,adapt_h,RS_h)


def sample_from_frontend(front:BAND.FrontendResult,*,still,still_time,still_attenuation,sigma_wave_sqrt):
    if not isinstance(front,BAND.FrontendResult): raise TypeError('finite band/statistics successor required')
    return WaveBandSample(front.current_tuner_frequency,front.variance_ready,front.accel_variance,
                          front.band_noise_sigma,still,still_time,still_attenuation,sigma_wave_sqrt)


def step_from_frontend(previous:TuneState,front:BAND.FrontendResult,cfg:CandidateConfig,*,
                       still,still_time,still_attenuation,sigma_wave_sqrt,
                       dt,time,last_adapt_time,spectral:SpectralWitness,ema:EmaWitness):
    sample=sample_from_frontend(front,still=still,still_time=still_time,
                                still_attenuation=still_attenuation,sigma_wave_sqrt=sigma_wave_sqrt)
    return step(previous,sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,spectral=spectral,ema=ema)


def step_from_wpe(previous:TuneState,wpe:WPE.UpdateResult,sample:WaveBandSample,cfg:CandidateConfig,*,
                  dt,time,last_adapt_time,spectral:SpectralWitness,ema:EmaWitness):
    """Conditional shortcut valid only when SeaStateAutoTuner's f clamp is inactive."""
    if not isinstance(wpe,WPE.UpdateResult): raise TypeError('finite WPE successor required')
    if wpe.frequency is None or wpe.period is None: raise ValueError('finite canonical WPE output required')
    if sample.frequency_hz != wpe.frequency: raise ValueError('conditional WPE shortcut has active/detached tuner frequency clamp')
    return step(previous,sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,spectral=spectral,ema=ema)


def readiness():
    return {
      'frequency_variance_to_tau_sigma_targets':True,
      'default_SpectralMSE_target_same_tau_sigma_cadence':True,
      'tau_sigma_and_RS_EMA_recurrence':True,
      'sample_vs_commit_cadence_separated':True,
      'band_variance_tuner_frequency_same_history_attached':True,
      'stillness_state_attached':False,
      'exp_sqrt_pow_binary32_enclosed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
