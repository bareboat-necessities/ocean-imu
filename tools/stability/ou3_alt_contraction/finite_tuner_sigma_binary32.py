"""Binary32 shipping sigma-target graph for the ALT deployment proof.

Shipping derives sigma from the wave-band variance before the common tau/sigma
EMA:

    var_noise = band_noise_sigma * band_noise_sigma
    var_total = var_ready ? max(0, accel_variance) : var_noise
    var_wave  = max(0, var_total - var_noise)
    if still: var_wave *= clamp(exp(-still_time),0,1)
    var_wave  = max(var_wave, 1e-6f)
    sigma_wave = sqrt(var_wave)
    sigma_target = min(sigma_wave * sigma_coeff_, max_sigma_a_)
    if !var_ready:
        sigma_target = max(sigma_target, max(0.05f, band_noise_sigma))

Every ordinary operation is exact binary32 here. ``sqrt`` and the optional
stillness ``exp`` are explicit witnesses bound to their SAME rounded arguments.
For sqrt, a deployed binary32 result is related to the rigorous exact-real root
through its exact RNE rounding cell; it is not falsely required to equal or lie
inside the narrow exact-real root interval. Platform-libm correspondence remains
open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
ZERO=B.rn32(0); ONE=B.rn32(1); VAR_FLOOR=B.rn32(F(1,10**6)); SIGMA_FLOOR=B.rn32(F(1,20))
QUALIFICATION='OU3_ALT_SIGMA_BINARY32_V2'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


def _positive_rne_cell(q:F):
    """Closed nearest-even rounding cell around a positive normal binary32.

    Closed tie boundaries are deliberately conservative: this helper proves
    that an exact-real root interval intersects the machine value's legal RNE
    cell. It does not assert which endpoint tie a platform sqrt implementation
    selects, and therefore does not close target-libm correspondence.
    """
    q=F(q)
    if q<=0 or not B.is_binary32(q): raise ValueError('positive normal binary32 required for RNE cell')
    e=B._floor_log2_positive(q); quantum=B._pow2(e-23)
    # Immediately below an exact power of two the predecessor lies in the
    # previous binade and has half the current-binade spacing.
    lower_step=quantum/2 if q==B._pow2(e) else quantum
    prev=q-lower_step; nxt=q+quantum
    return (prev+q)/2,(q+nxt)/2


def _root_interval_hits_rne_cell(lo,hi,rounded):
    clo,chi=_positive_rne_cell(F(rounded))
    return max(F(lo),clo)<=min(F(hi),chi)


@dataclass(frozen=True)
class Target:
    cfg:D.DeploymentConfig
    var_ready:bool
    accel_variance:F
    band_noise_sigma:F
    still:bool
    still_time:F
    attenuation:F
    var_noise:F
    var_total:F
    var_wave_pre_attenuation:F
    var_wave_attenuated:F
    var_wave:F
    sqrt_result:F
    sigma_wave:F
    scaled_sigma:F
    sigma_target:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
        if not isinstance(self.var_ready,bool) or not isinstance(self.still,bool): raise TypeError('literal readiness/stillness branches required')
        names=('accel_variance','band_noise_sigma','still_time','attenuation','var_noise','var_total',
               'var_wave_pre_attenuation','var_wave_attenuated','var_wave','sqrt_result','sigma_wave','scaled_sigma','sigma_target')
        for n in names: object.__setattr__(self,n,F(getattr(self,n)))
        if not all(B.is_binary32(getattr(self,n)) for n in names): raise ValueError('sigma graph stores binary32 values only')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong sigma target qualification')
        if self.accel_variance<0 or self.band_noise_sigma<0 or self.still_time<0: raise ValueError('nonnegative sigma source operands required')
        if not 0<=self.attenuation<=1: raise ValueError('stillness attenuation outside [0,1]')
        if self.var_wave<VAR_FLOOR or self.sqrt_result<=0 or self.sigma_target<=0: raise ValueError('invalid positive sigma target')


def target(cfg:D.DeploymentConfig,*,var_ready,accel_variance,band_noise_sigma,
           still,still_time=ZERO,still_exp_result=None,sqrt_result=None,bits=96):
    if not isinstance(cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
    if not isinstance(var_ready,bool) or not isinstance(still,bool): raise TypeError('literal readiness/stillness branches required')
    av=_q(accel_variance,'accel variance'); bn=_q(band_noise_sigma,'band noise sigma'); st=_q(still_time,'still time')
    if av<0 or bn<0 or st<0: raise ValueError('nonnegative sigma source operands required')
    vn=B.mul(bn,bn)
    vt=max(ZERO,av) if var_ready else vn
    pre=max(ZERO,B.sub(vt,vn))
    if still:
        if still_exp_result is None: raise TypeError('still branch requires exp attenuation witness')
        arg=st  # STILL_VAR_DECAY_SEC is literal 1.0f, so -still_time/1 has magnitude still_time.
        e=_q(still_exp_result,'stillness exp result')
        elo,ehi,_,_=EXP.enclosure(arg)
        if not elo<=e<=ehi: raise ValueError('stillness exp witness detached from SAME rounded argument')
        atten=min(max(e,ZERO),ONE)
        attenuated=B.mul(pre,atten)
    else:
        if still_exp_result is not None: raise ValueError('nonstill branch consumes no attenuation exp witness')
        atten=ONE; attenuated=pre
    vw=max(attenuated,VAR_FLOOR)
    if sqrt_result is None: raise TypeError('sigma target requires sqrt(var_wave) witness')
    sr=_q(sqrt_result,'sigma sqrt result')
    slo,shi=ROOT.sqrt_enclosure(vw,bits)
    if not _root_interval_hits_rne_cell(slo,shi,sr):
        raise ValueError('sigma sqrt witness detached from SAME rounded var_wave RNE cell')
    sw=sr
    scaled=B.mul(sw,cfg.sigma_coeff)
    sig=min(scaled,cfg.max_sigma) if cfg.clamp_enabled else scaled
    if not var_ready:
        sig=max(sig,max(SIGMA_FLOOR,bn))
    return Target(cfg,var_ready,av,bn,still,st,atten,vn,vt,pre,attenuated,vw,sr,sw,scaled,sig)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=('const float var_noise = band_noise_sigma * band_noise_sigma;',
      'float var_wave = var_total - var_noise;',
      'if (var_wave < 0.0f) var_wave = 0.0f;',
      'var_wave *= atten;',
      'var_wave = std::max(var_wave, 1e-6f);',
      'float sigma_wave = std::sqrt(var_wave);',
      'sigma_target_ = std::min(sigma_wave * sigma_coeff_,      max_sigma_a_);',
      'sigma_target_ = std::max(sigma_target_, std::max(0.05f, band_noise_sigma));')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_sigma_target_source_shape_matches':_source_shape_matches(),
      'noise_variance_subtraction_and_zero_floor_binary32_materialized':True,
      'optional_stillness_attenuation_binary32_multiply_materialized':True,
      'variance_1e_minus6_floor_binary32_materialized':True,
      'sigma_sqrt_witness_bound_to_same_rounded_var_wave':True,
      'sigma_sqrt_binary32_result_related_by_exact_RNE_cell_not_false_real_equality':True,
      'sigma_gain_max_clamp_and_unready_floor_binary32_materialized':True,
      'stillness_exp_target_libm_correspondence_closed':False,
      'sigma_sqrt_target_libm_correspondence_closed':False,
      'upstream_accel_variance_binary32_production_closed':False,
      'upstream_band_noise_sigma_binary32_production_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
