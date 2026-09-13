"""Binary32 shipping sigma-target graph for the ALT deployment proof.

Shipping derives sigma from the wave-band variance before the common tau/sigma
EMA. Every ordinary operation is exact binary32 here. ``sqrt`` and the optional
stillness ``exp`` are explicit witnesses bound to their SAME rounded arguments.
Both transcendental results are related to rigorous exact-real intervals through
their exact RNE rounding cells; machine values are never falsely identified with
exact-real roots. Platform-libm correspondence remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
ZERO=B.rn32(0); ONE=B.rn32(1); VAR_FLOOR=B.rn32(F(1,10**6)); SIGMA_FLOOR=B.rn32(F(1,20))
QUALIFICATION='OU3_ALT_SIGMA_BINARY32_V4'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


def _positive_rne_cell(q:F):
    """Closed nearest-even rounding cell around a positive normal binary32."""
    q=F(q)
    if q<=0 or not B.is_binary32(q): raise ValueError('positive normal binary32 required for RNE cell')
    e=B._floor_log2_positive(q); quantum=B._pow2(e-23)
    lower_step=quantum/2 if q==B._pow2(e) else quantum
    prev=q-lower_step; nxt=q+quantum
    return (prev+q)/2,(q+nxt)/2


def _interval_hits_rne_cell(lo,hi,rounded):
    clo,chi=_positive_rne_cell(F(rounded))
    return max(F(lo),clo)<=min(F(hi),chi)


def _exp_minus_unit_enclosure(x,terms=14):
    """Alternating-series enclosure of exp(-x) on 0<=x<=1."""
    x=F(x)
    if x<0 or x>1: raise ValueError('unit exp argument outside [0,1]')
    if not isinstance(terms,int) or terms<2 or terms%2: raise ValueError('even Taylor degree >=2 required')
    total=F(1); term=F(1); partial={0:total}
    for n in range(1,terms+2):
        term *= -x/F(n); total += term; partial[n]=total
    return partial[terms+1],partial[terms]


def exp_minus_enclosure(x,terms=14):
    """Rigorous rational enclosure of exp(-x) for 0<=x<=60.

    Range-reduce by a power of two: choose m=2^k with y=x/m<=1, enclose
    exp(-y) by the alternating Taylor series, then use

        exp(-x) = exp(-y)^m.

    Since the base interval is nonnegative, powering both endpoints preserves
    order exactly.  Repeated squaring is used only as exact rational algebra;
    this is a proof enclosure, not a model of the target libm implementation.
    """
    x=F(x)
    if x<0 or x>60: raise ValueError('stillness exp argument outside certified [0,60] domain')
    if x<=1: return _exp_minus_unit_enclosure(x,terms)
    m=1
    while x>m: m*=2
    lo,hi=_exp_minus_unit_enclosure(x/F(m),terms)
    power=m
    while power>1:
        lo*=lo; hi*=hi; power//=2
    return lo,hi


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
        e=_q(still_exp_result,'stillness exp result')
        elo,ehi=exp_minus_enclosure(st)
        if not _interval_hits_rne_cell(elo,ehi,e):
            raise ValueError('stillness exp witness detached from SAME rounded argument RNE cell')
        atten=min(max(e,ZERO),ONE)
        attenuated=B.mul(pre,atten)
    else:
        if still_exp_result is not None: raise ValueError('nonstill branch consumes no attenuation exp witness')
        atten=ONE; attenuated=pre
    vw=max(attenuated,VAR_FLOOR)
    if sqrt_result is None: raise TypeError('sigma target requires sqrt(var_wave) witness')
    sr=_q(sqrt_result,'sigma sqrt result')
    slo,shi=ROOT.sqrt_enclosure(vw,bits)
    if not _interval_hits_rne_cell(slo,shi,sr):
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
      'stillness_exp_tight_rational_enclosure_bound_to_same_argument':True,
      'stillness_exp_range_reduction_covers_full_0_to_60_second_machine_domain':True,
      'stillness_exp_binary32_result_related_by_exact_RNE_cell':True,
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
