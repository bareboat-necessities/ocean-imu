"""Deployment-faithful SpectralMSE tuner configuration for ALT.

The legacy ``finite_tuner_candidate.CandidateConfig`` is intentionally useful
for exact-rational component algebra, but it enforces

    qeff_pow**14 == 2 * accel_noise_density

as an exact rational identity.  A real shipping cache is a binary32 result of
``pow(2*r_a,1/14)`` and generally cannot satisfy that equality.  Consequently
that legacy type must not be used as the theorem representation of the shipping
SpectralMSE configuration.

This module carries the same ordinary tuner constants but represents q_eff with
its actual provenance: the configured binary32 ``r_a`` and an explicit cached
binary32 pow witness.  The exact mathematical 14th root remains available from
the witness's rigorous real enclosure; the cached machine value is a distinct
quantity exposed through ``qeff_pow`` for deployment graphs.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as Q

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_TUNER_DEPLOYMENT_CONFIG_V1'

# CURRENT SHIPPING DEFAULT compiled-float values.
MIN_FREQ=B.rn32(F(3,100)); MAX_FREQ=B.rn32(F(6,5)); TAU_COEFF=B.rn32(1)
SIGMA_COEFF=B.rn32(F(9,10)); MIN_TAU=B.rn32(F(1,50)); MAX_TAU=B.rn32(12)
MAX_SIGMA=B.rn32(4)
PSEUDO_NOMINAL=B.rn32(F(15,1000)); TAU_NOMINAL=B.rn32(F(11,10))
PSEUDO_RATIO=B.div(PSEUDO_NOMINAL,TAU_NOMINAL)
PSEUDO_MIN=B.div(B.rn32(1),B.rn32(200)); PSEUDO_MAX=B.rn32(F(15,100))
MIN_RS=B.rn32(F(15,100)); MAX_RS=B.rn32(100); RS_MSE_COEFF=B.rn32(F(538,10000))
ADAPT_TAU_SEC=B.rn32(F(18,10)); ADAPT_TAU_SEA_PERIODS=B.rn32(F(4,10))
ADAPT_RS_MULT=B.rn32(F(15,10)); ADAPT_RS_SLEW_LOG=B.rn32(0); ADAPT_EVERY_SEC=B.rn32(F(1,10))


@dataclass(frozen=True)
class DeploymentConfig:
    min_freq:F; max_freq:F
    tau_coeff:F; sigma_coeff:F
    min_tau:F; max_tau:F; max_sigma:F
    pseudo_tau_ratio:F; pseudo_min:F; pseudo_max:F
    min_RS:F; max_RS:F
    rs_mse_coeff:F
    accel_noise_density:F
    qeff_cache:Q.CacheWitness
    adapt_tau_sec:F; adapt_tau_sea_periods:F
    adapt_RS_mult:F; adapt_RS_slew_log:F
    adapt_every_sec:F
    clamp_enabled:bool=True
    qualification:str=QUALIFICATION
    def __post_init__(self):
        names=('min_freq','max_freq','tau_coeff','sigma_coeff','min_tau','max_tau','max_sigma',
               'pseudo_tau_ratio','pseudo_min','pseudo_max','min_RS','max_RS','rs_mse_coeff',
               'accel_noise_density','adapt_tau_sec','adapt_tau_sea_periods','adapt_RS_mult',
               'adapt_RS_slew_log','adapt_every_sec')
        for n in names: object.__setattr__(self,n,P.rational(getattr(self,n)))
        if not isinstance(self.qeff_cache,Q.CacheWitness): raise TypeError('explicit qeff cache witness required')
        if not isinstance(self.clamp_enabled,bool): raise TypeError('literal clamp branch required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong deployment tuner qualification')
        if self.min_freq<=0 or self.max_freq<self.min_freq or self.tau_coeff<=0 or self.sigma_coeff<=0:
            raise ValueError('invalid deployment frequency/tuner coefficients')
        if self.min_tau<=0 or self.max_tau<self.min_tau or self.max_sigma<0:
            raise ValueError('invalid deployment tau/sigma clamps')
        if self.pseudo_tau_ratio<=0 or self.pseudo_min<=0 or self.pseudo_max<self.pseudo_min:
            raise ValueError('invalid deployment pseudo cadence')
        if self.max_RS<self.min_RS or self.rs_mse_coeff<=0 or self.accel_noise_density<=0:
            raise ValueError('invalid deployment SpectralMSE configuration')
        if self.adapt_tau_sec<=0 or self.adapt_tau_sea_periods<0 or self.adapt_RS_mult<=0 or self.adapt_every_sec<0:
            raise ValueError('invalid deployment adaptation horizons')
        if self.adapt_RS_slew_log!=0: raise ValueError('current theorem covers deployed slew_log=0 only')
        if self.accel_noise_density!=self.qeff_cache.r_a:
            raise ValueError('deployment r_a detached from cached qeff pow source')

    @property
    def qeff_pow(self):
        """Actual cached binary32 value, for the machine arithmetic graph only."""
        return self.qeff_cache.result

    @property
    def qeff_exact_root_interval(self):
        """Exact mathematical (2 r_a)^(1/14), distinct from cached float."""
        return self.qeff_cache.true_lo,self.qeff_cache.true_hi


def shipping_defaults(*,qeff_pow_result):
    cache=Q.produce(r_a=Q.default_r_a_binary32(),pow_result=qeff_pow_result)
    return DeploymentConfig(
      MIN_FREQ,MAX_FREQ,TAU_COEFF,SIGMA_COEFF,MIN_TAU,MAX_TAU,MAX_SIGMA,
      PSEUDO_RATIO,PSEUDO_MIN,PSEUDO_MAX,MIN_RS,MAX_RS,RS_MSE_COEFF,
      cache.r_a,cache,ADAPT_TAU_SEC,ADAPT_TAU_SEA_PERIODS,ADAPT_RS_MULT,
      ADAPT_RS_SLEW_LOG,ADAPT_EVERY_SEC,True)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'constexpr float MIN_TUNE_FREQ_HZ = 0.03f;',
      'constexpr float MAX_TUNE_FREQ_HZ = 1.2f;',
      'constexpr float MIN_TAU_S   = 0.02f;',
      'constexpr float MAX_TAU_S   = 12.0f;',
      'constexpr float MAX_SIGMA_A = 4.0f;',
      'constexpr float MIN_R_S     = 0.15f;',
      'constexpr float MAX_R_S     = 100.0f;',
      'constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;',
      'float tau_coeff_    = 1.0f;',
      'float sigma_coeff_  = 0.9f;',
      'RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;',
      'rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_default_source_shape_matches':_source_shape_matches(),
      'legacy_exact_rational_qeff_identity_not_required':True,
      'configured_r_a_and_cached_qeff_machine_value_carried_distinctly':True,
      'exact_qeff_fourteenth_root_interval_carried_separately_from_cache':True,
      'all_current_shipping_SpectralMSE_scalar_defaults_materialized_as_binary32':True,
      'qeff_cache_target_libm_correspondence_closed':False,
      'per_sample_sqrt_pow_target_libm_correspondence_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
