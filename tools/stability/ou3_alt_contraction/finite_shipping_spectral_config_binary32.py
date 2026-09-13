"""Shipping-default binary32 qualification for the ALT SpectralMSE tuner path.

Tau qualification alone is insufficient: a CandidateConfig may carry the
correct frequency/tau clamps while silently substituting another pseudo cadence,
sigma coefficient, SpectralMSE coefficient, R_S clamp or cached q_eff root.
This module pins those fields to the CURRENT SHIPPING DEFAULT source spellings
and requires the cached q_eff value to descend from the explicit default-r_a
cache-production witness.

The cache witness's numerical pow correctness remains open.  Therefore this is
configuration/provenance qualification, not full deployment correspondence.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as Q

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
PSEUDO_NOMINAL=B.rn32(F(15,1000)); TAU_NOMINAL=B.rn32(F(11,10))
PSEUDO_RATIO=B.div(PSEUDO_NOMINAL,TAU_NOMINAL)
PSEUDO_MIN=B.div(B.rn32(1),B.rn32(200)); PSEUDO_MAX=B.rn32(F(15,100))
SIGMA_COEFF=B.rn32(F(9,10)); RS_MSE_COEFF=B.rn32(F(538,10000))
MIN_RS=B.rn32(F(15,100)); MAX_RS=B.rn32(100); MAX_SIGMA=B.rn32(4)
QUALIFICATION='OU3_ALT_SHIPPING_SPECTRAL_CONFIG_BINARY32_V1'


@dataclass(frozen=True)
class Qualified:
    cfg:C.CandidateConfig
    qeff_cache:Q.CacheWitness
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cfg,C.CandidateConfig) or not isinstance(self.qeff_cache,Q.CacheWitness):
            raise TypeError('CandidateConfig and qeff cache witness required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong shipping SpectralMSE qualification')
        qualify(self.cfg,self.qeff_cache)


def qualify(cfg:C.CandidateConfig,cache:Q.CacheWitness):
    if not isinstance(cfg,C.CandidateConfig): raise TypeError('CandidateConfig required')
    if not isinstance(cache,Q.CacheWitness): raise TypeError('default-r_a qeff cache witness required')
    required={
      'pseudo_tau_ratio':PSEUDO_RATIO,
      'pseudo_min':PSEUDO_MIN,
      'pseudo_max':PSEUDO_MAX,
      'sigma_coeff':SIGMA_COEFF,
      'rs_mse_coeff':RS_MSE_COEFF,
      'min_RS':MIN_RS,
      'max_RS':MAX_RS,
      'max_sigma':MAX_SIGMA,
    }
    for name,value in required.items():
        if F(getattr(cfg,name))!=F(value):
            raise ValueError(f'SpectralMSE {name} detached from shipping binary32 default')
    if cfg.clamp_enabled is not True:
        raise ValueError('shipping SpectralMSE qualification requires clamp-enabled branch')
    if cache.r_a!=Q.default_r_a_binary32():
        raise ValueError('qeff cache r_a detached from shipping default source expression')
    Q.qualify_candidate_cache(cfg.qeff_pow,cache)
    return Qualified.__new__(Qualified) if False else cfg


def bind(cfg:C.CandidateConfig,cache:Q.CacheWitness):
    qualify(cfg,cache)
    obj=object.__new__(Qualified)
    object.__setattr__(obj,'cfg',cfg); object.__setattr__(obj,'qeff_cache',cache)
    object.__setattr__(obj,'qualification',QUALIFICATION)
    return obj


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f;',
      'constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f;',
      'PSEUDO_UPDATE_PERIOD_NOMINAL_S / PSEUDO_UPDATE_TAU_NOMINAL_S;',
      'constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT;',
      'constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.15f;',
      'constexpr float MIN_R_S     = 0.15f;',
      'constexpr float MAX_R_S     = 100.0f;',
      'constexpr float MAX_SIGMA_A = 4.0f;',
      'constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;',
      'float sigma_coeff_  = 0.9f;',
      'RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;')
    return all(x in s for x in needles)


def readiness():
    q=Q.readiness()
    return {
      'qualification':QUALIFICATION,
      'shipping_SpectralMSE_default_source_shape_matches':_source_shape_matches(),
      'pseudo_cadence_compiled_binary32_defaults_pinned':True,
      'sigma_coefficient_compiled_binary32_default_pinned':True,
      'SpectralMSE_CJ_compiled_binary32_default_pinned':True,
      'RS_and_sigma_clamps_compiled_binary32_defaults_pinned':True,
      'default_r_a_source_expression_pinned':q['default_r_a_source_order_binary32_graph_materialized'],
      'candidate_qeff_cache_requires_explicit_default_r_a_cache_witness':True,
      'SpectralMSE_default_law_source_shape_pinned':_source_shape_matches(),
      'qeff_cache_target_libm_correspondence_closed':False,
      'per_sample_sqrt_pow_target_libm_correspondence_closed':False,
      'shipping_SpectralMSE_machine_target_fully_qualified':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
