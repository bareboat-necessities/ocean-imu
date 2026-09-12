"""Tau-specific shipping qualification for a structural 600-transition word.

Low-level finite runtime modules intentionally accept alternative configs for
component algebra tests.  The theorem-facing complete shipping word may not.
This layer checks the persisted tau-relevant ``CandidateConfig`` against the
actual compiled binary32 defaults, then imports the analytic 600-step tau-domain
and roundoff theorem.

It still does NOT claim complete deployment arithmetic: the current interleaved
event ledger does not retain a binary32 TauStep/roundoff certificate for every
IMU edge.  The next master-word step is to attach those per-event certificates.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as WORD
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as LIVE
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_tau_600_domain as DOMAIN

ADAPT_SEC=B.rn32(F(9,5))
ADAPT_PERIODS=B.rn32(F(2,5))
QUALIFICATION='OU3_ALT_COMPLETE_WORD_TAU_CONFIG_V1'


@dataclass(frozen=True)
class QualifiedTauWord:
    word: WORD.CompleteWord
    def __post_init__(self):
        if not isinstance(self.word,WORD.CompleteWord):
            raise TypeError('structurally complete 600-transition word required')
        qualify_runtime(self.word.state.prefix.live.live_word.runtime)
        if not DOMAIN.assert_induction_closes():
            raise AssertionError('canonical 600-step tau domain theorem did not close')


def qualify_runtime(runtime:LIVE.RuntimeConfig):
    if not isinstance(runtime,LIVE.RuntimeConfig): raise TypeError('persistent RuntimeConfig required')
    c=runtime.candidate_cfg
    required={
      'min_freq':TARGET.FLOOR,
      'max_freq':TARGET.CEIL,
      'tau_coeff':B.rn32(1),
      'min_tau':TARGET.TAU_MIN,
      'max_tau':TARGET.TAU_MAX,
      'adapt_tau_sec':ADAPT_SEC,
      'adapt_tau_sea_periods':ADAPT_PERIODS,
    }
    for name,value in required.items():
        if F(getattr(c,name)) != F(value):
            raise ValueError(f'persistent tuner {name} detached from shipping binary32 default')
    if c.clamp_enabled is not True:
        raise ValueError('shipping tau theorem requires default clamp-enabled branch')
    return runtime


def readiness():
    target=TARGET.readiness(); domain=DOMAIN.readiness()
    return {
      'qualification':QUALIFICATION,
      'tau_relevant_persistent_runtime_config_source_qualified':True,
      'compiled_float_operands_not_ideal_decimal_substitutes':True,
      'shipping_tau_target_binary32_cell_bound_available':target['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'],
      'shipping_tau_predecessor_domain_inductively_closed_for_600_step_word':domain['shipping_tau_predecessor_domain_inductively_closed_for_600_step_word'],
      'source_uniform_tau_roundoff_supply_bound_available_for_600_step_word':domain['source_uniform_tau_roundoff_supply_bound_closed_for_600_step_word'],
      'every_IMU_event_retains_binary32_tau_step_and_supply_certificate':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
