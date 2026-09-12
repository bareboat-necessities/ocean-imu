"""Carry the startup binary32 tau ledger unchanged across shipping ``goLive``.

``finite_startup_live_runtime_bridge`` already proves that goLive preserves the
frontend memories and commits the current exact TuneState into active OU/S
parameters.  It does not execute another SeaStateAutoTuner EMA.  Therefore the
machine tau ledger produced by ``finite_guarded_tuner_tau_deployment`` must be a
literal identity across this control edge.

This bridge makes that identity structural and returns the pair consumed by the
Live admitted tau interleaver.  It does not prove that every admitted startup
history reaches the supplied TunerReady frontend, nor WPE/libm correspondence.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_tau_deployment as START
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as LIVE
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER


@dataclass(frozen=True)
class Result:
    live:LIVE.Result
    tau:LEDGER.State
    def __post_init__(self):
        if not isinstance(self.live,LIVE.Result) or not isinstance(self.tau,LEDGER.State):
            raise TypeError('goLive result and tau deployment ledger required')


def bridge(startup:START.State,entry,fresh,**kwargs):
    if not isinstance(startup,START.State):
        raise TypeError('startup frontend tau product required')
    if startup.frontend.tuner.stage!='TunerReady':
        raise ValueError('tau goLive bridge requires TunerReady startup frontend')
    out=LIVE.bridge(entry,fresh,startup.frontend,**kwargs)
    return Result(out,startup.tau)


def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'startup_tau_frontend_product_consumed':True,
      'goLive_executes_no_extra_tau_EMA':True,
      'goLive_preserves_binary32_tau_ledger_by_identity':True,
      'frontend_memory_and_exact_TuneState_goLive_bridge_retained': bool(
          l['Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive'] and
          l['goLive_unconditional_active_parameters_derived_from_same_carried_TuneState']),
      'construction_seed_to_TunerReady_tau_ledger_shape_available':s['shipping_tau_ledger_begins_at_construction_seed'],
      'admitted_startup_reachability_with_tau_ledger_closed':False,
      'upstream_WPE_to_StoredFrequency_binary32_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'Live_admitted_product_consumes_this_exact_goLive_tau_ledger':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
