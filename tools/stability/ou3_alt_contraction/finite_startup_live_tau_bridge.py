"""Carry the startup binary32 tau ledger unchanged across shipping ``goLive``.

``finite_startup_live_runtime_bridge`` already proves that goLive preserves the
frontend memories and commits the current exact TuneState into active OU/S
parameters. It does not execute another SeaStateAutoTuner EMA. Therefore the
machine tau ledger produced by ``finite_guarded_tuner_tau_deployment`` is a
literal identity across this control edge.

The result now retains the complete startup tau product as provenance, not just
a naked ledger. Its constructor rechecks that goLive consumed exactly that
frontend and that the returned ledger is exactly the invariant-checked startup
ledger. This prevents a boundary caller from pairing an unrelated binary32 tau
state with an otherwise valid goLive result.
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
    startup:START.State
    def __post_init__(self):
        if not isinstance(self.live,LIVE.Result) or not isinstance(self.tau,LEDGER.State) or not isinstance(self.startup,START.State):
            raise TypeError('goLive result, tau deployment ledger and startup tau product required')
        if self.live.frontend_before != self.startup.frontend:
            raise ValueError('goLive result detached from retained startup frontend state')
        if self.tau != self.startup.tau:
            raise ValueError('goLive tau ledger detached from retained startup tau invariant')


def bridge(startup:START.State,entry,fresh,**kwargs):
    if not isinstance(startup,START.State):
        raise TypeError('startup frontend tau product required')
    if startup.frontend.tuner.stage!='TunerReady':
        raise ValueError('tau goLive bridge requires TunerReady startup frontend')
    out=LIVE.bridge(entry,fresh,startup.frontend,**kwargs)
    return Result(out,startup.tau,startup)


def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'startup_tau_frontend_product_consumed':True,
      'goLive_executes_no_extra_tau_EMA':True,
      'goLive_preserves_binary32_tau_ledger_by_identity':True,
      'goLive_result_retains_invariant_checked_startup_tau_product':True,
      'frontend_memory_and_exact_TuneState_goLive_bridge_retained': bool(
          l['Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive'] and
          l['goLive_unconditional_active_parameters_derived_from_same_carried_TuneState']),
      'construction_seed_to_TunerReady_tau_ledger_shape_available':s['shipping_tau_ledger_begins_at_same_1p1f_exact_frontend_seed'],
      'exact_frontend_machine_tau_envelope_carried_to_goLive':s['exact_frontend_vs_each_machine_tau_track_error_envelope_inductively_checked'],
      'admitted_startup_reachability_with_tau_ledger_closed':False,
      'upstream_WPE_to_StoredFrequency_binary32_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'Live_admitted_product_consumes_this_exact_goLive_tau_ledger':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
