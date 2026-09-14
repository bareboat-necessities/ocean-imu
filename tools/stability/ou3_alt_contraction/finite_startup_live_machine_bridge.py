"""Carry coherent startup WPE + whole machine TuneState unchanged across goLive.

``goLive`` executes no WPE/tau/sigma/R_S adaptation recurrence.  Therefore the
entire deployment history produced by ``finite_guarded_tuner_machine_deployment``
crosses the control handoff by identity, including the common pending bit.

The exact-real runtime bridge still constructs the shadow active parameters.
This module does NOT identify those with either machine compiler history; the
machine TuneState -> active-parameter binary32 commit correspondence remains an
explicit downstream obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_deployment as START
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as LIVE
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as MACHINE
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPE


@dataclass(frozen=True)
class Result:
    live:LIVE.Result
    machine:MACHINE.State
    wpe:WPE.State
    startup:START.State
    def __post_init__(self):
        if not isinstance(self.live,LIVE.Result) or not isinstance(self.startup,START.State):
            raise TypeError('goLive result and coherent startup machine product required')
        if not isinstance(self.machine,MACHINE.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('whole machine TuneState and WPE ledger required')
        if self.live.frontend_before != self.startup.base.frontend:
            raise ValueError('goLive result detached from coherent startup frontend')
        if self.machine != self.startup.machine or self.wpe != self.startup.base.wpe:
            raise ValueError('goLive machine histories detached from retained startup product')
        if self.machine.tau != self.startup.base.tau:
            raise ValueError('goLive machine tau detached from startup WPE/tau ancestry')
        if self.live.frontend_before.tuner.pending != self.machine.pending:
            raise ValueError('goLive exact and machine common pending bits disagree')


def bridge(startup:START.State,entry,fresh,**kwargs):
    if not isinstance(startup,START.State):
        raise TypeError('coherent startup whole-machine product required')
    if startup.base.frontend.tuner.stage!='TunerReady':
        raise ValueError('coherent machine goLive bridge requires TunerReady startup frontend')
    out=LIVE.bridge(entry,fresh,startup.base.frontend,**kwargs)
    return Result(out,startup.machine,startup.base.wpe,startup)


def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'coherent_startup_WPE_tau_sigma_RS_product_consumed':True,
      'goLive_executes_no_WPE_or_machine_TuneState_adaptation_recurrence':True,
      'goLive_preserves_WPE_and_whole_machine_TuneState_by_identity':True,
      'goLive_preserves_common_pending_bit_into_first_Live_boundary':True,
      'frontend_memory_and_exact_shadow_goLive_bridge_retained': bool(
          l['Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive'] and
          l['goLive_unconditional_active_parameters_derived_from_same_carried_TuneState']),
      'construction_to_TunerReady_whole_machine_history_shape_available':
          s['startup_global_compiler_track_coherence_includes_WPE_tau_sigma_RS'],
      'machine_TuneState_to_goLive_active_parameters_binary32_correspondence_closed':False,
      'every_admitted_startup_history_reaches_this_goLive_product':False,
      'all_target_libm_correspondence_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
