"""Carry the coherent whole machine TuneState through the shipping goLive edge.

Shipping ``goLive`` unconditionally applies the current TuneState to the MEKF,
but it executes no tuner EMA update and does NOT consume
``online_tune_apply_pending_``.  Therefore the machine tau/sigma/R_S histories
and their one pending bit cross the control edge by identity, while each global
compiler history also produces its own whole-TuneState goLive commit.

This deliberately differs from the ordinary next-IMU pending boundary: goLive
applies even when the pending bit is false and preserves that bit afterward.
The exact-real frontend goLive bridge remains the shadow/control relation.  The
machine commits are kept distinct from its active parameters; bounding and
injecting that exact-vs-machine active-parameter displacement into the Live word
is a subsequent deployment obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_tunestate_deployment as START
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as LIVE
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPE

QUALIFICATION='OU3_ALT_STARTUP_LIVE_MACHINE_TUNESTATE_BRIDGE_V1'


@dataclass(frozen=True)
class Result:
    live:LIVE.Result
    machine:PRODUCT.State
    wpe:WPE.State
    separate_commit:C.CommitResult
    fma_commit:C.CommitResult
    separate_active:ACTIVE.ActiveParameters
    fma_active:ACTIVE.ActiveParameters
    startup:START.State
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.live,LIVE.Result) or not isinstance(self.startup,START.State):
            raise TypeError('goLive result and coherent whole-machine startup product required')
        if not isinstance(self.machine,PRODUCT.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('whole machine TuneState and WPE deployment states required')
        if not isinstance(self.separate_commit,C.CommitResult) or not isinstance(self.fma_commit,C.CommitResult):
            raise TypeError('both global compiler goLive commits required')
        if not isinstance(self.separate_active,ACTIVE.ActiveParameters) or not isinstance(self.fma_active,ACTIVE.ActiveParameters):
            raise TypeError('both global compiler active parameter states required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong whole-machine goLive qualification')
        if self.live.frontend_before!=self.startup.lower.frontend:
            raise ValueError('goLive result detached from whole-machine startup frontend')
        if self.machine!=self.startup.machine or self.wpe!=self.startup.lower.wpe:
            raise ValueError('goLive machine histories detached from retained startup histories')
        if self.live.frontend_live.tuner.pending!=self.machine.pending:
            raise ValueError('goLive failed to preserve common exact/machine pending bit')
        if self.separate_active!=ACTIVE.ActiveParameters.from_commit(self.separate_commit):
            raise ValueError('separate active parameters detached from goLive machine commit')
        if self.fma_active!=ACTIVE.ActiveParameters.from_commit(self.fma_commit):
            raise ValueError('FMA active parameters detached from goLive machine commit')


def _tune(machine:PRODUCT.State,mode):
    if mode=='separate':
        return C.TuneState(machine.tau.separate,machine.sigma.separate,machine.rs.separate)
    if mode=='fma':
        return C.TuneState(machine.tau.fma,machine.sigma.fma,machine.rs.fma)
    raise ValueError('unknown global compiler mode')


def bridge(startup:START.State,entry,fresh,*,
           separate_machine_rs_sqrt_scale=None,fma_machine_rs_sqrt_scale=None,
           **kwargs):
    if not isinstance(startup,START.State): raise TypeError('whole machine startup product required')
    if startup.lower.frontend.tuner.stage!='TunerReady':
        raise ValueError('whole-machine goLive bridge requires TunerReady startup frontend')
    exact=LIVE.bridge(entry,fresh,startup.lower.frontend,**kwargs)
    cfg=kwargs.get('commit_cfg')
    if not isinstance(cfg,C.CommitConfig): raise TypeError('goLive CommitConfig required')
    nf=exact.band_noise_floor_sigma
    rs_scale=kwargs.get('rs_scale',1)
    # goLive applies current TuneState unconditionally, independent of the
    # online pending bit, and synchronizes posterior P_aw to stationary Sigma.
    sc=C.commit(_tune(startup.machine,'separate'),cfg,pending=True,live=True,
                band_noise_floor_sigma=nf,rs_sqrt_scale=separate_machine_rs_sqrt_scale,
                sync_covariance=True,rs_scale=rs_scale)
    fc=C.commit(_tune(startup.machine,'fma'),cfg,pending=True,live=True,
                band_noise_floor_sigma=nf,rs_sqrt_scale=fma_machine_rs_sqrt_scale,
                sync_covariance=True,rs_scale=rs_scale)
    if not isinstance(sc,C.CommitResult) or not isinstance(fc,C.CommitResult):
        raise AssertionError('unconditional machine goLive commit missing')
    return Result(exact,startup.machine,startup.lower.wpe,sc,fc,
                  ACTIVE.ActiveParameters.from_commit(sc),ACTIVE.ActiveParameters.from_commit(fc),startup)


def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'whole_machine_startup_product_consumed_at_goLive':s['startup_frontend_machine_TuneState_product_attached'],
      'goLive_executes_no_machine_tuner_EMA_update':True,
      'goLive_preserves_tau_sigma_RS_and_common_pending_bit_by_identity':True,
      'goLive_unconditionally_commits_each_global_compiler_whole_TuneState':True,
      'machine_goLive_commits_synchronize_aw_covariance_like_shipping':True,
      'exact_frontend_memory_and_control_goLive_bridge_retained':l['Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive'],
      'goLive_carries_whole_machine_TuneState_product':True,
      'machine_active_parameters_identified_with_exact_shadow_active_parameters':False,
      'machine_active_parameter_displacement_bound_closed':False,
      'compiler_FP_contraction_mode_qualified':False,
      'startup_frontend_binary32_correspondence_closed':False,
      'all_target_libm_correspondence_closed':False,
      'every_admitted_startup_history_reaches_this_goLive_product':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
