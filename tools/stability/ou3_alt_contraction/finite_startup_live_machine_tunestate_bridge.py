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
from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_tunestate_deployment as START
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as LIVE
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as BOUND
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPE

from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as MF

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
    arithmetic:BOUND.Boundary
    noise_floors:tuple
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
        if not isinstance(self.arithmetic,BOUND.Boundary): raise TypeError('binary32 goLive transaction required')
        if self.arithmetic.arithmetic.before!=replace(self.machine,pending=True):
            raise ValueError('goLive transaction detached from same machine snapshot')
        if not self.arithmetic.arithmetic.live or not self.arithmetic.sync_covariance:
            raise ValueError('goLive must apply Live R_S and synchronize a_w')
        if (self.separate_commit,self.fma_commit)!=(self.arithmetic.separate_commit,self.arithmetic.fma_commit):
            raise ValueError('goLive commits detached from binary32 transaction')
        MF.require_boundary(self.startup.frontends,self.noise_floors,self.arithmetic)
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

    @property
    def frontends(self): return self.startup.frontends


def bridge(startup:START.State,entry,fresh,*,
           separate_noise_sqrt_gain=None,fma_noise_sqrt_gain=None,
           **kwargs):
    if not isinstance(startup,START.State): raise TypeError('whole machine startup product required')
    exact=LIVE.bridge(entry,fresh,startup.lower.frontend,**kwargs)
    cfg=kwargs.get('commit_cfg')
    if not isinstance(cfg,C.CommitConfig): raise TypeError('goLive CommitConfig required')
    if kwargs.get('rs_scale',1)!=1:
        raise ValueError('shipping goLive uses literal R_S scale one')
    # Force only the local commit selector: goLive applies even if the real
    # online pending bit is false.  Return the ORIGINAL machine state below;
    # its pending flag and all scalar histories cross goLive by identity.
    floors=MF.boundary_floors(startup.frontends,bench_noise_sigma=kwargs['bench_noise_sigma'],required=True,
        separate_sqrt_gain=separate_noise_sqrt_gain,fma_sqrt_gain=fma_noise_sqrt_gain)
    arithmetic=BOUND.imu_boundary(replace(startup.machine,pending=True),cfg,live=True,
        separate_band_noise_floor_sigma=floors[0].noise_sigma,
        fma_band_noise_floor_sigma=floors[1].noise_sigma,sync_covariance=True)
    sc=arithmetic.separate_commit; fc=arithmetic.fma_commit
    return Result(exact,startup.machine,startup.lower.wpe,sc,fc,
                  ACTIVE.ActiveParameters.from_commit(sc),ACTIVE.ActiveParameters.from_commit(fc),
                  startup,arithmetic,floors)



def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'whole_machine_startup_product_consumed_at_goLive':s['startup_frontend_machine_TuneState_product_attached'],
      'goLive_executes_no_machine_tuner_EMA_update':True,
      'goLive_preserves_tau_sigma_RS_and_common_pending_bit_by_identity':True,
      'goLive_unconditionally_commits_each_global_compiler_whole_TuneState':True,
      'machine_goLive_commits_synchronize_aw_covariance_like_shipping':True,
      'goLive_applied_parameters_come_from_binary32_commit_graph':True,
      'goLive_preserves_startup_machine_band_stats_histories':True,
      'goLive_noise_floors_derived_from_startup_machine_band':True,
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
