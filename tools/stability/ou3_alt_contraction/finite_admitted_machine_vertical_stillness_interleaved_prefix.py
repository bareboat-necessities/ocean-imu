"""Admitted Live TuneState word with machine guard/Mahony/LPF/stillness ancestry.

The lower word already carries one admitted BRMM/BIAS/ISS history, exact source
and shipping event, coherent WPE/tau/TuneState histories, and separate/FMA
machine band/statistics/sigma targets. This layer independently materializes the
shipping float ``update`` boundary and AccelVibrationGuard recurrence before the
private Mahony observer; it never borrows the exact-model guarded successor.

One common binary32 guard and one common private-Mahony successor feed both
compiler histories. Separate/FMA histories begin only at the tracker LPF and
stillness arithmetic, where permitted contraction choices may diverge.

The finite-real theorem runtime stores physical configuration values. Their
private-Mahony gains/gravity/settling time are explicitly projected through the
binary32 configuration store before entering the strict binary32 Mahony graph.
Mutable configuration setter provenance remains open; the projection is not a
claim that arbitrary exact-real values are executed directly by shipping C++.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAHONY
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as STARTUP
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_VERTICAL_STILLNESS_INTERLEAVER_V3'


def _runtime(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('machine TuneState state required')
    return base.base.base.prefix.prefix.live.live_word.runtime

def _machine_vertical_cfg(runtime):
    """Actual float configuration store seen by the private Mahony graph."""
    c=runtime.vertical_cfg
    return VERT.Config(B.rn32(c.two_kp),B.rn32(c.two_ki),B.rn32(c.gravity),B.rn32(c.settle_sec))

def _machine_imu_ordinal(base:LOWER.State):
    return base.machine.tau.updates-base.live_entry_machine_updates

@dataclass(frozen=True)
class State:
    base:LOWER.State
    guard:GUARD.State
    guard_cfg:GUARD.Config
    separate_source:VS.State
    fma_source:VS.State
    entry_machine_ordinal:int
    source_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('machine TuneState state required')
        if not isinstance(self.guard,GUARD.State) or not isinstance(self.guard_cfg,GUARD.Config): raise TypeError('persistent machine guard state/config required')
        if not isinstance(self.separate_source,VS.State) or not isinstance(self.fma_source,VS.State): raise TypeError('both machine vertical/stillness source histories required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted vertical/stillness qualification')
        if not isinstance(self.entry_machine_ordinal,int) or not isinstance(self.source_steps,int) or self.entry_machine_ordinal<0 or self.source_steps<0: raise ValueError('nonnegative machine source counters required')
        if _machine_imu_ordinal(self.base)-self.entry_machine_ordinal!=self.source_steps: raise ValueError('machine vertical/stillness count detached from TuneState IMU ordinal')
        n=self.base.frontends.samples
        if self.guard.samples!=n or self.separate_source.samples!=n or self.fma_source.samples!=n: raise ValueError('machine guard/vertical/stillness sample count detached from machine frontend history')
        if self.separate_source.vertical!=self.fma_source.vertical: raise ValueError('separate/FMA histories cannot duplicate or diverge private-Mahony state')

@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    guard:GUARD.Result
    separate_source:VS.Result
    fma_source:VS.Result
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult) or not isinstance(self.guard,GUARD.Result): raise TypeError('admitted guarded vertical/stillness result malformed')
        if not isinstance(self.separate_source,VS.Result) or not isinstance(self.fma_source,VS.Result): raise TypeError('both machine vertical/stillness results required')
        if self.separate_source.mahony!=self.fma_source.mahony: raise ValueError('compiler histories diverged in common private-Mahony source')
        if self.state.guard!=self.guard.state: raise ValueError('persistent machine guard detached from same-event result')
        if self.state.separate_source!=self.separate_source.state or self.state.fma_source!=self.fma_source.state: raise ValueError('persistent machine source state detached from same-event result')

@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord): raise TypeError('vertical/stillness complete word requires lower complete word')
        if self.lower.state!=self.state.base: raise ValueError('vertical/stillness complete word detached from lower word')
        if self.state.source_steps!=SOURCE.TRANSITIONS: raise ValueError('machine guard/vertical/stillness ancestry not attached to all 600 IMU edges')

def begin(base:LOWER.State,*,guard:GUARD.State,guard_cfg:GUARD.Config,separate_source:VS.State,fma_source:VS.State):
    return State(base,guard,guard_cfg,separate_source,fma_source,_machine_imu_ordinal(base),0)

def _source_step(source,guarded:GUARD.Result,runtime,h,*,lpf_alpha_exp,lpf_successor,still_energy_successor,still_attenuation_exp):
    hq=B.rn32(h); vcfg=_machine_vertical_cfg(runtime)
    mah=STARTUP.step(source.vertical,vcfg,dt=hq,gyro=guarded.raw_gyro,acc=guarded.conditioned_acc)
    band_input=B.rn32(mah.vertical.vertical_accel)
    lp=VS.lpf_step(source.lpf,x=band_input,dt=hq,alpha_exp=lpf_alpha_exp,successor=lpf_successor)
    st=STILL.step(source.stillness,runtime.still_cfg,vertical_lp=lp.state.value,dt=hq,energy_successor=still_energy_successor,attenuation_exp=still_attenuation_exp)
    nxt=VS.State(mah.vertical.state,lp.state,st.state,source.samples+1)
    return VS.Result(source,nxt,mah,lp,st,band_input)

def imu_step(state:State,*,machine_dt,machine_gyro_body,machine_acc_body,
             guard_lp_alpha_exp=None,guard_lp_successors=None,guard_detect_gamma_exp=None,guard_detect_successors=None,
             guard_removed_beta_exp=None,guard_removed_ms_successor=None,guard_removed_rms_sqrt=None,guard_slew_exp=None,guard_weight_successor=None,guard_output_successor=None,
             separate_lpf_alpha_exp=None,separate_lpf_successor=None,separate_still_energy_successor,separate_still_attenuation_exp=None,
             fma_lpf_alpha_exp=None,fma_lpf_successor=None,fma_still_energy_successor,fma_still_attenuation_exp=None,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    restricted=kwargs.get('restricted')
    if restricted is None: raise TypeError('same admitted physical restriction required')
    lower=LOWER.imu_step(state.base,**kwargs); live=LOWER._live_result(lower.lower); runtime=_runtime(state.base); h=restricted.segment.h
    hq=GUARD.api_scalar(h,machine_dt,'filter update dt')
    gyro=GUARD.api_vec(live.guarded.raw_gyro_body,machine_gyro_body,'filter update gyro')
    acc=GUARD.api_vec(live.guarded.raw_accel_body,machine_acc_body,'filter update accel')
    guarded=GUARD.step(state.guard,state.guard_cfg,raw_gyro=gyro,raw_acc=acc,dt=hq,
        lp_alpha_exp=guard_lp_alpha_exp,lp_successors=guard_lp_successors,detect_gamma_exp=guard_detect_gamma_exp,detect_successors=guard_detect_successors,
        removed_beta_exp=guard_removed_beta_exp,removed_ms_successor=guard_removed_ms_successor,removed_rms_sqrt=guard_removed_rms_sqrt,
        slew_exp=guard_slew_exp,weight_successor=guard_weight_successor,output_successor=guard_output_successor)
    sep=_source_step(state.separate_source,guarded,runtime,h,lpf_alpha_exp=separate_lpf_alpha_exp,lpf_successor=separate_lpf_successor,still_energy_successor=separate_still_energy_successor,still_attenuation_exp=separate_still_attenuation_exp)
    fma=_source_step(state.fma_source,guarded,runtime,h,lpf_alpha_exp=fma_lpf_alpha_exp,lpf_successor=fma_lpf_successor,still_energy_successor=fma_still_energy_successor,still_attenuation_exp=fma_still_attenuation_exp)
    if sep.mahony!=fma.mahony: raise ValueError('separate/FMA machine histories lost common private-Mahony source')
    VS.require_frontend_input(sep,lower.separate_frontend); VS.require_frontend_input(fma,lower.fma_frontend)
    VS.require_sigma_stillness(sep,lower.separate_sigma_join.machine); VS.require_sigma_stillness(fma,lower.fma_sigma_join.machine)
    nxt=State(lower.state,guarded.state,state.guard_cfg,sep.state,fma.state,state.entry_machine_ordinal,state.source_steps+1)
    return ImuResult(nxt,lower,guarded,sep,fma)

def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.guard,state.guard_cfg,state.separate_source,state.fma_source,state.entry_machine_ordinal,state.source_steps),event

def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.guard,state.guard_cfg,state.separate_source,state.fma_source,state.entry_machine_ordinal,state.source_steps),event

def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))

def readiness():
    low=LOWER.readiness(); guard=GUARD.readiness(); src=VS.readiness(); mah=MAHONY.readiness()
    return {
      'admitted_machine_TuneState_word_consumed':low['Live_600_step_machine_TuneState_product_attached'],
      'one_common_binary32_guard_feeds_private_Mahony_and_band_history':True,
      'machine_guard_displacement_already_injected_into_Racc':False,
      'shipping_float_API_dt_gyro_accel_projection_is_mandatory':True,
      'same_machine_private_Mahony_vertical_feeds_band_and_tracker_LPF':True,
      'machine_band_input_bound_to_same_frontend_consumed_by_sigma_target':True,
      'machine_tracker_LPF_and_stillness_bound_to_same_sigma_target':True,
      'complete_word_requires_machine_source_on_all_600_IMU_edges':True,
      'guard_runtime_config_ancestry_closed':False,
      'private_Mahony_runtime_config_and_startup_ancestry_closed':False,
      'tracker_LPF_runtime_config_ancestry_closed':False,
      'guard_exp_sqrt_compiler_correspondence_closed':guard['target_libm_and_Eigen_expression_correspondence_closed'],
      'private_Mahony_compiler_profile_qualified':mah['actual_target_compiler_profile_qualified'],
      'tracker_LPF_exp_target_libm_correspondence_closed':src['target_libm_and_compiler_profile_correspondence_closed'],
      'stillness_exp_sqrt_target_libm_correspondence_closed':src['target_libm_and_compiler_profile_correspondence_closed'],
      'source_uniform_machine_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
