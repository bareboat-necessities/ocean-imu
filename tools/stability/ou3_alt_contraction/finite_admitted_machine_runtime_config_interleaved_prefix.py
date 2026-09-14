"""Bind machine guard/LPF configuration to the carried admitted runtime.

The lower admitted machine vertical/stillness word already executes the binary32
filter API, vibration guard, private Mahony observer, tracker LPF and stillness
recurrence on one admitted IMU history.  Its remaining configuration gap is
that the machine guard ``Config`` and tracker-LPF cutoff can still be supplied
independently when that lower word is constructed.

This layer removes those two independent choices for the current theorem scope:

* the machine vibration-guard configuration is exactly the binary32 projection
  of ``RuntimeConfig.guard_cfg`` already carried by the source-owning word;
* the tracker LPF cutoff is the shipping default-constructor/reset value,
  ``MAX_FREQ_HZ == 6.0f``.  The represented theorem word has no configuration
  setter event, so both machine compiler histories must carry that same stored
  cutoff for the whole finite word.

This does not prove arbitrary calls to ``setFreqInputCutoffHz`` or mutable
private-Mahony gain setters.  Those remain separate startup/configuration
ancestry obligations.  It also does not close target-libm/compiler/Eigen
correspondence or authorize storage search.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_RUNTIME_CONFIG_INTERLEAVER_V1'
FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
DEFAULT_TRACKER_CUTOFF=B.rn32(6)


def _runtime(base:LOWER.State):
    if not isinstance(base,LOWER.State):
        raise TypeError('machine vertical/stillness state required')
    return LOWER._runtime(base.base)


def _machine_guard_cfg(runtime):
    """Binary32 object actually consumed by the machine guard graph."""
    c=runtime.guard_cfg
    return GUARD.Config(
        cutoff_hz=B.rn32(c.cutoff_hz),
        poles=c.poles,
        detect_hz=B.rn32(c.detect_hz),
        engage_lo=B.rn32(c.engage_lo),
        engage_hi=B.rn32(c.engage_hi),
        slew_tau=B.rn32(c.slew_tau),
        removed_rms_hz=B.rn32(c.removed_rms_hz),
        weight_epsilon=B.rn32(F(1,10000)),
    )


def _default_tracker_cutoff_source_matches():
    text=FILTER.read_text()
    return all(x in text for x in (
        'constexpr float MAX_FREQ_HZ = 6.0f;',
        'freq_input_lpf_.setCutoff(max_freq_hz_);',
        'float max_freq_hz_            = MAX_FREQ_HZ;',
        'void setFreqInputCutoffHz(float fc) {',
        'freq_input_lpf_.setCutoff(fc);',
    ))


def _validate_configuration(base:LOWER.State):
    runtime=_runtime(base)
    expected=_machine_guard_cfg(runtime)
    if base.guard_cfg!=expected:
        raise ValueError('machine guard configuration detached from carried RuntimeConfig.guard_cfg')
    if not _default_tracker_cutoff_source_matches():
        raise RuntimeError('shipping tracker-LPF default/reset source shape changed')
    for mode,source in (('separate',base.separate_source),('fma',base.fma_source)):
        if source.lpf.cutoff_hz!=DEFAULT_TRACKER_CUTOFF:
            raise ValueError(mode+' tracker LPF cutoff detached from shipping default-constructor/reset value')
    return runtime


@dataclass(frozen=True)
class State:
    base:LOWER.State
    entry_source_steps:int
    config_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State):
            raise TypeError('machine vertical/stillness state required')
        if self.qualification!=QUALIFICATION:
            raise ValueError('wrong machine runtime-config qualification')
        if not isinstance(self.entry_source_steps,int) or not isinstance(self.config_steps,int):
            raise TypeError('integer runtime-config counters required')
        if self.entry_source_steps<0 or self.config_steps<0:
            raise ValueError('nonnegative runtime-config counters required')
        if self.base.source_steps-self.entry_source_steps!=self.config_steps:
            raise ValueError('runtime-config validation count detached from machine source IMU count')
        _validate_configuration(self.base)


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('machine runtime-config IMU result malformed')
        if self.lower.state!=self.state.base:
            raise ValueError('runtime-config result detached from lower machine source result')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('runtime-config complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('runtime-config complete word detached from lower word')
        if self.state.config_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine runtime configuration not checked on all 600 admitted IMU edges')


def begin(base:LOWER.State):
    _validate_configuration(base)
    return State(base,base.source_steps,0)


def imu_step(state:State,**kwargs):
    if not isinstance(state,State):
        raise TypeError('machine runtime-config State required')
    lower=LOWER.imu_step(state.base,**kwargs)
    nxt=State(lower.state,state.entry_source_steps,state.config_steps+1)
    return ImuResult(nxt,lower)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State):
        raise TypeError('machine runtime-config State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.entry_source_steps,state.config_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State):
        raise TypeError('machine runtime-config State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.entry_source_steps,state.config_steps),event


def complete(state:State):
    return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness()
    return {
      'machine_vertical_stillness_word_consumed':low['complete_word_requires_machine_source_on_all_600_IMU_edges'],
      'machine_guard_config_binary32_projection_bound_to_carried_RuntimeConfig':True,
      'machine_guard_runtime_config_ancestry_closed_for_current_word':True,
      'tracker_LPF_shipping_default_constructor_reset_source_shape_checked':_default_tracker_cutoff_source_matches(),
      'tracker_LPF_cutoff_bound_to_default_constructor_reset_for_current_word':True,
      'separate_and_FMA_tracker_LPF_cutoff_cannot_diverge':True,
      'complete_word_requires_runtime_config_check_on_all_600_IMU_edges':True,
      'tracker_LPF_mutable_setter_ancestry_closed':False,
      'private_Mahony_mutable_config_setter_ancestry_closed':False,
      'startup_machine_runtime_config_history_attached':False,
      'target_libm_Eigen_and_compiler_profile_correspondence_closed':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
