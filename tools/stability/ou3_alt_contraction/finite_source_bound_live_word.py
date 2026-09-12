"""Source-owning finite Live product for the ALT proof.

The lower ``finite_live_interleave`` layer composes exact shipping event
relations but accepts source-qualified transitions and per-call runtime operands.
For the universal word that is still too weak in two ways: independently valid
source transitions could be spliced, and persistent shipping configuration could
be silently changed between events.

This module carries a checked necessary BRMM/BIAS outer continuation and the
static runtime configuration in the theorem product state. Every IMU event
appends exactly the next source transition and injects the same carried Qbase,
nominal Racc, guard/WPE/band/stillness/tuner configs and bench-noise constants.
Those values are no longer theorem-step inputs. Dynamic transcendental/solver
witnesses remain per-event, but their component modules check them against the
same carried state expressions.

Magnetic and external-hold events do not advance the physical source. The
initial Live endpoint is inherited from startup. The source checks are only
necessary conditions, not full generator/O^601_BRMM membership; that sample-zero
handoff in O^601_BRMM remains explicitly open. Quantitative sensor ISS bounds
and deployment arithmetic correspondence also remain open, so storage and all
ALT theorem gates remain false.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_live_interleave as LIVE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT

R=M.rational
STATIC_KEYS=frozenset({
    'commit_cfg','boundary_bench_noise_sigma','guard_cfg','vertical_cfg',
    'band_cfg','racc_cfg','nominal_racc_std','Qbase','wpe_cfg','stats_cfg',
    'bench_noise_sigma','still_cfg','candidate_cfg',
})


@dataclass(frozen=True)
class RuntimeConfig:
    """Persistent shipping/model configuration that may not vary per IMU event."""
    commit_cfg: COMMIT.CommitConfig
    boundary_bench_noise_sigma: F
    guard_cfg: GUARD.Config
    vertical_cfg: VERT.Config
    band_cfg: BAND.BandConfig
    racc_cfg: RACC.Config
    nominal_racc_std: tuple
    Qbase: tuple
    wpe_cfg: WPE.WPEConfig
    stats_cfg: BAND.StatsConfig
    bench_noise_sigma: F
    still_cfg: STILL.Config
    candidate_cfg: CAND.CandidateConfig
    def __post_init__(self):
        if not isinstance(self.commit_cfg,COMMIT.CommitConfig): raise TypeError('persistent CommitConfig required')
        if not isinstance(self.guard_cfg,GUARD.Config): raise TypeError('persistent guard Config required')
        if not isinstance(self.vertical_cfg,VERT.Config): raise TypeError('persistent private-Mahony Config required')
        if not isinstance(self.band_cfg,BAND.BandConfig): raise TypeError('persistent band Config required')
        if not isinstance(self.racc_cfg,RACC.Config): raise TypeError('persistent Racc Config required')
        if not isinstance(self.wpe_cfg,WPE.WPEConfig): raise TypeError('persistent WPE Config required')
        if not isinstance(self.stats_cfg,BAND.StatsConfig): raise TypeError('persistent statistics Config required')
        if not isinstance(self.still_cfg,STILL.Config): raise TypeError('persistent stillness Config required')
        if not isinstance(self.candidate_cfg,CAND.CandidateConfig): raise TypeError('persistent tuner candidate Config required')
        object.__setattr__(self,'boundary_bench_noise_sigma',R(self.boundary_bench_noise_sigma))
        object.__setattr__(self,'bench_noise_sigma',R(self.bench_noise_sigma))
        n=tuple(M.vec(self.nominal_racc_std,3));
        if min(n)<=0: raise ValueError('positive persistent nominal Racc std required')
        object.__setattr__(self,'nominal_racc_std',n)
        q=tuple(tuple(x for x in row) for row in M.mat(self.Qbase,6,6))
        if q != tuple(tuple(x for x in row) for row in M.transpose(q)):
            raise ValueError('persistent Qbase must be symmetric')
        object.__setattr__(self,'Qbase',q)
    @classmethod
    def from_step_kwargs(cls,kw):
        missing=STATIC_KEYS-set(kw)
        if missing: raise ValueError('cannot extract persistent runtime config; missing '+repr(sorted(missing)))
        return cls(**{k:kw[k] for k in STATIC_KEYS})
    def inject(self,dynamic):
        overlap=STATIC_KEYS & set(dynamic)
        if overlap: raise TypeError('theorem IMU call cannot override carried runtime config '+repr(sorted(overlap)))
        out=dict(dynamic)
        for k in STATIC_KEYS: out[k]=getattr(self,k)
        return out


def dynamic_kwargs(kw):
    """Remove persistent configuration from a component-test operand bundle."""
    return {k:v for k,v in kw.items() if k not in STATIC_KEYS}


@dataclass(frozen=True)
class State:
    live: LIVE.State
    source: SOURCE.Continuation
    sensor_root: SOURCE.SensorDisturbanceRoot
    bias_history_id: str
    runtime: RuntimeConfig
    def __post_init__(self):
        if not isinstance(self.live,LIVE.State) or not isinstance(self.source,SOURCE.Continuation):
            raise TypeError('interleaved Live state and persistent source continuation required')
        if not isinstance(self.sensor_root,SOURCE.SensorDisturbanceRoot):
            raise TypeError('persistent sensor disturbance root required')
        if not isinstance(self.runtime,RuntimeConfig): raise TypeError('persistent runtime configuration required')
        if not isinstance(self.bias_history_id,str) or not self.bias_history_id:
            raise ValueError('persistent concrete bias-history identity required')
        if self.sensor_root.source_root != self.source.root:
            raise ValueError('sensor histories detached from carried source continuation')
        core=self.live.live.live.mekf; root=self.source.root
        if core.reference.history_id != root.history_id:
            raise ValueError('carried source history detached from current MEKF physical reference')
        if core.reference.live_origin != root.live_origin:
            raise ValueError('carried source root restarted one-time Live origin')
        if core.reference.bias_family != root.bias_family:
            raise ValueError('current physical bias family detached from carried analytic BIAS contract')
        if core.reference.bias_root != self.bias_history_id:
            raise ValueError('concrete physical bias history restarted inside carried source word')
        if self.source.steps:
            if self.source.steps[-1].segment.after != core.reference:
                raise ValueError('Live product state detached from last carried source endpoint')
        elif core.reference.time != root.live_origin:
            raise ValueError('empty source continuation must sit at fresh Live origin')


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_live(live:LIVE.State, root:SOURCE.SourceRoot,
              sensor_root:SOURCE.SensorDisturbanceRoot,runtime:RuntimeConfig):
    if not isinstance(live,LIVE.State) or not isinstance(root,SOURCE.SourceRoot):
        raise TypeError('Live product and source root required')
    if sensor_root.source_root != root:
        raise ValueError('sensor histories detached from source root')
    ref=live.live.live.mekf.reference
    if ref.bias_family != root.bias_family:
        raise ValueError('fresh Live bias family detached from selected BIAS source family')
    return State(live,SOURCE.begin(root),sensor_root,ref.bias_root,runtime)


def imu_step(state:State, *, witness:SOURCE.StepWitness,
             segment:PHYS.PhysicalSegment, raw:SENSOR.RawImuSample,
             packet_id:str, **dynamic):
    """Append next source transition and execute IMU with carried static config."""
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    if witness.ordinal != state.source.next_ordinal:
        raise ValueError('IMU event must consume exactly the next source ordinal')
    if 'dt' in dynamic:
        raise TypeError('theorem IMU dt is owned by the qualified physical segment')
    nxt_source=SOURCE.append(state.source,witness=witness,segment=segment)
    qualified=nxt_source.steps[-1]
    packet=SOURCE.qualify_raw_imu(qualified,state.sensor_root,raw,packet_id)
    operands=state.runtime.inject(dynamic); operands['dt']=segment.h
    event=LIVE.imu_step_source_qualified(state.live,packet,**operands)
    return Result(State(event.state,nxt_source,state.sensor_root,state.bias_history_id,state.runtime),event)


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    if state.source.steps:
        endpoint=SOURCE.endpoint(state.source.steps[-1],'after')
    else:
        # A magnetic call may occur at the fresh Live origin before the first
        # IMU transition.  Check the ACTUAL fresh Reference directly; do not
        # manufacture a predecessor segment or advance source ordinal/time.
        endpoint=SOURCE.origin_endpoint(state.source.root,state.live.live.live.mekf.reference)
    event=LIVE.mag_step_source_qualified(state.live,endpoint,**kwargs)
    return Result(State(event.state,state.source,state.sensor_root,state.bias_history_id,state.runtime),event)


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    event=LIVE.set_hold(state.live,hold=hold)
    return Result(State(event.state,state.source,state.sensor_root,state.bias_history_id,state.runtime),event)


def finite_storage_status():
    return {
      'map_representation':'finite_physical_descriptor',
      'finite_error_identity_for_every_event':False,
      'physical_reference_forcing_retained':True,
      # Persistent/static model coefficients are now state-owned, but dynamic
      # exp/trig/sqrt/solver products and a few frontend outputs still need a
      # complete same-source/deployment closure before this can become true.
      'all_coefficient_product_graphs_retained':False,
      'all_configured_branches_bound_to_finite_graph':False,
      'zero_wind_heel_scope_enforced':True,
    }


def readiness():
    lower=LIVE.readiness(); src=SOURCE.readiness(); master=finite_storage_status()
    return {
      'source_continuation_is_part_of_theorem_product_state':True,
      'physical_caps_and_coupled_moments_checked_before_IMU_execution':src['physical_vector_caps_and_joint_moment_IQC_checked'],
      'correlated_prefix_energy_derived_from_same_physical_segments':src['prefix_moments_and_derived_energy_not_independent_supply'],
      'full_source_membership_follows_from_history_tokens':False,
      'every_IMU_event_appends_exactly_next_source_ordinal':True,
      'source_cell_and_physical_primitive_chains_advance_with_filter_state':True,
      'raw_IMU_sensor_histories_advance_with_same_source_transition':True,
      'analytic_BIAS_family_token_and_concrete_bias_history_both_persist':True,
      'persistent_static_runtime_configuration_carried_in_product_state':True,
      'theorem_IMU_event_cannot_override_static_runtime_configuration':True,
      'theorem_IMU_dt_owned_by_qualified_physical_segment':True,
      'magnetic_and_hold_events_preserve_current_source_endpoint':True,
      'qualified_post_first_IMU_magnetic_entry_available':lower['source_checked_async_magnetic_endpoint_entry_available'],
      'sample_zero_magnetic_entry_uses_checked_fresh_physical_origin':src['sample_zero_checked_outer_endpoint_available_without_fake_transition'],
      'sample_zero_magnetic_entry_advances_no_source_ordinal':True,
      'sample_zero_full_source_membership_proved':False,
      'one_bias_parameter_token_carried_over_word':src['one_bias_family_parameter_token_over_word_required'],
      'physical_reference_forcing_retained_in_finite_master_status':master['physical_reference_forcing_retained'],
      'sample_zero_startup_to_checked_outer_endpoint_bridge_closed':True,
      'sample_zero_startup_to_COMPLETE_BRMM_endpoint_bridge_closed':False,
      'quantitative_sensor_residual_ISS_envelope_attached':False,
      'finite_estimator_coefficients_bound_to_same_source_continuation':False,
      'all_runtime_arithmetic_deployment_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
