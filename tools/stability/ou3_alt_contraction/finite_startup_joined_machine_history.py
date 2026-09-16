"""Construction-to-goLive guard/Mahony/LPF/stillness/TuneState relation.

The lower startup product owns the one exact frontend event. This product adds
arithmetic relations to THAT event, not a replay or a second physical history.
All upstream machine memories start at reset, advance even during Cold, and
are substituted unchanged into goLive and the joined Live Racc word.

This is a conditional finite graph, not a startup reachability certificate.
The ordinary FromTwoVectors seed is represented; the nearly antiparallel SVD,
source admission, timeout/ungauged handoff, clocks and deployment profile remain
open. The shared conditional sample ceiling is bookkeeping, not a startup time bound.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_tunestate_deployment as LOWER
from tools.stability.ou3_alt_contraction import finite_startup_live_machine_tunestate_bridge as GO
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as T
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as VERTICAL
from tools.stability.ou3_alt_contraction import finite_admitted_machine_runtime_config_interleaved_prefix as CONFIG
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SIGMA
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as RUNTIME
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as LPF
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as STILL
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as MSTILL
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR_CONTRACT
from tools.stability.ou3_alt_contraction import finite_machine_startup_core as CORE_INIT

STARTUP_CONFIG_KEYS=frozenset(('guard_cfg','vertical_cfg','wpe_cfg','band_cfg',
    'stats_cfg','bench_noise_sigma','still_cfg','candidate_cfg'))


@dataclass(frozen=True)
class State:
    base: LOWER.State
    runtime: RUNTIME.RuntimeConfig
    deployment_cfg: D.DeploymentConfig
    guard: GUARD.State
    separate_source: VS.State
    fma_source: VS.State
    racc: RACC.State
    last_raw: SENSOR.RawImuSample | None = None
    sensor_history: SENSOR_CONTRACT.History | None = None
    core_construction: CORE_INIT.Construction | None = None

    def __post_init__(self):
        if not isinstance(self.base,LOWER.State) or not isinstance(self.runtime,RUNTIME.RuntimeConfig):
            raise TypeError('startup TuneState product and carried RuntimeConfig required')
        if not isinstance(self.deployment_cfg,D.DeploymentConfig):
            raise TypeError('carried deployment configuration required')
        if not SIGMA.configs_match(self.runtime.candidate_cfg,self.deployment_cfg):
            raise ValueError('startup deployment config detached from carried tuner config')
        if B.rn32(self.runtime.boundary_bench_noise_sigma)!=B.rn32(self.runtime.bench_noise_sigma):
            raise ValueError('startup boundary/candidate bench sigma detached')
        if not isinstance(self.guard,GUARD.State): raise TypeError('persistent machine guard required')
        n=self.base.frontends.samples
        if bool(n)!=(self.last_raw is not None):
            raise ValueError('startup physical packet memory detached from sample count')
        if self.last_raw is not None and not isinstance(self.last_raw,SENSOR.RawImuSample):
            raise TypeError('same last raw physical startup packet required')
        if self.sensor_history is not None:
            if not isinstance(self.sensor_history,SENSOR_CONTRACT.History):
                raise TypeError('persistent commissioned startup sensor history required')
            if B.rn32(self.runtime.vertical_cfg.gravity)!=B.rn32(SENSOR_CONTRACT.G):
                raise ValueError('startup observer gravity detached from selected sensor domain')
            if self.last_raw is not None:
                SENSOR_CONTRACT.check_packet(self.sensor_history,self.last_raw,ordinal=n)
        if self.guard.samples!=n: raise ValueError('startup guard count detached from machine frontend')
        for source in (self.separate_source,self.fma_source):
            if not isinstance(source,VS.State): raise TypeError('persistent machine vertical/stillness source required')
            if source.samples!=n: raise ValueError('startup source count detached from machine frontend')
            if source.lpf.cutoff_hz!=CONFIG.DEFAULT_TRACKER_CUTOFF:
                raise ValueError('startup tracker cutoff detached from shipping constructor')
        if self.separate_source.vertical!=self.fma_source.vertical:
            raise ValueError('startup compiler histories diverged private Mahony')
        if self.racc!=RACC.State():
            raise ValueError('bootstrap drive_mekf=false cannot advance Racc inflation state')
        if self.base.lower.frontend.tuner.stage=='Live':
            raise ValueError('startup product must cross goLive, not accept a Live snapshot')
        if self.core_construction is not None and not isinstance(self.core_construction,CORE_INIT.Construction):
            raise TypeError('source-owned default CORE construction required')

    @property
    def guard_cfg(self): return CONFIG._machine_guard_cfg(self.runtime)


def initial(runtime:RUNTIME.RuntimeConfig,deployment_cfg:D.DeploymentConfig,*,sensor_history=None):
    """Literal default wrapper reset; no caller-provided frontend/observer seed."""
    if not CONFIG._default_tracker_cutoff_source_matches():
        raise RuntimeError('shipping tracker default/reset source shape changed')
    from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
    from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
    from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
    from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
    tune=COMMIT.TuneState(B.rn32(F(11,10)),B.rn32(F(1,100)),B.rn32(F(1,2)))
    frontend=FRONT.State(G.State(),T.State(V.State(),W.WPEState(),BAND.BandState(),
        BAND.StatsState(),LPF.LPFState(),
        STILL.State(),tune,stage='Cold',warmup_sec=F(10)))
    src=VS.State(V.State(),VS.LPFState(cutoff_hz=CONFIG.DEFAULT_TRACKER_CUTOFF),MSTILL.State())
    return State(LOWER.initial(frontend),runtime,deployment_cfg,GUARD.State(),src,src,RACC.State(),
                 sensor_history=sensor_history,core_construction=CORE_INIT.literal_reset())


@dataclass(frozen=True)
class StepResult:
    state: State
    lower: LOWER.StepResult
    guard: GUARD.Result
    separate_source: VS.Result
    fma_source: VS.Result
    def __post_init__(self):
        if self.state.base!=self.lower.state or self.state.guard!=self.guard.state:
            raise ValueError('startup successor detached from same executed event')
        if self.state.separate_source!=self.separate_source.state or self.state.fma_source!=self.fma_source.state:
            raise ValueError('startup successor detached from machine source event')
        if self.separate_source.mahony!=self.fma_source.mahony:
            raise ValueError('startup common Mahony source diverged')
        VS.require_frontend_input(self.separate_source,self.lower.separate_frontend)
        VS.require_frontend_input(self.fma_source,self.lower.fma_frontend)
        if self.lower.machine is not None:
            VS.require_sigma_stillness(self.separate_source,self.lower.separate_sigma_join.machine)
            VS.require_sigma_stillness(self.fma_source,self.lower.fma_sigma_join.machine)


def _with_base(state,base):
    return State(base,state.runtime,state.deployment_cfg,state.guard,
                 state.separate_source,state.fma_source,state.racc,state.last_raw,state.sensor_history,state.core_construction)


def boundary(state:State,**witnesses):
    """Apply the pending transaction; no guard/observer/filter-memory restart."""
    if not isinstance(state,State): raise TypeError('joined startup state required')
    if {'cfg','bench_noise_sigma','sync_covariance','rs_scale'} & set(witnesses):
        raise TypeError('startup boundary configuration is carried, not event-local')
    out=LOWER.boundary(state.base,state.runtime.commit_cfg,
        bench_noise_sigma=state.runtime.boundary_bench_noise_sigma,**witnesses)
    return _with_base(state,out.state),out


def step(state:State,raw,*,dt,machine_dt,machine_gyro_body,machine_acc_body,
         guard_witnesses, separate_source_witnesses, fma_source_witnesses, **witnesses):
    if not isinstance(state,State): raise TypeError('joined startup state required')
    if not isinstance(raw,SENSOR.RawImuSample): raise TypeError('same physical raw startup packet required')
    if state.sensor_history is not None:
        SENSOR_CONTRACT.check_packet(state.sensor_history,raw,ordinal=state.guard.samples+1)
        if F(dt)!=SENSOR_CONTRACT.DT:
            raise ValueError('commissioned startup profile requires canonical 5ms sampling')
        if state.last_raw is not None:
            rate=F(SENSOR_CONTRACT.domain()['true_gyro_bias_rate_upper_rad_s2'])
            delta=tuple(a-b for a,b in zip(raw.physical.gyro_bias,state.last_raw.physical.gyro_bias))
            if SENSOR_CONTRACT.norm2(delta)>(rate*F(dt))**2:
                raise ValueError('startup true gyro bias increment exceeds profile')
    physical=raw.physical
    if physical.time!=state.base.lower.frontend.tuner.time:
        raise ValueError('startup raw packet time detached from sample-entry clock')
    for name in ('history_id','bias_root','bias_family'):
        value=getattr(physical,name,None)
        if not isinstance(value,str) or not value:
            raise TypeError('startup raw packet requires concrete physical/bias history identity')
        if state.last_raw is not None and value!=getattr(state.last_raw.physical,name):
            raise ValueError('startup physical/bias history identity changed between samples')
    if STARTUP_CONFIG_KEYS & set(witnesses) or 'deployment_cfg' in witnesses:
        raise TypeError('startup event cannot override carried runtime configuration')
    h=GUARD.api_scalar(dt,machine_dt,'startup dt')
    gyro=GUARD.api_vec(raw.raw_gyro_body,machine_gyro_body,'startup gyro')
    acc=GUARD.api_vec(raw.raw_accel_body,machine_acc_body,'startup accel')
    forbidden={'raw_gyro','raw_acc','dt'} & set(guard_witnesses)
    if forbidden: raise TypeError('machine guard operands belong to the same raw startup packet')
    guard=GUARD.step(state.guard,state.guard_cfg,raw_gyro=gyro,raw_acc=acc,dt=h,**guard_witnesses)
    if state.sensor_history is not None and (guard.conditioned_acc!=acc or guard.state.weight!=0 or
            guard.removed_rms>state.guard_cfg.engage_lo):
        raise ValueError('commissioned startup profile requires the retained dormant guard branch')
    sep=VERTICAL._source_step(state.separate_source,guard,state.runtime,dt,**separate_source_witnesses)
    fma=VERTICAL._source_step(state.fma_source,guard,state.runtime,dt,**fma_source_witnesses)
    if state.sensor_history is not None:
        cert=SENSOR_CONTRACT.vertical_supply_certificate(state.sensor_history.profile)
        for source in (sep,fma):
            if SENSOR_CONTRACT.norm2(source.state.vertical.q)>cert['computed_quaternion_norm2_upper']:
                raise AssertionError('same scalar Mahony normalization exceeded proved shell')
            if abs(source.band_input)>cert['vertical_abs_upper_after_defined_Mahony_update']:
                raise AssertionError('same Mahony projection exceeded proved vertical bound')
    config={k:getattr(state.runtime,k) for k in STARTUP_CONFIG_KEYS}
    if state.base.lower.frontend.tuner.stage=='Cold':
        config.pop('candidate_cfg')  # shipping returns before the candidate read
    # Exactly one lower event; bind the already-consumed band/sigma operands.
    out=LOWER.step(state.base,raw,dt=dt,deployment_cfg=state.deployment_cfg,**config,**witnesses)
    nxt=State(out.state,state.runtime,state.deployment_cfg,guard.state,sep.state,fma.state,state.racc,raw,state.sensor_history,state.core_construction)
    return StepResult(nxt,out,guard,sep,fma)


@dataclass(frozen=True)
class GoLive:
    startup: State
    lower: GO.Result
    scheduler_before: POST.Scheduler
    scheduler_park: ACTIVE.NextafterParkWitness | None = None
    def __post_init__(self):
        if not isinstance(self.startup,State) or not isinstance(self.lower,GO.Result):
            raise TypeError('joined startup and same goLive result required')
        if self.lower.startup!=self.startup.base:
            raise ValueError('goLive detached from joined startup machine history')
        if self.lower.live.state.racc!=self.startup.racc:
            raise ValueError('goLive restarted startup Racc state')
        if not isinstance(self.scheduler_before,POST.Scheduler):
            raise TypeError('same pre-goLive scheduler required')
        if ACTIVE.retarget_scheduler(self.lower.live.active,self.scheduler_before,park=self.scheduler_park)!=self.lower.live.state.scheduler:
            raise ValueError('goLive scheduler detached from carried predecessor')
        last=self.startup.last_raw
        if last is not None:
            ref=self.lower.live.state.mekf.reference
            if any(getattr(last.physical,k)!=getattr(ref,k) for k in ('history_id','bias_root','bias_family')):
                raise ValueError('goLive physical/bias history detached from startup packets')
        # These identity joins do not prove full startup BRMM/BIAS admission.
        # GO.Result checks machine/WPE/frontends/pending against its startup.


def go_live(state:State,entry,fresh,**witnesses):
    """Conditional quality/timeout handoff; no universal capture inferred."""
    if not isinstance(state,State): raise TypeError('joined startup state required')
    if {'racc','commit_cfg','bench_noise_sigma'} & set(witnesses):
        raise TypeError('goLive cannot replace carried Racc/configuration history')
    go=GO.bridge(state.base,entry,fresh,racc=state.racc,
        commit_cfg=state.runtime.commit_cfg,bench_noise_sigma=state.runtime.bench_noise_sigma,**witnesses)
    return GoLive(state,go,witnesses['scheduler'],witnesses.get('scheduler_park'))



def admitted_live(go:GoLive,magnetic,origin,bias_history,*,
                  separate_scheduler_park=None,fma_scheduler_park=None,**source_witnesses):
    """Construct the joined admitted Live word from the same conditional goLive.

    No caller-selected Live guard, observer, frontend, sigma, Racc or scheduler
    snapshot is accepted. The lower factory checks the magnetic dual-clock
    ancestry and the same admitted BRMM/BIAS origin and bounded ISS histories.
    The admitted-origin quantifier is not proof of startup reachability.
    """
    if not isinstance(go,GoLive): raise TypeError('joined startup goLive result required')
    if 'runtime' in source_witnesses:
        raise TypeError('admitted Live runtime is owned by joined startup')
    from tools.stability.ou3_alt_contraction import finite_startup_live_wpe_tau_bridge as WG
    from tools.stability.ou3_alt_contraction import finite_admitted_startup_wpe_tau_word as ADMIT
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MT
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_scheduler_interleaved_prefix as SCHED
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as PRED
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as MEAS
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_racc_supply_interleaved_prefix as RC
    from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
    # Project a proved subrelation. Do not execute goLive a second time.
    wg=WG.Result(go.lower.live,go.lower.machine.tau,go.lower.wpe,go.startup.base.lower)
    base=ADMIT.build(wg,magnetic,origin,bias_history,runtime=go.startup.runtime,**source_witnesses)
    mt=MT.begin_from_goLive(base,go.lower,go.startup.deployment_cfg)
    pred=SCHED.LOWER.begin(mt)
    sched=SCHED.begin_from_goLive(pred,go.lower,scheduler_before=go.scheduler_before,
        exact_scheduler_park=go.scheduler_park,separate_scheduler_park=separate_scheduler_park,
        fma_scheduler_park=fma_scheduler_park)
    return JOIN.begin_from_goLive(RC.begin(MEAS.begin(PRED.begin(sched))),go)

def readiness():
    return {
        'literal_reset_roots_guard_observer_LPF_stillness_band_and_TuneState':True,
        'Cold_and_postCold_machine_sources_joined_to_one_executed_frontend_event':True,
        'bootstrap_Racc_identity_follows_drive_mekf_false':True,
        'conditional_goLive_preserves_joined_machine_source_history':True,
        'admitted_Live_factory_substitutes_same_startup_without_new_snapshots':True,
        'same_physical_bias_history_identity_and_sample_clock_carried':True,
        'commissioned_sensor_profile_checked_on_same_startup_packets':True,
        'sensor_profile_preserved_by_pending_boundary_and_goLive':True,
        'peak_sensor_checks_prove_temporal_direction_membership':False,
        'history_identity_alone_proves_admission':False,
        'startup_source_membership_and_clock_reachability_closed':False,
        'nearly_antiparallel_seed_SVD_closed':False,
        'ungauged_timeout_machine_handoff_closed':False,
        'source_uniform_complete_600_step_word_qualified':False,
        'storage_search_allowed':False,
        'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
