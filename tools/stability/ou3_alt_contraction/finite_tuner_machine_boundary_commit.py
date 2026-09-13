"""Binary32 next-sample commit of the coherent machine TuneState product.

Shipping smooths tau/sigma/R_S on sample k and commits the carried TuneState at
the beginning of sample k+1 when ``online_tune_apply_pending_`` is set.  This
module consumes ``finite_tuner_machine_tunestate_product.State`` directly and
materializes, per coherent global compiler history:

* ``set_aw_time_constant`` command, including its 1e-3 floor;
* tau-scaled pseudo-S period multiply and clamps;
* sigma floor, horizontal scale and per-axis stationary-a_w covariance squares;
* in Live only, the already-source-shaped SpectralMSE R_S commit.

The separate and FMA histories may carry different tau/sigma/R_S values and may
therefore receive different source-qualified band-noise floors.  No exact-real
TuneState is substituted at this boundary.  Eigen setter execution/vectorization
and upstream production of the band-noise-floor float remain explicit open
deployment facts.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as M
from tools.stability.ou3_alt_contraction import finite_tuner_rs_commit_binary32 as RS

FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
MEKF=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
ONE=B.rn32(1); TAU_FLOOR=B.rn32(F(1,1000)); SIGMA_FLOOR=B.rn32(F(1,20))
QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_BOUNDARY_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must compile to an actual binary32 value')
    return q


def _cfg32(cfg:C.CommitConfig):
    if not isinstance(cfg,C.CommitConfig): raise TypeError('CommitConfig required')
    names=('pseudo_tau_ratio','pseudo_period_min','pseudo_period_max','pseudo_fixed_period',
           'min_R_S','max_R_S','S_factor','R_S_x_factor','R_S_y_factor')
    return {n:B.rn32(getattr(cfg,n)) for n in names}


@dataclass(frozen=True)
class ModeCommit:
    mode:str
    stored_tau:F
    stored_sigma:F
    stored_RS:F
    band_noise_floor_sigma:F
    tau_command:F
    pseudo_requested:F
    pseudo_period:F
    sigma_floor:F
    sigma_z:F
    sigma_h:F
    aw_std:tuple
    aw_covariance_diag:tuple
    rs:RS.Commit|None
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('literal compiler mode required')
        for n in ('stored_tau','stored_sigma','stored_RS','band_noise_floor_sigma','tau_command',
                  'pseudo_requested','pseudo_period','sigma_floor','sigma_z','sigma_h'):
            object.__setattr__(self,n,F(getattr(self,n)))
            if not B.is_binary32(getattr(self,n)): raise ValueError(f'{n} is not binary32')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine boundary qualification')
        a=tuple(F(x) for x in self.aw_std); c=tuple(F(x) for x in self.aw_covariance_diag)
        object.__setattr__(self,'aw_std',a); object.__setattr__(self,'aw_covariance_diag',c)
        if not all(B.is_binary32(x) for x in a+c): raise ValueError('a_w commit vectors must be binary32')
        if c != tuple(B.mul(x,x) for x in a): raise ValueError('a_w covariance diagonal detached from std square graph')
        if self.rs is not None and not isinstance(self.rs,RS.Commit): raise TypeError('R_S commit relation required')


@dataclass(frozen=True)
class Result:
    before:M.State
    state:M.State
    separate:ModeCommit|None
    fma:ModeCommit|None
    live:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.before,M.State) or not isinstance(self.state,M.State): raise TypeError('machine TuneState boundary states required')
        if not isinstance(self.live,bool): raise TypeError('literal Live branch required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine boundary qualification')
        if self.before.pending:
            if self.separate is None or self.fma is None: raise ValueError('pending machine TuneState requires both compiler-mode commits')
            if self.state.pending: raise ValueError('shipping pending boundary must clear common pending bit')
            if (self.state.tau,self.state.sigma,self.state.rs)!=(self.before.tau,self.before.sigma,self.before.rs):
                raise ValueError('boundary commit must not advance candidate ledgers')
            if self.live != (self.separate.rs is not None and self.fma.rs is not None):
                raise ValueError('R_S commit presence detached from Live branch')
        else:
            if self.separate is not None or self.fma is not None: raise ValueError('no-pending boundary consumes no commit witnesses')
            if self.state != self.before: raise ValueError('no-pending machine boundary is exact identity')


def _mode(mode,*,tau,sigma,rs,band_noise_floor_sigma,cfg:C.CommitConfig,live:bool):
    q=_cfg32(cfg); bn=_q(band_noise_floor_sigma,f'{mode} band noise floor sigma')
    if bn<0: raise ValueError('band noise floor sigma must be nonnegative')
    t=_q(tau,f'{mode} stored tau'); s=_q(sigma,f'{mode} stored sigma'); r=_q(rs,f'{mode} stored RS')
    tau_cmd=max(TAU_FLOOR,t)
    if cfg.tau_scaled_cadence:
        requested=B.mul(q['pseudo_tau_ratio'],t)
        period=min(max(requested,q['pseudo_period_min']),q['pseudo_period_max'])
    else:
        requested=q['pseudo_fixed_period']; period=requested
    floor=max(SIGMA_FLOOR,bn); sz=max(floor,s); sh=B.mul(sz,q['S_factor'])
    aw=(sh,sh,sz); cov=tuple(B.mul(x,x) for x in aw)
    rs_commit=None
    if live:
        if cfg.cubic_rs_law: raise ValueError('current machine boundary covers deployed SpectralMSE law, not Cubic')
        rs_commit=RS.commit(r,min_RS=q['min_R_S'],max_RS=q['max_R_S'],rs_scale=ONE,
                            x_factor=q['R_S_x_factor'],y_factor=q['R_S_y_factor'])
    return ModeCommit(mode,t,s,r,bn,tau_cmd,requested,period,floor,sz,sh,aw,cov,rs_commit)


def apply(state:M.State,cfg:C.CommitConfig,*,live,separate_band_noise_floor_sigma=None,
          fma_band_noise_floor_sigma=None):
    """Apply the pending TuneState before the next sample reaches the filter."""
    if not isinstance(state,M.State) or not isinstance(cfg,C.CommitConfig): raise TypeError('machine TuneState and CommitConfig required')
    if not isinstance(live,bool): raise TypeError('literal Live branch required')
    if not state.pending:
        if separate_band_noise_floor_sigma is not None or fma_band_noise_floor_sigma is not None:
            raise ValueError('no-pending boundary consumes no band-noise-floor witnesses')
        return Result(state,state,None,None,live)
    if separate_band_noise_floor_sigma is None or fma_band_noise_floor_sigma is None:
        raise ValueError('pending machine boundary requires source-qualified band-noise floor for each compiler track')
    sep=_mode('separate',tau=state.tau.separate,sigma=state.sigma.separate,rs=state.rs.separate,
              band_noise_floor_sigma=separate_band_noise_floor_sigma,cfg=cfg,live=live)
    fma=_mode('fma',tau=state.tau.fma,sigma=state.sigma.fma,rs=state.rs.fma,
              band_noise_floor_sigma=fma_band_noise_floor_sigma,cfg=cfg,live=live)
    nxt=M.State(state.tau,state.sigma,state.rs,False)
    return Result(state,nxt,sep,fma,live)


def _source_shape_matches():
    f=FILTER.read_text(); k=MEKF.read_text()
    return all(x in f for x in (
      'apply_ou_tune_(false);','if (startup_stage_ == StartupStage::Live) {','apply_RS_tune_();',
      'mekf_->set_aw_time_constant(tune_.tau_applied);','const float requested = pseudo_update_tau_ratio_ * tau;',
      'const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());',
      'const float sZ = std::max(sigma_floor, tune_.sigma_applied);','const float sH = sZ * S_factor_;',
      'mekf_->set_aw_stationary_std(aw_std);')) and all(x in k for x in (
      'tau_aw = std::max(T(1e-3), tau_seconds);','Sigma_aw_stat = s.array().square().matrix().asDiagonal();'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_pending_common_TuneState_boundary_source_shape_matches':_source_shape_matches(),
      'boundary_consumes_no_exact_real_TuneState':True,
      'same_mode_tau_sigma_RS_snapshot_committed_together':True,
      'tau_setter_floor_and_pseudo_period_binary32_graph_materialized':True,
      'sigma_floor_horizontal_scale_and_covariance_square_binary32_graph_materialized':True,
      'Live_SpectralMSE_RS_commit_composed_from_same_machine_snapshot':True,
      'preLive_boundary_commits_OU_but_not_RS':True,
      'pending_bit_cleared_without_advancing_candidate_ledgers':True,
      'separate_and_FMA_boundaries_may_consume_distinct_band_noise_floor_values':True,
      'upstream_band_noise_floor_binary32_production_closed':False,
      'Eigen_aw_stationary_std_execution_correspondence_closed':False,
      'Eigen_RS_noise_execution_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
