"""Binary32 common pending commit for the coherent machine TuneState.

All three candidate scalars cross a single pending boundary.  The existing
binary32 operation graph, not the exact-real commit formulas, now produces the
applied tau, Sigma_aw, pseudo period and Live R_S.  Fractions in CommitResult
are exact representations of these rounded values, not recomputed real squares.
Each global compiler mode supplies its own explicit machine band-noise floor;
there is no fallback to the exact shadow's value.

The source of those floor witnesses, Eigen execution, finite/nonfinite branch
coverage and injection of active-parameter displacements into Live dynamics
remain open.  This scalar composition does not authorize storage.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as M
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary_commit as BOUND

QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_BOUNDARY_V2'


def runtime_commit(mode:BOUND.ModeCommit,*,sync_covariance=False):
    """Package already-rounded scalar outputs without executing real formulas."""
    if not isinstance(mode,BOUND.ModeCommit): raise TypeError('binary32 ModeCommit required')
    if not isinstance(sync_covariance,bool): raise TypeError('literal covariance sync branch required')
    def diag(v): return tuple(tuple(v[i] if i==j else 0 for j in range(3)) for i in range(3))
    sigma=diag(mode.aw_covariance_diag)
    rs=None if mode.rs is None else diag(mode.rs.covariance_diag)
    return C.CommitResult(mode.tau_command,sigma,mode.pseudo_period,rs,
                          sigma if sync_covariance else None,False)


@dataclass(frozen=True)
class Boundary:
    state:M.State
    separate_commit:C.CommitResult|None
    fma_commit:C.CommitResult|None
    consumed:bool
    arithmetic:BOUND.Result
    sync_covariance:bool=False
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.state,M.State) or not isinstance(self.consumed,bool):
            raise TypeError('machine TuneState successor and literal consumed flag required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine TuneState boundary qualification')
        if not isinstance(self.arithmetic,BOUND.Result) or not isinstance(self.sync_covariance,bool):
            raise TypeError('retained binary32 transaction and literal sync branch required')
        if self.state!=self.arithmetic.state or self.consumed!=self.arithmetic.before.pending:
            raise ValueError('common boundary detached from binary32 transaction')
        if self.consumed:
            if self.state.pending: raise ValueError('consumed common pending bit must clear')
            expected=(runtime_commit(self.arithmetic.separate,sync_covariance=self.sync_covariance),
                      runtime_commit(self.arithmetic.fma,sync_covariance=self.sync_covariance))
            if (self.separate_commit,self.fma_commit)!=expected:
                raise ValueError('applied parameters detached from rounded binary32 outputs')
        elif self.separate_commit is not None or self.fma_commit is not None:
            raise ValueError('unconsumed boundary cannot carry commit results')


def imu_boundary(state:M.State,cfg:C.CommitConfig,*,live,
                 separate_band_noise_floor_sigma=None,fma_band_noise_floor_sigma=None,
                 sync_covariance=False):
    """Consume the pending snapshot using each mode's own machine floor."""
    if not isinstance(sync_covariance,bool): raise TypeError('literal sync branch required')
    arithmetic=BOUND.apply(state,cfg,live=live,
        separate_band_noise_floor_sigma=separate_band_noise_floor_sigma,
        fma_band_noise_floor_sigma=fma_band_noise_floor_sigma)
    sep=fma=None
    if state.pending:
        sep=runtime_commit(arithmetic.separate,sync_covariance=sync_covariance)
        fma=runtime_commit(arithmetic.fma,sync_covariance=sync_covariance)
    return Boundary(arithmetic.state,sep,fma,state.pending,arithmetic,sync_covariance)


def mag_or_hold(state:M.State):
    """Asynchronous non-IMU event: literal identity on all tuner machine state."""
    if not isinstance(state,M.State): raise TypeError('coherent machine TuneState required')
    return state


def readiness():
    lower=M.readiness(); commit=BOUND.readiness()
    return {
      'coherent_machine_TuneState_product_consumed':lower['shipping_tau_sigma_RS_persistent_machine_states_composed'],
      'one_common_pending_bit_controls_tau_sigma_RS_transaction':True,
      'separate_global_compiler_track_committed_atomically':True,
      'FMA_global_compiler_track_committed_atomically':True,
      'same_carried_tau_drives_OU_and_pseudoS_cadence':True,
      'same_carried_sigma_drives_stationary_aw':True,
      'same_carried_RS_enters_Live_RS_commit':True,
      'preLive_pending_consumes_OU_transaction_without_RS_write':True,
      'MAG_and_HOLD_preserve_whole_machine_TuneState_and_pending_by_identity':True,
      'next_boundary_common_machine_commit_attached':True,
      'applied_parameters_packaged_from_binary32_operation_graph':True,
      'mode_specific_machine_noise_floors_required_without_shadow_fallback':True,
      'binary32_scalar_boundary_source_shape_matches':commit['shipping_pending_common_TuneState_boundary_source_shape_matches'],
      'compiler_FP_contraction_mode_qualified':False,
      'exact_commit_to_binary32_shipping_correspondence_closed':False,
      'target_libm_correspondence_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
