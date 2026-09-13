"""Common pending boundary for the coherent machine TuneState product.

Shipping does not commit tau, sigma, and R_S independently.  Sample k updates
one ``TuneState`` and raises one pending bit; the beginning of the next IMU
sample consumes that complete snapshot.  This layer joins the deployment
machine product to the already-certified finite commit transaction without
allowing scalar histories to be spliced.

Two global compiler histories are retained until FP contraction mode is
qualified.  For each history the carried machine values are converted to one
``finite_tuner_commit.TuneState`` and committed atomically.  MAG/HOLD events do
not execute an IMU boundary and preserve the entire product by identity.

This closes topology/ancestry of the common boundary only.  It deliberately
does not claim target-libm correctness, compiler-mode selection, binary32
correspondence of the exact-real commit formulas, startup attachment, or the
600-step theorem word.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as M
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C

QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_BOUNDARY_V1'


@dataclass(frozen=True)
class Boundary:
    state:M.State
    separate_commit:C.CommitResult|None
    fma_commit:C.CommitResult|None
    consumed:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.state,M.State) or not isinstance(self.consumed,bool):
            raise TypeError('machine TuneState successor and literal consumed flag required')
        if self.qualification!=QUALIFICATION:
            raise ValueError('wrong machine TuneState boundary qualification')
        if self.consumed:
            if not isinstance(self.separate_commit,C.CommitResult) or not isinstance(self.fma_commit,C.CommitResult):
                raise TypeError('consumed common boundary requires both compiler commits')
            if self.state.pending:
                raise ValueError('consumed common pending bit must clear')
        elif self.separate_commit is not None or self.fma_commit is not None:
            raise ValueError('unconsumed boundary cannot carry commit results')


def _tune(state:M.State,mode:str):
    if not isinstance(state,M.State): raise TypeError('coherent machine TuneState required')
    if mode=='separate':
        return C.TuneState(state.tau.separate,state.sigma.separate,state.rs.separate)
    if mode=='fma':
        return C.TuneState(state.tau.fma,state.sigma.fma,state.rs.fma)
    raise ValueError('unknown global compiler mode')


def imu_boundary(state:M.State,cfg:C.CommitConfig,*,live,band_noise_floor_sigma,
                 separate_rs_sqrt_scale=None,fma_rs_sqrt_scale=None,
                 sync_covariance=False,rs_scale=1):
    """Consume one carried whole-TuneState pending snapshot at the next IMU edge."""
    if not isinstance(state,M.State) or not isinstance(cfg,C.CommitConfig):
        raise TypeError('machine TuneState and CommitConfig required')
    if not isinstance(live,bool) or not isinstance(sync_covariance,bool):
        raise TypeError('literal live/sync branches required')
    if not state.pending:
        if separate_rs_sqrt_scale is not None or fma_rs_sqrt_scale is not None:
            raise ValueError('nonpending boundary consumes no R_S sqrt witnesses')
        return Boundary(state,None,None,False)

    sep=C.commit(_tune(state,'separate'),cfg,pending=True,live=live,
                 band_noise_floor_sigma=band_noise_floor_sigma,
                 rs_sqrt_scale=separate_rs_sqrt_scale,
                 sync_covariance=sync_covariance,rs_scale=rs_scale)
    fma=C.commit(_tune(state,'fma'),cfg,pending=True,live=live,
                 band_noise_floor_sigma=band_noise_floor_sigma,
                 rs_sqrt_scale=fma_rs_sqrt_scale,
                 sync_covariance=sync_covariance,rs_scale=rs_scale)
    if sep is None or fma is None:
        raise AssertionError('pending common boundary unexpectedly produced no commit')
    return Boundary(M.State(state.tau,state.sigma,state.rs,False),sep,fma,True)


def mag_or_hold(state:M.State):
    """Asynchronous non-IMU event: literal identity on all tuner machine state."""
    if not isinstance(state,M.State): raise TypeError('coherent machine TuneState required')
    return state


def readiness():
    lower=M.readiness(); commit=C.readiness()
    return {
      'coherent_machine_TuneState_product_consumed':lower['shipping_tau_sigma_RS_persistent_machine_states_composed'],
      'one_common_pending_bit_controls_tau_sigma_RS_transaction':True,
      'separate_global_compiler_track_committed_atomically':True,
      'FMA_global_compiler_track_committed_atomically':True,
      'same_carried_tau_drives_OU_and_pseudoS_cadence':commit['same_tau_drives_OU_and_S_cadence'],
      'same_carried_sigma_drives_stationary_aw':commit['same_sigma_drives_stationary_aw_target'],
      'same_carried_RS_enters_Live_RS_commit':commit['live_RS_from_same_TuneState_and_realized_period'],
      'preLive_pending_consumes_OU_transaction_without_RS_write':True,
      'MAG_and_HOLD_preserve_whole_machine_TuneState_and_pending_by_identity':True,
      'next_boundary_common_machine_commit_attached':True,
      'compiler_FP_contraction_mode_qualified':False,
      'exact_commit_to_binary32_shipping_correspondence_closed':False,
      'target_libm_correspondence_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
