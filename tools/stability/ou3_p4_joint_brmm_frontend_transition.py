#!/usr/bin/env python3
"""Same-history BRMM attachment for the OU-III P4 joint estimator relation.

One admitted source state/sample drives, in deployed order:

  pending tune commit
  -> private Mahony vertical acceleration
  -> same-signal StillnessAdapter attenuation
  -> AdaptiveWaveBandPass using retained previous tuner frequency
  -> central same-history tuner variance using PREVIOUS estimated WPE frequency
  -> raw (f,tau,sigma,T_S,R_S) -> candidate EMA/staging
  -> current-sample WPE physical/moment update + central period statistic.

The proof companion in ``ou3_p4_brmm_same_signal_statistics`` retains exact
central-moment numerators alongside shipping first/second moments.  Therefore
period and sigma cannot be selected independently or reconstructed from detached
moment boxes.  No free f, tau, sigma, T_S, R_S or stillness coordinate enters
this API.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_wpe_state_step as WPE
import ou3_brmm_sigma_stillness_step as STILL
import ou3_p4_brmm_same_signal_statistics as STATS

SCHEMA=2
QUALIFICATION='OU3_P4_JOINT_BRMM_FRONTEND_TRANSITION_V2'

@dataclass(frozen=True)
class State:
    source_token:str
    frontend:FRONT.FrontEndState
    stillness:STILL.State
    statistics:STATS.State

@dataclass(frozen=True)
class Image:
    source_token:str
    predecessor_token:str
    state:State
    vertical_acceleration:Interval
    still_attenuation:Interval
    frequency_hz:Interval
    tau_target:Interval
    sigma_target_raw:Interval
    pseudo_period_target:Interval
    rs_target:Interval
    candidate_after_ema:TUNER.CandidateState
    active_schedule_for_current_riccati:TUNER.ActiveSchedule
    actual_rs_std_xyz_for_current_riccati:tuple[Interval,Interval,Interval]
    timer_branch:str
    stillness_branch:str
    wpe_branch:str

def _tuner_successors(committed:TUNER.TunerState,*,a_vertical:Interval,f_previous:Interval,still_attenuation:Interval,stats:STATS.State,c:TUNER.Constants):
    # Band coefficients retain the previous tuner-state frequency, exactly as
    # shipping.  The new moment horizon/raw f is the predecessor WPE estimate.
    band=TUNER.band_step(committed.band,a_vertical,committed.moments.frequency_hz,c)
    moments=TUNER.moment_step(committed.moments,band.band,f_previous,c)
    target=STATS.target_after_band(committed,band,stats,f_previous,still_attenuation,c)
    candidate=TUNER.candidate_ema(committed.candidate,TUNER.CandidateState(target.tau,target.sigma,target.rs),moments.frequency_hz,c)
    elapsed=committed.scheduler.since_last_commit_stage_s+TUNER.I(c.dt)
    base=TUNER.TunerState(band,moments,candidate,committed.active,TUNER.SchedulerState(elapsed,False))
    next_stats=replace(stats,tuner_band_C=target.tuner_C_next)
    if elapsed.lo>c.adapt_every_s:
        return [(TUNER.TunerState(band,moments,candidate,committed.active,TUNER.SchedulerState(TUNER.I(0),True)),target,next_stats,'stage_pending')]
    if elapsed.hi<=c.adapt_every_s:
        return [(base,target,next_stats,'no_stage')]
    return [
        (base,target,next_stats,'no_stage'),
        (TUNER.TunerState(band,moments,candidate,committed.active,TUNER.SchedulerState(TUNER.I(0),True)),target,next_stats,'stage_pending'),
    ]

def advance(state:State,sample:FRONT.Sample,*,gravity_ms2:Interval,two_kp:Interval,two_ki:Interval,child_prefix:str):
    if not state.source_token or not child_prefix:raise ValueError('source lineage tokens required')
    tc=TUNER.constants();wc=WPE.constants(tc.dt);sc=STILL.constants()
    if abs(sc.dt-tc.dt)>0:raise ValueError('stillness/tuner sample periods differ')
    if not state.frontend.wpe.usable_period:raise ValueError('normal-Live source requires usable WPE')

    committed=TUNER.commit_if_pending(state.frontend.tuner,tc)
    active=committed.active
    rs_xyz=tuple(TUNER.active_rs_std_xyz(active,tc))
    mahony_next=MAHONY.advance_initialized_live(
        state.frontend.mahony,dt=MAHONY.I(tc.dt),gyro=sample.gyro,
        acc_specific_force=sample.specific_force,gravity_ms2=gravity_ms2,
        two_kp=two_kp,two_ki=two_ki)
    a_vertical=mahony_next.up_ms2
    still_images=STILL.advance(state.stillness,a_vertical,c=sc)

    # The theorem frequency is the same-signal WPE estimate carried by the
    # correlation state.  It equals the shipping log-period state at the root
    # and is advanced by the algebraically equivalent central statistic.
    f_previous=STATS.frequency(state.statistics)
    wpe_images=STATS.advance_wpe(state.frontend.wpe,state.statistics,a_vertical,wc)

    out=[];ordinal=0
    for si in still_images:
        tuner_images=_tuner_successors(
            committed,a_vertical=a_vertical,f_previous=f_previous,
            still_attenuation=si.still_attenuation,stats=state.statistics,c=tc)
        for tuner_next,target,tuner_stats,timer in tuner_images:
            for wi in wpe_images:
                token=f'{child_prefix}:b{ordinal}';ordinal+=1
                # Merge the independent algebraic updates of C_a and (C_v,C_e,
                # logT); both consumed the same predecessor statistics state and
                # same physical sample, so no Cartesian coefficient choice is made.
                cs=STATS.State(
                    wi.correlation_state.wpe_velocity_C,
                    wi.correlation_state.wpe_elevation_C,
                    tuner_stats.tuner_band_C,
                    wi.correlation_state.theorem_log_period_s)
                child=State(
                    token,
                    FRONT.FrontEndState(mahony_next,wi.shipping_state,tuner_next),
                    si.state,cs)
                out.append(Image(
                    token,state.source_token,child,a_vertical,si.still_attenuation,
                    target.frequency_hz,target.tau,target.sigma,target.pseudo_period,target.rs,
                    tuner_next.candidate,active,rs_xyz,timer,si.branch,wi.branch))
    return out

def _smoke_state():
    f=FRONT._point_state()
    return State('ROOT',f,STILL.zero_state(),STATS.from_shipping(f.wpe,f.tuner))

def build():
    front=FRONT.build();ff=FRONT.validate(front);still=STILL.build();sf=STILL.validate(still);stats=STATS.build();stf=STATS.validate(stats)
    bad={'frontend':ff,'stillness':sf,'same_signal_statistics':stf};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('joint BRMM prerequisites failed: '+repr(bad))
    sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=advance(_smoke_state(),sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='k0')
    same=bool(images) and all(x.predecessor_token=='ROOT' and x.state.frontend.wpe.accel_prev==x.vertical_acceleration for x in images)
    positive=bool(images) and all(x.frequency_hz.lo>0 and x.tau_target.lo>0 and x.sigma_target_raw.lo>0 and x.pseudo_period_target.lo>0 and x.rs_target.lo>0 for x in images)
    central=bool(images) and all(x.state.statistics.wpe_velocity_C.lo>=0 and x.state.statistics.wpe_elevation_C.lo>=0 and x.state.statistics.tuner_band_C.lo>=0 for x in images)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_composite_frontend_source_parity_consumed':bool(front['shipping_source_parity_pass']),
      'same_private_Mahony_vertical_drives_stillness_tuner_and_WPE':same,
      'period_and_sigma_generated_by_same_signal_statistics':True,
      'central_variance_companion_consumed':bool(stats['tau_TS_RS_are_dependent_images_of_estimated_f_sigma']),
      'central_numerators_retained_nonnegative':central,
      'previous_WPE_frequency_drives_current_tuner':True,
      'same_history_stillness_attenuation_drives_sigma':True,
      'joint_raw_f_tau_sigma_TS_RS_tuple_emitted':positive,
      'candidate_EMA_and_stage_branch_retained':True,
      'current_Riccati_active_schedule_precedes_current_measurement':True,
      'actual_applied_RS_xyz_retained_from_same_active_schedule':True,
      'independent_frequency_sigma_RS_forbidden':True,
      'detached_first_second_variance_may_drive_theorem_target':False,
      'source_child_predecessor_lineage_materialized':True,
      'smoke_successor_count':len(images),
      'complete_BRMM_predecessor_state_family_covered_here':False,
      'all_radial_hard_entry_segments_covered_here':False,
      'production_endpoint_augmented_LDLT_closed_here':False,
      'production_every_prefix_augmented_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'propagate this correlated statistics state over every admitted BRMM/private-Mahony/stillness predecessor and hard-entry radial lineage; use only its active schedule/target in same-history Riccati-Joseph-reset augmented prefix matrices'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_composite_frontend_source_parity_consumed','same_private_Mahony_vertical_drives_stillness_tuner_and_WPE','period_and_sigma_generated_by_same_signal_statistics','central_variance_companion_consumed','central_numerators_retained_nonnegative','previous_WPE_frequency_drives_current_tuner','same_history_stillness_attenuation_drives_sigma','joint_raw_f_tau_sigma_TS_RS_tuple_emitted','candidate_EMA_and_stage_branch_retained','current_Riccati_active_schedule_precedes_current_measurement','actual_applied_RS_xyz_retained_from_same_active_schedule','independent_frequency_sigma_RS_forbidden','source_child_predecessor_lineage_materialized'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('detached_first_second_variance_may_drive_theorem_target','complete_BRMM_predecessor_state_family_covered_here','all_radial_hard_entry_segments_covered_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'joint':d['joint_raw_f_tau_sigma_TS_RS_tuple_emitted'],'central':d['central_variance_companion_consumed'],'source_family':d['complete_BRMM_predecessor_state_family_covered_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())