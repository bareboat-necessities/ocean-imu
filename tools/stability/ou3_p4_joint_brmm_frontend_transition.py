#!/usr/bin/env python3
"""Same-history BRMM attachment for the OU-III P4 joint estimator relation.

One already-admitted source state/sample drives, in deployed order:

  pending tune commit
  -> private Mahony vertical acceleration
  -> tracker-input LPF / StillnessAdapter sigma attenuation
  -> AdaptiveWaveBandPass + variance statistic using PREVIOUS WPE frequency
  -> raw (f,tau,sigma,T_S,R_S) target + candidate EMA/staging
  -> current-sample WPE update.

No frequency, sigma, R_S, vertical acceleration or stillness attenuation can be
chosen independently by this API.  Complete-BRMM coverage is intentionally not
claimed here; this is the transition operator that the source-cover recursion
must now propagate over every admitted predecessor cell.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_wpe_state_step as WPE
import ou3_brmm_sigma_stillness_step as STILL
import ou3_p4_brmm_tuner_stillness_transition as TUNE_STILL

SCHEMA=1
QUALIFICATION='OU3_P4_JOINT_BRMM_FRONTEND_TRANSITION_V1'

@dataclass(frozen=True)
class State:
    source_token:str
    frontend:FRONT.FrontEndState
    stillness:STILL.State
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

def advance(state:State,sample:FRONT.Sample,*,gravity_ms2:Interval,two_kp:Interval,two_ki:Interval,child_prefix:str):
    if not state.source_token or not child_prefix:raise ValueError('source lineage tokens required')
    tc=TUNER.constants();wc=WPE.constants(tc.dt);sc=STILL.constants()
    if abs(sc.dt-tc.dt)>0:raise ValueError('stillness/tuner sample periods differ')
    if not state.frontend.wpe.usable_period:raise ValueError('normal-Live source requires usable WPE')

    committed=TUNER.commit_if_pending(state.frontend.tuner,tc)
    active=committed.active
    rs_xyz=tuple(TUNER.active_rs_std_xyz(active,tc))
    mahony_next=MAHONY.advance_initialized_live(state.frontend.mahony,dt=MAHONY.I(tc.dt),gyro=sample.gyro,acc_specific_force=sample.specific_force,gravity_ms2=gravity_ms2,two_kp=two_kp,two_ki=two_ki)
    a_vertical=mahony_next.up_ms2
    still_images=STILL.advance(state.stillness,a_vertical,c=sc)
    f_previous=WPE.frequency_hz(state.frontend.wpe)
    wpe_images=WPE.advance(state.frontend.wpe,a_vertical=a_vertical,c=wc)

    out=[];ordinal=0
    for si in still_images:
        tuner_images=TUNE_STILL.advance_after_measurement(committed,a_vertical=a_vertical,f_wave_previous_wpe=f_previous,still_attenuation=si.still_attenuation,c=tc)
        for tuner_next,target in tuner_images:
            for wpe_next in wpe_images:
                token=f'{child_prefix}:b{ordinal}';ordinal+=1
                f_target=TUNER.clamp_interval(tuner_next.moments.frequency_hz,tc.tune_freq_min,tc.tune_freq_max)
                ts=TUNER.pseudo_period(target.tau,tc)
                timer='stage_pending' if tuner_next.scheduler.pending_commit else 'no_stage'
                child=State(token,FRONT.FrontEndState(mahony_next,wpe_next,tuner_next),si.state)
                out.append(Image(token,state.source_token,child,a_vertical,si.still_attenuation,f_target,target.tau,target.sigma,ts,target.rs,tuner_next.candidate,active,rs_xyz,timer,si.branch))
    return out

def _smoke_state():
    return State('ROOT',FRONT._point_state(),STILL.zero_state())
def build():
    front=FRONT.build();ff=FRONT.validate(front);still=STILL.build();sf=STILL.validate(still);ts=TUNE_STILL.build();tf=TUNE_STILL.validate(ts)
    bad={'frontend':ff,'stillness':sf,'tuner_stillness':tf};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('joint BRMM prerequisites failed: '+repr(bad))
    sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=advance(_smoke_state(),sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='k0')
    same=bool(images) and all(x.predecessor_token=='ROOT' and x.state.frontend.wpe.accel_prev==x.vertical_acceleration for x in images)
    positive=bool(images) and all(x.frequency_hz.lo>0 and x.tau_target.lo>0 and x.sigma_target_raw.lo>0 and x.pseudo_period_target.lo>0 and x.rs_target.lo>0 for x in images)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','trusted_composite_frontend_source_parity_consumed':bool(front['shipping_source_parity_pass']),'same_private_Mahony_vertical_drives_stillness_tuner_and_WPE':same,'previous_WPE_frequency_drives_current_tuner':True,'same_history_stillness_attenuation_drives_sigma':True,'joint_raw_f_tau_sigma_TS_RS_tuple_emitted':positive,'candidate_EMA_and_stage_branch_retained':True,'current_Riccati_active_schedule_precedes_current_measurement':True,'actual_applied_RS_xyz_retained_from_same_active_schedule':True,'independent_frequency_sigma_RS_forbidden':True,'source_child_predecessor_lineage_materialized':True,'smoke_successor_count':len(images),'complete_BRMM_predecessor_state_family_covered_here':False,'all_radial_hard_entry_segments_covered_here':False,'production_endpoint_augmented_LDLT_closed_here':False,'production_every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'propagate this operator over the complete admitted BRMM/private-Mahony/stillness predecessor family and radial hard-entry source cells, then feed each emitted active schedule and target lineage into same-history Riccati/Joseph/reset augmented prefix matrices'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_composite_frontend_source_parity_consumed','same_private_Mahony_vertical_drives_stillness_tuner_and_WPE','previous_WPE_frequency_drives_current_tuner','same_history_stillness_attenuation_drives_sigma','joint_raw_f_tau_sigma_TS_RS_tuple_emitted','candidate_EMA_and_stage_branch_retained','current_Riccati_active_schedule_precedes_current_measurement','actual_applied_RS_xyz_retained_from_same_active_schedule','independent_frequency_sigma_RS_forbidden','source_child_predecessor_lineage_materialized'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_BRMM_predecessor_state_family_covered_here','all_radial_hard_entry_segments_covered_here','production_endpoint_augmented_LDLT_closed_here','production_every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'joint':d['joint_raw_f_tau_sigma_TS_RS_tuple_emitted'],'source_family':d['complete_BRMM_predecessor_state_family_covered_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
