#!/usr/bin/env python3
"""Reset-complete execution facade for the literal complete-SEA3 word.

The base assembler already packs the exact shipping H18/A21 prediction and
Joseph matrices.  This facade inserts the missing shipping covariance reset
immediately after every accepted S=0, accelerometer and magnetometer correction.
The correction injection dtheta is supplied by the same source/event word; it is
never invented here.  The reset primitive is nonsingular for every finite
injection, so no small-angle gate is added.
"""
from __future__ import annotations
import json
from pathlib import Path
import ou3_sea3_full_normal_live_word as BASE
import ou3_sea3_full_word_reset_congruence as RESET

SCHEMA=1
QUALIFICATION="OU3_SEA3_LITERAL_FULL_NORMAL_LIVE_WORD_WITH_SHIPPING_RESETS"
CANONICAL_SOURCE="COMPLETE_SEA3_NORMAL_LIVE_WORD"
DEFAULT_DOMAIN=BASE.DEFAULT_DOMAIN

# Re-export the exact matrix packers/state type so this remains one literal API.
LiteralWordState=BASE.LiteralWordState
state_dimension=BASE.state_dimension
initialize_word=BASE.initialize_word
pack_prediction=BASE.pack_prediction
H_S_zero=BASE.H_S_zero
R_S_zero=BASE.R_S_zero
H_accelerometer=BASE.H_accelerometer
H_magnetometer=BASE.H_magnetometer
diagonal_R=BASE.diagonal_R
aw_floor_increment=BASE.aw_floor_increment
BACKEND=BASE.BACKEND
OFF_TH,OFF_BG,OFF_V,OFF_P,OFF_S,OFF_AW,OFF_BA=(BASE.OFF_TH,BASE.OFF_BG,BASE.OFF_V,BASE.OFF_P,BASE.OFF_S,BASE.OFF_AW,BASE.OFF_BA)

def _reset(word:LiteralWordState,dtheta,label:str):
    if dtheta is None or len(dtheta)!=3: raise ValueError(f"{label} correction requires same-word dtheta")
    G=RESET.reset_matrix(dtheta,word.dimension)
    RESET.apply_joint_reset(word.riccati,G)
    word.event_log.append(f"left_reset_{label}")

def apply_imu_sample(word:LiteralWordState,*,F,Q,f_cog_body,R_wb,Racc,due_S,rs_std_xyz,Delta_aw,S_reset_dtheta=None,acc_reset_dtheta=None):
    """Execute prediction/floor/S-reset/accelerometer-reset in shipping order."""
    BACKEND.predict(word.riccati,F,Q); word.event_log.append("prediction")
    if Delta_aw is not None:
        BACKEND.add_psd_floor(word.riccati,aw_floor_increment(word.mode,Delta_aw)); word.aw_floor_applications+=1; word.event_log.append("aw_floor")
    if due_S:
        if rs_std_xyz is None: raise ValueError("due S update requires actual applied R_S")
        BACKEND.joseph_measurement(word.riccati,H_S_zero(word.mode),R_S_zero(rs_std_xyz)); word.S_updates+=1; word.event_log.append("S_zero")
        _reset(word,S_reset_dtheta,"S")
    BACKEND.joseph_measurement(word.riccati,H_accelerometer(word.mode,f_cog_body,R_wb),Racc); word.accel_updates+=1; word.event_log.append("accelerometer")
    _reset(word,acc_reset_dtheta,"acc")
    word.imu_samples+=1

def apply_magnetometer(word:LiteralWordState,*,m_body,Rmag,reset_dtheta):
    BACKEND.joseph_measurement(word.riccati,H_magnetometer(word.mode,m_body),Rmag); word.mag_updates+=1; word.event_log.append("magnetometer")
    _reset(word,reset_dtheta,"mag")

def certify_literal_endpoint(word:LiteralWordState,delta=BASE.USEFUL_GATE): return BASE.certify_literal_endpoint(word,delta)

def _self_test(mode):
    from ou3_interval import Interval,matrix_identity,matrix_point
    n=state_dimension(mode); P0=matrix_point([[2.0 if i==j else 0.0 for j in range(n)] for i in range(n)]); w=initialize_word(mode,P0)
    Faa=matrix_identity(6); Qaa=matrix_point([[0.02 if i==j else 0.0 for j in range(6)] for i in range(6)]); Fll=matrix_identity(12)
    h=0.005
    for a in range(3): Fll[a][9+a]=Interval.point(h); Fll[3+a][a]=Interval.point(h); Fll[6+a][3+a]=Interval.point(h)
    Qll=matrix_point([[0.01 if i==j else 0.0 for j in range(12)] for i in range(12)])
    if mode=="A": F,Q=pack_prediction(mode,Faa,Qaa,Fll,Qll,phi_ba=Interval.point(0.999999),Q_ba=matrix_point([[1e-8 if i==j else 0.0 for j in range(3)] for i in range(3)]))
    else: F,Q=pack_prediction(mode,Faa,Qaa,Fll,Qll)
    d=[Interval.point(0.01),Interval.point(-0.005),Interval.point(0.002)]
    apply_imu_sample(w,F=F,Q=Q,f_cog_body=[Interval.point(0),Interval.point(0),Interval.point(-9.80665)],R_wb=matrix_identity(3),Racc=diagonal_R([0.2]*3),due_S=True,rs_std_xyz=[Interval.point(0.72),Interval.point(0.72),Interval.point(1.0)],Delta_aw=matrix_point([[0.001 if i==j else 0.0 for j in range(3)] for i in range(3)]),S_reset_dtheta=d,acc_reset_dtheta=d)
    apply_magnetometer(w,m_body=[Interval.point(20),Interval.point(0),Interval.point(40)],Rmag=diagonal_R([0.3]*3),reset_dtheta=d)
    return {"mode":mode,"decomposition_identity_enclosed":BACKEND.decomposition_identity_enclosed(w.riccati),"event_log":w.event_log,
            "all_three_measurements_followed_immediately_by_reset":w.event_log==["prediction","aw_floor","S_zero","left_reset_S","accelerometer","left_reset_acc","magnetometer","left_reset_mag"]}

def build(domain_path:Path=DEFAULT_DOMAIN):
    b=BASE.build(domain_path); bf=BASE.validate(b); rf=RESET.validate()
    if bf or rf: raise RuntimeError(f"reset-complete literal prerequisites: base={bf}, reset={rf}")
    H=_self_test("H"); A=_self_test("A"); parity=RESET.source_parity()
    out=dict(b); out["schema"]=SCHEMA; out["qualification"]=QUALIFICATION; out["base_qualification"]=b["qualification"]
    # BASE is the literal complete-SEA3 word facade; this tag records that source
    # identity for downstream composition and does not create or replace a source.
    out["canonical_source"]=CANONICAL_SOURCE
    out["shipping_reset_source_parity"]=parity; out["shipping_reset_source_parity_pass"]=all(parity.values())
    out["reset_injection_supplied_by_same_source_word"]=True; out["reset_small_angle_bound_required"]=False
    out["S_Joseph_immediately_followed_by_left_reset"]=True; out["accelerometer_Joseph_immediately_followed_by_left_reset"]=True; out["magnetometer_Joseph_immediately_followed_by_left_reset"]=True
    out["literal_reset_execution_complete"]=True; out["H_reset_execution_self_test"]=H; out["A_reset_execution_self_test"]=A
    api=dict(b["numeric_execution_api"]); api["left_error_reset"]="RESET.reset_matrix + RESET.apply_joint_reset immediately after every accepted Joseph correction"; out["numeric_execution_api"]=api
    return out

def validate(d):
    probe=dict(d); probe["schema"]=BASE.SCHEMA; probe["qualification"]=d.get("base_qualification"); f=BASE.validate(probe)
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION: f.append("facade schema/qualification")
    if d.get("canonical_source")!=CANONICAL_SOURCE: f.append("canonical source")
    for k in ("shipping_reset_source_parity_pass","reset_injection_supplied_by_same_source_word","S_Joseph_immediately_followed_by_left_reset","accelerometer_Joseph_immediately_followed_by_left_reset","magnetometer_Joseph_immediately_followed_by_left_reset","literal_reset_execution_complete"):
        if d.get(k) is not True: f.append(k)
    if d.get("reset_small_angle_bound_required") is not False: f.append("reset small-angle gate introduced")
    for mode in ("H","A"):
        s=d.get(f"{mode}_reset_execution_self_test",{})
        if s.get("decomposition_identity_enclosed") is not True or s.get("all_three_measurements_followed_immediately_by_reset") is not True: f.append(f"{mode} reset self-test")
    return list(dict.fromkeys(f))

def main():
    d=build(); f=validate(d); print(json.dumps({"qualification":d["qualification"],"canonical_source":d["canonical_source"],"H_log":d["H_reset_execution_self_test"]["event_log"],"A_log":d["A_reset_execution_self_test"]["event_log"],"failures":f},indent=2)); return 0 if not f else 2
if __name__=="__main__": raise SystemExit(main())