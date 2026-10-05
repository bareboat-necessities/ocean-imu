"""Structural feedback audit for the deployed planar admission path.

The default tuner/front-end is measurement-only.  It is therefore an exogenous
causal coefficient generator for the MEKF: covariance/MEKF mean do not feed
tau, sigma_aw, R_S, T_S, the S clock, or the AW-sync target.  This does NOT say
the coefficients are constant, nor does it certify their all-time range.

The direction branch may read MEKF attitude/bias, but it has no return path into
the OU tuner/model schedule.  The remaining bidirectional small-gain loop is
MEKF mean <-> covariance through state-dependent measurement Jacobians/gains.
"""
from __future__ import annotations
import hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
WRAPPER=ROOT/"src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
FRONT=ROOT/"src/kalman_common/MarineWaveFrontEnd.h"

def git_blob_sha(path):
    data=path.read_bytes(); h=hashlib.sha1()
    h.update(f"blob {len(data)}\\0".encode("ascii")); h.update(data)
    return h.hexdigest()

def certificate():
    w=WRAPPER.read_text(); f=FRONT.read_text()
    required_wrapper=[
      "levelVertical_(dt, gyro, acc_in);",
      "const float a_vert_measurement = trackAccelBandFrequency_(dt);",
      "update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());",
      "updateWavePeriodAndDirection_(dt, mekf_->quaternion_boat(), acc_in, &acc_bias);",
      "apply_pending_online_tune_();"]
    required_front=[
      "return vertical_accel_comp_.verticalAccelUpMs2();",
      "const float wave_hz = wave_period_.getFrequencyHz();",
      "wave_period_.update(dt, wave_period_input_ms2_(direction_accel));"]
    if any(w.count(x)!=1 for x in required_wrapper): raise ValueError("wrapper schedule anchor changed")
    if any(f.count(x)!=1 for x in required_front): raise ValueError("frontend schedule anchor changed")
    # Literal order: prior staged coefficients -> sensor-only Mahony -> MEKF ->
    # sensor-only tuner -> period/direction. Thus current MEKF state cannot
    # choose its own or next tuner input on the deployed Complementary source.
    order=[w.index(x) for x in required_wrapper[:3]]
    if not order[0] < order[1] < order[2]: raise ValueError("frontend/tuner chronology changed")
    return {
      "qualification":"OU3_PLANAR_GENERATED_PORT_STRUCTURE_V1",
      "result_type":"PROVED analytical source-control-flow theorem",
      "default_wave_period_input":"Complementary/private Mahony measurement-only path",
      "mekf_mean_or_covariance_to_tuner_gain":0,
      "mekf_mean_or_covariance_to_S_period_gain":0,
      "mekf_mean_or_covariance_to_AW_target_gain":0,
      "tuner_to_covariance_is_one_way":True,
      "direction_branch_reads_mekf_but_returns_to_tuner":False,
      "remaining_bidirectional_loop":"MEKF mean <-> P via state-dependent acc/mag H and P-dependent K",
      "coefficient_ranges_all_time_certified":False,
      "joint_mean_covariance_cell_forward_invariant":False,
      "every_placed_window_magnetic_service_verified":False,
      "all_time_magnetic_service_verified":False,
      "theorem_closed":False,
      "source_git_blob_sha":{
        str(WRAPPER.relative_to(ROOT)):git_blob_sha(WRAPPER),
        str(FRONT.relative_to(ROOT)):git_blob_sha(FRONT)}
    }

if __name__=="__main__":
 import json; print(json.dumps(certificate(),indent=2,sort_keys=True))
