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
PROXY=ROOT/"src/kalman_common/ProxyStartupFusion.h"
MAG_COMMON=ROOT/"src/kalman_common/MagneticStartupCommon.h"
HARD_IRON=ROOT/"src/tuner/ContinuousMagHardIronEstimator.h"
MAG_TUNER=ROOT/"src/tuner/MagAutoTuner.h"
DEFAULTS=ROOT/"src/kalman_common/SeaStateFusionDefaults.h"


def locked_live_reference_audit():
    """Source anchors for app:locked-live-causal-reference's control-flow proof.

    This checks the proved implementation, not arbitrary edits or global
    branch retention. In particular the Live gravity gate still reads MEKF.
    """
    p, m, h = PROXY.read_text(), MAG_COMMON.read_text(), HARD_IRON.read_text()
    defaults=DEFAULTS.read_text()
    if ('bool mag_estimate_hard_iron = false;' not in p
            or 'MAG_HI_MAX_BIAS_FRACTION    = 0.35f;' not in defaults
            or 'MAG_HI_APPLY_FRACTION       = 1.0f;' not in defaults):
        raise ValueError('default reference amplitude qualification changed')
    begin=p.index('void maybeRefineMagReference_(')
    end=p.index('// Hand the bootstrap attitude',begin)
    refinement=p[begin:end]
    for anchor in ['if (!mag_ref_set_) return;',
                   'if (stage_ != Stage::Live) return;',
                   'observation ? yawRemovedBoatQuat(observation->proxy_bw)',
                   ': impl_.startupProxyTiltQuat();',
                   'mag_body_ned - mag_hard_iron_body_uT_',
                   'setMagWorldRef_(mag_world_ref_uT);',
                   'writeMekfYaw_(wrapPi(-mag_tilt_yaw_rad), mag_tilt_yaw_rad);']:
        if anchor not in refinement:
            raise ValueError('reference refinement control flow changed')
    if 'attitudeReferenceQuat_()' in refinement or 'gravity_gate_' in refinement:
        raise ValueError('refinement acquired an unproved MEKF/gate dependency')
    for anchor in ['gravity_gate_.step(attitudeReferenceQuat_()',
                   'if (!mag_ref_set_) {',
                   'observation ? yawRemovedBoatQuat(observation->proxy_bw) : impl_.startupProxyTiltQuat(),',
                   'maybeRefineMagReference_(mag_body_ned, observation, sample_t);',
                   'maybeApplyContinuousHardIron_();',
                   'last_mag_applied_ = correct(mag_body_ned - mag_hard_iron_body_uT_);']:
        if anchor not in p:
            raise ValueError('locked Live magnetic chronology changed')
    for anchor in ['estimator.update(dt_mag, q_tilt_bw, mag_body_uT);',
                   'anchor_world_ref_uT = current_world_ref_uT;',
                   'applied_body_uT + alpha * (target - applied_body_uT)',
                   'estimator.levelReferenceForBias(applied, level_new)',
                   'estimator.levelReferenceForBias(anchor_bias_body_uT, level_anchor)',
                   'anchor_world_ref_uT.x() + (h_new - h_anchor)',
                   'anchor_world_ref_uT.z() + (level_new.z() - level_anchor.z())']:
        if anchor not in m:
            raise ValueError('paired hard-iron/reference map changed')
    for anchor in ['(level_mag_sum_ - rot_sum_ * applied_bias_body_uT.cast<double>())',
                   'weight_sum_ *= lambda;', 'rot_sum_ *= lambda;',
                   'weight_sum_ += 1.0;', 'rot_sum_ += R;',
                   'b.norm() > double(cfg_.max_bias_fraction) * field_scale']:
        if anchor not in h:
            raise ValueError('hard-iron statistics or accepted-bias bound changed')
    return {
        'proof': 'app:locked-live-causal-reference',
        'locked_Live_reference_reverse_MEKF_port': 0,
        'locked_Live_corrected_mag_reverse_MEKF_port': 0,
        'same_input_finite_MEKF_root_changes_preserve_reference_execution': True,
        'timed_input_endogenous_transport_covered': False,
        'scope': 'legacy updateMag interface only; timed magnetic preprocessing adds gyro-bias/covariance ports requiring separate bounds. Magnetically locked Live continuation, inherited auxiliary state and delivered history fixed; default Complementary path, no Cold/reacquisition/external reconfiguration',
        'live_gravity_gate_reads_MEKF': True,
        'live_gravity_gate_controls_locked_refinement_or_slew': False,
        'refinement_yaw_write_remains_actual_estimator_jump': True,
        'source_reference_derivatives_discarded': False,
        'inherited_auxiliary_root_correlations_discarded': False,
        'actual_correlated_lift_map': 'M_lift=M+R_a*Lambda if da_0=Lambda*v_0+zeta is proved; otherwise retain augmented lift',
        'canonical_reference_map_global_Lipschitz_constant': 1,
        'reference_increment_bound': '|B_ref-B_anchor| <= |A (b-b_anchor)| <= |b-b_anchor| at the SAME statistics',
        'default_real_reference_norm_bound': '(1+7/20)*M_mag after the actual acquisition/refinement, startup hard iron disabled and initial applied offset zero',
        'float_reference_accumulation_qualified': False,
        'uniform_complete_net_work_domination': False,
    }

def git_blob_sha(path):
    data=path.read_bytes(); h=hashlib.sha1()
    h.update(f"blob {len(data)}\0".encode("ascii")); h.update(data)
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
      "locked_live_reference":locked_live_reference_audit(),
      "coefficient_ranges_all_time_certified":False,
      "joint_mean_covariance_cell_forward_invariant":False,
      "every_placed_window_magnetic_service_verified":False,
      "all_time_magnetic_service_verified":False,
      "theorem_closed":False,
      "source_git_blob_sha":{
        str(FRONT.relative_to(ROOT)):git_blob_sha(FRONT),
        str(WRAPPER.relative_to(ROOT)):git_blob_sha(WRAPPER),
        **{str(p.relative_to(ROOT)):git_blob_sha(p)
           for p in (PROXY,MAG_COMMON,HARD_IRON,MAG_TUNER,DEFAULTS)}}
    }

if __name__=="__main__":
 import json; print(json.dumps(certificate(),indent=2,sort_keys=True))
