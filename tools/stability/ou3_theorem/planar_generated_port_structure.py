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
from fractions import Fraction
import hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
WRAPPER=ROOT/"src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
FRONT=ROOT/"src/kalman_common/MarineWaveFrontEnd.h"
PROXY=ROOT/"src/kalman_common/ProxyStartupFusion.h"
WAVE_STARTUP=ROOT/"src/kalman_common/WaveStartup.h"
MAG_COMMON=ROOT/"src/kalman_common/MagneticStartupCommon.h"
HARD_IRON=ROOT/"src/tuner/ContinuousMagHardIronEstimator.h"
GYRO_MAG=ROOT/"src/tuner/GyroMagCalibration.h"
GYRO_MAG_IMPL=ROOT/"src/tuner/GyroMagCalibrationImpl.h"
GYRO_MAG_CPP=ROOT/"src/tuner/GyroMagCalibration.cpp"
MAG_TUNER=ROOT/"src/tuner/MagAutoTuner.h"
DEFAULTS=ROOT/"src/kalman_common/SeaStateFusionDefaults.h"


def wave_handover_reference_factor():
    """Real arithmetic on an accepted seed, before the refined reference takes over."""
    source = WAVE_STARTUP.read_text()
    for anchor in (
            'if (time_ + dt > 35.001) reset();',
            'b.norm() > 0.004f',
            'j.col(k) = y.cross(a.col(k))',
            'integrated_rotation_ += (0.5f*dt)*(old_r+r)',
            'const V field = relative_ * body;',
            'mean_j_ += alpha*dj',
            'const V field = mean_field_ + mean_j_*b;',
            'frame*field, b'):
        if anchor not in source:
            raise ValueError('wave handover reference argument changed: ' + anchor)
    # The decimal upper bounds enclose the represented time and bias guards.
    # The transported sample and its first-order correction are orthogonal.
    upper = Fraction(101, 100)
    if 1 + (Fraction("35.002") * Fraction("0.004001"))**2 >= upper**2:
        raise ValueError('insufficient reference norm bound')
    return upper


def locked_live_reference_audit():
    """Source anchors for app:locked-live-causal-reference's control-flow proof.

    This checks the proved implementation, not arbitrary edits or global
    branch retention. The Live gravity gate reads MEKF during acquisition.
    """
    p, m, h = PROXY.read_text(), MAG_COMMON.read_text(), HARD_IRON.read_text()
    g = GYRO_MAG.read_text() + GYRO_MAG_IMPL.read_text()
    defaults=DEFAULTS.read_text()
    if ('bool mag_estimate_hard_iron = false;' not in p
            or 'bool  mag_gyro_calibration       = true;' not in p
            or 'MAG_HI_MAX_BIAS_FRACTION    = 0.35f;' not in defaults
            or 'MAG_HI_APPLY_FRACTION       = 1.0f;' not in defaults):
        raise ValueError('default reference amplitude qualification changed')
    begin=p.index('void maybeRefineMagReference_(')
    end=p.index('// Hand the bootstrap attitude',begin)
    refinement=p[begin:end]
    for anchor in ['if (!mag_ref_set_) return;',
                   'if (stage_ != Stage::Live) return;',
                   'const Eigen::Quaternionf q_tilt_bw = impl_.startupProxyTiltQuat();',
                   'mag_body_ned - mag_hard_iron_body_uT_',
                   'setMagWorldRef_(mag_world_ref_uT);',
                   'writeMekfYaw_(wrapPi(-mag_tilt_yaw_rad), mag_tilt_yaw_rad);']:
        if anchor not in refinement:
            raise ValueError('reference refinement control flow changed')
    if 'attitudeReferenceQuat_()' in refinement or 'gravity_gate_' in refinement:
        raise ValueError('refinement acquired an unproved MEKF/gate dependency')
    # In Live the attitude read is evaluated only inside the acquisition guard.
    proxy = ' '.join(p.split())
    guarded_read = ('if (stage_ != Stage::Live || (cfg_.with_mag && !mag_ref_set_)) { '
                    'gravity_gate_.step(attitudeReferenceQuat_(),')
    attitude_source = ('Eigen::Quaternionf attitudeReferenceQuat_() const { '
                       'if (stage_ != Stage::Live) { return impl_.startupProxyQuat(); } '
                       'return impl_.mekf().quaternion_boat(); }')
    if (proxy.count('gravity_gate_.step(') != 1 or guarded_read not in proxy
            or attitude_source not in proxy):
        raise ValueError('gravity-gate attitude source or acquisition guard changed')
    for anchor in ['if (!mag_ref_set_) {',
                   'setMagWorldRef_(learned.field);',
                   'if (cfg_.mag_refine_enabled && !mag_refine_done_) return;',
                   'impl_.startupProxyTiltQuat(), mag_body_ned);',
                   'maybeRefineMagReference_(mag_body_ned);',
                   'maybeApplyContinuousHardIron_(mag_body_ned);',
                   'float mag_hi_heading_transport_fraction = 0.f;',
                   'const float fraction = cfg_.mag_hi_heading_transport_fraction;',
                   'const bool transport = fraction > 0.f && fraction <= 1.f && raw.allFinite();',
                   'if (transport) {',
                   'hard_iron_.predict(dt, gyro_body_ned);',
                   'impl_.updateMag(hard_iron_.correct(mag_body_ned));']:
        if anchor not in p:
            raise ValueError('locked Live magnetic chronology changed')
    for anchor in ['estimator.update(dt_mag, q_tilt_bw, mag_body_uT);',
                   'anchor_world_ref_uT = current_world_ref_uT;',
                   'applied_body_uT + alpha * (target - applied_body_uT)',
                   'applied_correction + alpha*(target_correction-applied_correction)',
                   'estimator.levelReferenceForCalibration(correction, startup_body_uT + applied, level_new)',
                   'estimator.levelReferenceForCalibration(anchor_correction, anchor_bias_body_uT, level_anchor)',
                   'return applied_correction * (mag - startup_body_uT - applied_body_uT);',
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
    for anchor in ['(A - Matrix3::Identity()).norm() > 0.25f',
                   'next.segment<3>(3).norm() > 0.35f*mag.norm()',
                   'Matrix3 correction() const { return matrix_().inverse(); }']:
        if anchor not in g:
            raise ValueError('gyro magnetic calibration amplitude qualification changed')
    return {
        'proof': 'app:locked-live-causal-reference',
        'locked_Live_reference_reverse_MEKF_port': 0,
        'locked_Live_corrected_mag_reverse_MEKF_port': 0,
        'same_input_finite_MEKF_root_changes_preserve_reference_execution': True,
        'scope': 'magnetically locked Live continuation, inherited auxiliary state and delivered history fixed; default Complementary path, no Cold/reacquisition/external reconfiguration',
        'live_gravity_gate_MEKF_read_condition': 'with_mag && !mag_ref_set_',
        'locked_Live_gravity_gate_reads_MEKF': False,
        'live_gravity_gate_controls_locked_refinement_or_slew': False,
        'refinement_yaw_write_remains_actual_estimator_jump': True,
        'source_reference_derivatives_discarded': False,
        'inherited_auxiliary_root_correlations_discarded': False,
        'actual_correlated_lift_map': 'M_lift=M+R_a*Lambda if da_0=Lambda*v_0+zeta is proved; otherwise retain augmented lift',
        'canonical_reference_map_global_Lipschitz_constant': 1,
        'reference_increment_bound': '|B_ref-B_anchor| <= ||K-K_anchor||*M_mag + |K*b-K_anchor*b_anchor| at the SAME statistics',
        'default_real_correction_norm_bound': '||K|| <= 4/3 and ||K-I|| <= 1/3 from ||A-I||_F <= 1/4 and convex slew',
        'real_wave_handover_reference_norm_over_delivered_field': str(wave_handover_reference_factor()),
        'default_real_reference_norm_bound': '(9/5)*M_mag after actual acquisition/refinement, startup hard iron disabled, initial K=I and applied offset zero',
        'default_real_corrected_measurement_norm_bound': '(9/5)*M_mag',
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
           for p in (PROXY,WAVE_STARTUP,MAG_COMMON,HARD_IRON,GYRO_MAG,GYRO_MAG_IMPL,GYRO_MAG_CPP,MAG_TUNER,DEFAULTS)}}
    }

if __name__=="__main__":
 import json; print(json.dumps(certificate(),indent=2,sort_keys=True))
