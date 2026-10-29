#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    The W3dSimulationRunner adapter of the experimental three-axis CARMA INS,
    shared by carma3-sim (CARMA3 variants and the paired TFG-CARMA3 research
    variants) and tests/kalman_tfg_carma/tfg_carma3-sim (the gated TFG-CARMA3
    replay behind doc/kalman_tfg_carma).  The configuration is the one
    tests/kalman_ou_iii/kalman_ou_iii-sim.cpp gives the shipping OU-III; see
    carma3-sim.cpp.
*/

#include <memory>

#include "util/W3dSimCommon.h"
#include "FullAdapters.h"
#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"

namespace carma3_sim {

using Eigen::Quaternionf;
using Eigen::Vector3f;

inline bool add_noise = true;

template<int NHC, bool SH, int NZC, bool SZ,
         class MekfT = ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, NHC, SH, NZC, SZ>>
class CarmaAdapter final : public IW3dFusionAdapter {
public:
    static constexpr float SIGMA_A_RESCALE = 0.71f;
    static constexpr float SIGMA_G_RESCALE = 0.01f;
    static constexpr float SIGMA_M_RESCALE = 8.0f;

    CarmaAdapter(bool with_mag, const Vector3f& sigma_a_init, const Vector3f& sigma_g,
                 const Vector3f& sigma_m, const ocean_imu::kalman::carma::CarmaShape* hshape,
                 int tfg_reset_mode = -1) {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init * SIGMA_A_RESCALE;
        cfg_.sigma_g = sigma_g * SIGMA_G_RESCALE;
        cfg_.sigma_m = sigma_m * SIGMA_M_RESCALE;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();
        wave_direction::VesselRaoEqualizer::Config direction_rao;
        direction_rao.enabled = true;
        direction_rao.horizontal_x_tau_s = 1.0;
        direction_rao.horizontal_y_tau_s = 0.7;
        filter.setDirectionRao(direction_rao);
        wave_direction::VesselRaoNoiseWeighting low_wave_noise;
        low_wave_noise.response = direction_rao;
        low_wave_noise.max_std_scale = 4.0f;
        filter.setLowWaveNoiseWeighting(low_wave_noise);
        filter.enableTuner(true);
        filter.enableClamp(true);
        if (hshape) filter.setHorizontalShape(*hshape);
        if constexpr (requires { typename MekfT::ResetMode; }) {
            if (tfg_reset_mode >= 0)
                filter.mekf().set_reset_mode(static_cast<typename MekfT::ResetMode>(tfg_reset_mode));
        }
    }

    void updateMag(const Vector3f& mag_body_ned) override { fusion_.updateMag(mag_body_ned); }

    void update(float dt, const Vector3f& gyr, const Vector3f& acc, float temperature_c) override {
        fusion_.update(dt, gyr, acc, temperature_c);
    }

    // For a harness that configures the filter beyond what this adapter does
    // (carma_lever_arm-sim's filter-side lever arm).  A class-template member,
    // so the replays that never call it compile exactly as before.
    auto& fusion() noexcept { return fusion_; }

    FilterSnapshot snapshot() const override {
        const auto& filter = fusion_.raw();
        const auto& d = filter.dir();
        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(filter.mekf().get_position());
        s.vel_est_zu = ned_to_zu(filter.mekf().get_velocity());
        s.acc_est_zu = ned_to_zu(filter.mekf().get_world_accel());
        const Quaternionf q_bw_ned = filter.mekf().quaternion_boat().normalized();
        float roll_deg = 0.0f, pitch_deg = 0.0f, yaw_deg = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);
        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));
        s.acc_bias_est_ned = filter.mekf().get_acc_bias();
        s.gyro_bias_est_ned = filter.mekf().gyroscope_bias();
        s.mag_bias_est_ned_uT = fusion_.magHardIronBodyUT();
        s.tau_target = filter.getTauTarget();
        s.sigma_target = filter.getSigmaTarget();
        s.tuning_target = filter.getRSTarget();
        s.tau_applied = filter.getTauApplied();
        s.sigma_applied = filter.getSigmaApplied();
        s.tuning_applied = filter.getVerticalKappa();
        s.freq_hz = filter.getFreqHz();
        s.wave_period_sec = filter.getWavePeriodSec();
        s.period_sec = filter.getPeriodSec();
        s.accel_variance = filter.getAccelVariance();
        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps = filter.getVerticalSpeedEnvelopeMps(true);
        s.direction.phase = d.getPhase();
        s.direction.direction_deg = d.getAxisDegrees();
        s.direction.apparent_to_deg = filter.getApparentWaveDirectionToDeg();
        s.direction.apparent_from_deg = filter.getApparentWaveDirectionFromDeg();
        s.direction.sense_coherence = filter.getDirSenseCoherence();
        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getAxis());
        s.direction.uncertainty_deg = d.getAxisUncertaintyDegrees();
        s.direction.confidence = d.getLastStableConfidence();
        s.direction.amplitude = d.getAmplitude();
        s.direction.direction_vec = d.getAxis();
        s.direction.filtered_signal = d.getFilteredSignal();
        return s;
    }

private:
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, NHC, SH, NZC, SZ, MekfT>;
    mutable Fusion fusion_;
    typename Fusion::Config cfg_{};
};

#define TFG_ADAPTER(NAME, RESET)                                                        \
    class NAME final : public IW3dFusionAdapter {                                       \
    public:                                                                             \
        using Core = ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float, 4>;             \
        NAME(bool m, const Vector3f& a, const Vector3f& g, const Vector3f& s)           \
            : impl_(m, a, g, s, nullptr, static_cast<int>(Core::ResetMode::RESET)) {}   \
        void updateMag(const Vector3f& b) override { impl_.updateMag(b); }              \
        void update(float dt, const Vector3f& g, const Vector3f& a, float t) override { \
            impl_.update(dt, g, a, t);                                                  \
        }                                                                               \
        FilterSnapshot snapshot() const override { return impl_.snapshot(); }          \
    private:                                                                            \
        CarmaAdapter<4, false, 4, false, Core> impl_;                                   \
    };

TFG_ADAPTER(Adapter_T41, Full)
TFG_ADAPTER(Adapter_T41R0, None)

}  // namespace carma3_sim
