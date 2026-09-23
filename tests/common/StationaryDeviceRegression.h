#pragma once
// Full, uninterrupted device-rate replays. Truth generates sensors only; it is
// never injected into the estimator. Check raw position, before any detrender.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#if defined(STATIONARY_DEVICE_TFG)
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#elif defined(STATIONARY_DEVICE_OU2)
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#else
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#endif
extern const float g_std = 9.80665f;
namespace stationary_device {
using V = Eigen::Vector3f;
#if defined(STATIONARY_DEVICE_TFG)
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
constexpr const char* name = "TFG";
constexpr bool deployed_fixed_cadence = false;
#else
constexpr bool deployed_fixed_cadence = true;
#if defined(STATIONARY_DEVICE_OU2)
using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
constexpr const char* name = "OU-II";
#else
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
constexpr const char* name = "OU-III";
#endif
#endif
constexpr float dt = 0.005f, pi = 3.14159265358979323846f;
auto& raw(Fusion& f) {
#if defined(STATIONARY_DEVICE_TFG)
    return f;
#else
    return f.raw();
#endif
}
Eigen::Quaternionf attitude(Fusion& f) {
#if defined(STATIONARY_DEVICE_TFG)
    return f.mekf().quaternion();
#else
    return f.raw().mekf().quaternion_boat();
#endif
}
float period(Fusion& f) {
#if defined(STATIONARY_DEVICE_TFG)
    return f.pseudoUpdatePeriodSec();
#else
    return f.raw().mekf().get_pseudo_update_period_s();
#endif
}
void begin(Fusion& f, bool fixed_cadence) {
    Fusion::Config cfg;
    cfg.sigma_a = V::Constant(0.12f);
    cfg.sigma_m = V::Constant(0.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
#if defined(STATIONARY_DEVICE_TFG)
    cfg.gyro_noise_density = 0.00135f;
#else
    cfg.sigma_g = V::Constant(0.00135f);
    cfg.enable_displacement_detrend = false;
#endif
    f.begin(cfg);
    raw(f).enableTuner(true);
    raw(f).setAccNoiseFloorSigma(0.12f);
    // fixed_cadence selects the fixed 15 ms cadence the OU-II/OU-III sketches
    // deploy.  The TFG sketch deploys the tau-scaled cadence, so for TFG the
    // dense/sparse comparison below is an ablation of the fixed option.
#if defined(STATIONARY_DEVICE_TFG)
    raw(f).setTauScaledPseudoCadence(!fixed_cadence);
#else
    raw(f).setTauScaledPseudoUpdateCadence(!fixed_cadence);
#endif
}
struct Metrics {
    bool finite = true;
    float live_time = -1.0f;
    float max_yaw_deg = 0.0f, max_step_m = 0.0f, max_position_m = 0.0f;
    float final_period_s = 0.0f, quiet_max_step_m = 0.0f;
    float handoff_yaw_deg = 0.0f, tail_yaw_deg = 0.0f;
    double wave_sq = 0.0, quiet_sq = 0.0;
    int wave_n = 0, quiet_n = 0;
    double wave_rms() const { return std::sqrt(wave_sq / std::max(1, wave_n)); }
    double quiet_rms() const { return std::sqrt(quiet_sq / std::max(1, quiet_n)); }
};
struct Envelope { float value, d1, d2; };
Envelope ramp(float t, float start, float duration) {
    const float u = std::clamp((t - start) / duration, 0.0f, 1.0f);
    const float u2 = u*u, u3 = u2*u;
    return {u3*(10.0f-15.0f*u+6.0f*u2),
            30.0f*u2*(1.0f-u)*(1.0f-u)/duration,
            60.0f*u*(1.0f-3.0f*u+2.0f*u2)/(duration*duration)};
}
Metrics replay(bool dense, float bias, float gyro_bias, float noise_std, bool waves) {
    Fusion f;
    begin(f, dense);
    Metrics m;
    std::mt19937 rng(1234);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    float previous_p = 0.0f;
    const int count = waves ? 150000 : 60000;
    for (int k = 0; k < count; ++k) {
        const float t = static_cast<float>(k) * dt;
        const auto on = ramp(t, 150.0f, 20.0f), off = ramp(t, 450.0f, 20.0f);
        const float e = waves ? on.value - off.value : 0.0f;
        const float ed = waves ? on.d1 - off.d1 : 0.0f;
        const float edd = waves ? on.d2 - off.d2 : 0.0f;
        constexpr float wr = 2.0f*pi*0.17f, wp = 2.0f*pi*0.13f, w = 2.0f*pi*0.2f;
        const float r = 0.15f + 0.08f*e*std::sin(wr*t);
        const float p = 0.1f + 0.05f*e*std::sin(wp*t);
        const float rd = 0.08f*(ed*std::sin(wr*t) + e*wr*std::cos(wr*t));
        const float pd = 0.05f*(ed*std::sin(wp*t) + e*wp*std::cos(wp*t));
        const Eigen::Quaternionf q = Eigen::AngleAxisf(0.5f,V::UnitZ()) *
            Eigen::AngleAxisf(p,V::UnitY()) * Eigen::AngleAxisf(r,V::UnitX());
        const float z = 0.3f*e*std::sin(w*t);
        const float az = 0.3f*((edd - e*w*w)*std::sin(w*t) + 2.0f*ed*w*std::cos(w*t));
        V a = q.conjugate()*V(0,0,az-g_std) + V(0,0,bias);
        V gyro(rd, pd*std::cos(r), -pd*std::sin(r)+gyro_bias);
        for (int j=0; j<3; ++j) {
            a[j] += noise_std*normal(rng);
            gyro[j] += (noise_std > 0.0f ? 0.00135f : 0.0f)*normal(rng);
        }
        f.update(dt, gyro, a, 35.0f);
        if (k % 8 == 0) {
            V mag = q.conjugate()*V(20,0,43);
            for (int j=0; j<3; ++j) mag[j] += (noise_std > 0.0f ? 0.1f : 0.0f)*normal(rng);
            f.updateMag(mag);
        }
        if (!f.isLive()) continue;
        const bool first_live = m.live_time < 0.0f;
        if (first_live) m.live_time = t;
        const float position = raw(f).mekf().get_position().z();
        const float step = std::abs(position - previous_p);
        previous_p = position;
        const Eigen::Matrix3f R = attitude(f).toRotationMatrix();
        const float yaw_error = std::abs(std::remainder(std::atan2(R(1,0),R(0,0))-0.5f,2.0f*pi))*180.0f/pi;
        m.finite = m.finite && std::isfinite(position) && std::isfinite(yaw_error) &&
            raw(f).mekf().covariance_full().allFinite();
        m.max_yaw_deg = std::max(m.max_yaw_deg,yaw_error);
        if (first_live) m.handoff_yaw_deg = yaw_error;
        if (!waves && t >= 240.0f) m.tail_yaw_deg = std::max(m.tail_yaw_deg,yaw_error);
        m.max_step_m = std::max(m.max_step_m,step);
        m.max_position_m = std::max(m.max_position_m,std::abs(position));
        if (waves && t >= 270.0f && t < 430.0f) {
            const double error = static_cast<double>(position-z);
            m.wave_sq += error*error; ++m.wave_n;
        }
        if (waves && t >= 650.0f) {
            m.quiet_sq += static_cast<double>(position)*position; ++m.quiet_n;
            m.quiet_max_step_m = std::max(m.quiet_max_step_m,step);
        }
        if (!m.finite) break;
    }
    m.final_period_s = period(f);
    return m;
}
void print(const char* label, const Metrics& m) {
    std::cout << name << ' ' << label << " live_s=" << m.live_time
              << " max_yaw_deg=" << m.max_yaw_deg << " max_step_m=" << m.max_step_m
              << " max_position_m=" << m.max_position_m << " period_s=" << m.final_period_s
              << " handoff_yaw_deg=" << m.handoff_yaw_deg << " tail_yaw_deg=" << m.tail_yaw_deg
              << " wave_rms_m=" << m.wave_rms() << " quiet_rms_m=" << m.quiet_rms()
              << " quiet_max_step_m=" << m.quiet_max_step_m << '\n';
}
int run() {
    int failures = 0;
    auto check = [&](bool ok, const char* text) {
        if (!ok) { std::cerr << "FAIL: " << name << ' ' << text << '\n'; ++failures; }
    };
    const auto north = replay(deployed_fixed_cadence,0.03f,0.025f,0.0f,false);
    print("stationary-north",north);
    check(north.finite && north.live_time >= 0.0f,"stationary north replay did not stay finite/enter Live");
    check(north.handoff_yaw_deg < 3.0f,"stationary handoff lost magnetic north");
    check(north.tail_yaw_deg < 0.5f,"stationary magnetic heading did not settle");
    // OU-II/III retain their existing short gyro-bias-learning transient. The
    // TFG regression additionally rejects drift during acquisition/refinement;
    // the old gyro-only gauge produces a 90-degree error in this same replay.
#if defined(STATIONARY_DEVICE_TFG)
    check(north.max_yaw_deg < 3.0f,"stationary gyro bias corrupted magnetic heading at handoff/refinement");
#else
    check(north.max_yaw_deg < 6.0f,"OU stationary heading transient grew unexpectedly");
#endif
    const auto sparse = replay(false,0.2f,0.01f,0.12f,false);
    const auto dense = replay(true,0.2f,0.01f,0.12f,false);
    print("stationary-tau-scaled",sparse); print("stationary-fixed-cadence",dense);
    check(sparse.finite && dense.finite && dense.live_time >= 0.0f,"stationary heave replay failed");
    check(std::abs(dense.final_period_s-0.015f)<1e-6f,"fixed pseudo cadence is not 15 ms");
    check(dense.max_step_m < 0.05f,"device still has a large stationary raw-position tooth");
    check(dense.max_step_m <= 0.25f*sparse.max_step_m + 0.001f,"frequent virtual constraints did not reduce stationary teeth");
    const auto wave_sparse = replay(false,0.03f,0.001f,0.0148f,true);
    const auto wave_dense = replay(true,0.03f,0.001f,0.0148f,true);
    print("rest-wave-rest-tau-scaled",wave_sparse); print("rest-wave-rest-fixed-cadence",wave_dense);
    check(wave_sparse.finite && wave_dense.finite && wave_dense.wave_n > 0 && wave_dense.quiet_n > 0,
          "uninterrupted rest-wave-rest replay failed");
    check(wave_dense.wave_rms() <= 1.10*wave_sparse.wave_rms()+0.001,
          "device cadence materially degraded wave tracking");
    check(wave_dense.quiet_max_step_m < 0.01f,"return to rest retained large raw-position teeth");
    return failures ? 1 : 0;
}
} // namespace stationary_device
