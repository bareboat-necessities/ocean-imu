/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Exercise the magnetic/no-GNSS sketch path with sample delays, IMU noise,
  residual biases, and a still -> waves -> still history. Truth generates
  the sensors only. Check the observer states as well as the displayed heave;
  the shared detrender cannot hide an unstable translational observer.
*/
#define EIGEN_NON_ARDUINO
#include "nlo/TimeVarGainNLO_Adapter.h"
#include "detrend/AdaptiveWaveDetrender.h"
#include "util/ImuLoopHelpers.h"
#include "AtomS3R/AtomS3R_ImuUnits.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

extern const float g_std = 9.80665f;
using Fusion = TimeVarGainNloAdapter<false, NloMagType::Magnetometer, float>;
using Vec3 = Eigen::Vector3f;
constexpr float gravity = atoms3r_ical::ImuCalCfg::g_cal_local;
constexpr float pi = 3.14159265358979323846f;

static bool check(bool ok, const char* message) {
    if (!ok) std::fprintf(stderr, "FAIL: %s\n", message);
    return ok;
}

// With gyro propagation, bias learning and position aiding absent, attitude
// feedback must preserve R*f_b + xi and must not create linear acceleration.
static bool forceCorrection() {
    using Core = TimeVaryingGainNLO<false, NloMagType::Magnetometer, double>;
    Core::Config cfg;
    cfg.use_time_varying_attitude_gains = false;
    cfg.use_time_varying_tmo_gain = false;
    cfg.use_triad_style_force_injection = false;
    cfg.kI_nominal = 0.0;
    cfg.K_xiz_p0z = cfg.K_vz_p0z = cfg.K_pz_p0z = 0.0;
    Core f(cfg);
    Core::Aux aux;
    aux.mag_valid = true;
    aux.mag_b = Core::Vec3(0.0, 20.0, 45.0);
    const Core::Vec3 force(2.0, 0.0, -cfg.gravity_mps2);
    const Core::Vec3 accel(2.0, 0.0, 0.0);
    double max_force_error = 0.0, max_velocity_error = 0.0;
    constexpr double dt = 0.005;
    for (int k = 0; k < 400; ++k) {
        f.update(dt, Core::Vec3::Zero(), force, aux);
        max_force_error = std::max(max_force_error, (f.specificForceNED() - force).norm());
        max_velocity_error = std::max(max_velocity_error,
            (f.velocityNED() - (k + 1)*dt*accel).norm());
    }
    std::printf("NLO force feedback force error=%.9g velocity error=%.9g\n",
                max_force_error, max_velocity_error);
    return check(max_force_error < 1e-10 && max_velocity_error < 1e-10,
                 "attitude feedback generated fictitious force/velocity");
}

struct Ramp { float value, d1, d2; };
static Ramp ramp(float t, float start, float duration) {
    const float u = std::clamp((t-start)/duration, 0.0f, 1.0f);
    return {u*u*u*(10.0f-15.0f*u+6.0f*u*u),
            30.0f*u*u*(1.0f-u)*(1.0f-u)/duration,
            60.0f*u*(1.0f-3.0f*u+2.0f*u*u)/(duration*duration)};
}

static bool replay(const char* name, bool delayed, bool jittered, bool waves,
                   float vertical_bias = 0.08f) {
    Fusion f;
    f.config().gravity_mps2 = gravity;
    f.config().filter.gravity_mps2 = gravity;
    f.reset();
    ocean_imu::ins::SampleDtTracker clock(0.005f, 0.05f);
    AdaptiveWaveDetrender detrender(defaultHeaveDetrenderConfig(0.30f));
    detrender.reset(0.0f);
    std::mt19937 rng(7);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    auto noise = [&](float sigma) {
        Vec3 v;
        for (int i = 0; i < 3; ++i) v(i) = sigma*normal(rng);
        return v;
    };
    float max_raw = 0.0f, max_report = 0.0f, max_display = 0.0f;
    float max_xi = 0.0f, max_tilt_error = 0.0f;
    float max_startup_heave = 0.0f;
    double wave_sq = 0.0, wave_time = 0.0, wave_sin = 0.0, wave_cos = 0.0;
    float freq = 0.30f;
    std::uint32_t sample_us = 0;
    int samples = 0;
    while (sample_us < 600000000u) {
        const float t = sample_us*1e-6f;
        const float dt = clock.update(sample_us);
        const auto on = ramp(t, 30.0f, 10.0f), off = ramp(t, 400.0f, 10.0f);
        const float e = waves ? on.value-off.value : 0.0f;
        const float ed = waves ? on.d1-off.d1 : 0.0f;
        const float edd = waves ? on.d2-off.d2 : 0.0f;
        constexpr float w = 2.0f*pi*0.2f;
        const float z = 0.3f*e*std::sin(w*t);
        const float az = 0.3f*((edd-e*w*w)*std::sin(w*t)+2.0f*ed*w*std::cos(w*t));
        const float pitch = 0.1f*e*std::sin(0.7f*t);
        const float pitch_rate = 0.1f*(ed*std::sin(0.7f*t)+0.7f*e*std::cos(0.7f*t));
        const Eigen::Quaternionf q(Eigen::AngleAxisf(pitch, Vec3::UnitY()));
        const Vec3 mag = q.conjugate()*Vec3(20.0f, 0.0f, 45.0f) + noise(0.4f);
        const Vec3 gyro = Vec3(0.005f, pitch_rate-0.005f, 0.005f) + noise(0.0016f);
        const Vec3 acc = q.conjugate()*Vec3(0.0f, 0.0f, az-gravity)
                         + Vec3(0.05f, -0.05f, vertical_bias) + noise(0.18f);
        f.setMagBody(mag, true);
        f.update(dt, gyro, acc);
        const auto s = f.snapshot();
        if (!check(s.q_nb.coeffs().allFinite() && s.position_ned.allFinite()
                   && s.velocity_ned.allFinite() && s.tvg.xi_n.allFinite(),
                   "device replay has non-finite observer state")) return false;
        if (s.tvg.wave_freq_hz > 1e-6f) freq = s.tvg.wave_freq_hz;
        const auto display = detrender.update(s.disp_zu.z(), dt, freq, freq > 1e-6f);
        max_raw = std::max(max_raw, std::abs(f.filter().positionNED().z()));
        max_report = std::max(max_report, std::abs(s.disp_zu.z()));
        if (t < 150.0f) max_startup_heave = std::max(max_startup_heave, std::abs(s.disp_zu.z()));
        max_display = std::max(max_display, std::abs(display.wave_clean));
        max_xi = std::max(max_xi, s.tvg.xi_norm);
        max_tilt_error = std::max(max_tilt_error,
            std::hypot(s.euler_rad.x(), s.euler_rad.y()-pitch)*180.0f/pi);
        if (waves && t > 250.0f && t < 390.0f) {
            const double err = s.position_ned.z()-z;
            wave_sq += dt*err*err;
            wave_time += dt;
            wave_sin += dt*s.position_ned.z()*std::sin(w*t);
            wave_cos += dt*s.position_ned.z()*std::cos(w*t);
        }
        ++samples;
        sample_us += delayed ? (jittered && samples%20 != 0 ? 5000u : 50000u) : 5000u;
    }
    const double wave_rms = std::sqrt(wave_sq/std::max(wave_time, 1.0));
    const double sin_amplitude = 2.0*wave_sin/std::max(wave_time, 1.0);
    const double cos_amplitude = 2.0*wave_cos/std::max(wave_time, 1.0);
    std::printf("NLO device %-15s raw=%.3f report=%.3f display=%.3f m "
                "xi=%.3f tilt_error=%.3f deg startup=%.3f m "
                "wave_rms=%.3f m wave_sin/cos=%.3f/%.3f m\n",
                name, max_raw, max_report, max_display, max_xi, max_tilt_error,
                max_startup_heave, wave_rms, sin_amplitude, cos_amplitude);
    bool ok = check(f.initialized(), "device replay never initialized");
    ok &= check(max_raw < 10.0f && max_report < 8.0f && max_display < 2.0f,
                "device heave drifted by metres beyond the noise/bias transient");
    ok &= check(max_startup_heave < 4.0f,
                "wave scheduling weakened capture before the residual was learned");
    ok &= check(max_xi < std::abs(vertical_bias)+0.22f && max_tilt_error < 5.0f,
                "delayed updates destabilized attitude/force feedback");
    // The noisy low-frequency residual is included in the printed raw RMS.
    // Check the response at the injected wave frequency independently, so
    // zero/frozen/clipped heave cannot pass a mere boundedness check.
    if (waves) ok &= check(std::abs(sin_amplitude-0.3) < 0.05
                          && std::abs(cos_amplitude) < 0.08,
                          "delayed updates lost wave amplitude/phase");
    return ok;
}

int main() {
    bool ok = forceCorrection();
    ok &= replay("200 Hz still", false, false, false);
    ok &= replay("20 Hz still", true, false, false);
    ok &= replay("jittered still", true, true, false);
    ok &= replay("jittered waves", true, true, true);
    ok &= replay("200 Hz 0.3 bias", false, false, false, 0.3f);
    return ok ? 0 : 1;
}
