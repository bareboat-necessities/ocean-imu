#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <optional>
#include <cstdlib>
#include <cmath>
#include <numbers>
#include <algorithm>

/*
    Copyright (c), 2026  Mikhail Grushinskiy
*/

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "nlo/TimeVaryingGainNLO.h"   // adjust path if needed

using Eigen::Vector3f;

bool add_noise = true;

static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 22.0f,
    .err_limit_percent_z_pmstokes  = 25.0f,
    .err_limit_yaw_deg             = 15.0f,
    .err_limit_percent_3d_jonswap  = 9999.0f,
    .err_limit_percent_3d_pmstokes = 9999.0f,
    .acc_z_bias_percent            = 9999.0f,
    .bias_3d_percent               = 9999.0f,
};

class FusionAdapterTimeVarGainNLO_NoGnssNoMag final
    : public IW3dFusionAdapterTyped<TvgNloFilterSnapshot> {
public:
    using Snapshot = TvgNloFilterSnapshot;
    using Filter = TimeVaryingGainNLO<false, NloMagType::None, float>;

    FusionAdapterTimeVarGainNLO_NoGnssNoMag(const Vector3f& sigma_a_init,
                                            const Vector3f& sigma_g)
        : filter_(make_config_(sigma_a_init, sigma_g))
    {
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        (void)mag_body_ned;
        // Mag=None: no-op.
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        (void)temperature_c;

        /*
          Runner-facing convention:
            gyr_meas_ned: body-frame NED axes, rad/s
            acc_meas_ned: body-frame NED specific force, m/s^2

          TimeVaryingGainNLO<false, None> expects the same body-NED convention.
          At rest, level:
            acc_meas_ned ≈ [0, 0, -g]
        */

        if (!initialized_) {
            init_acc_sum_ += acc_meas_ned;
            init_time_s_ += dt;
            ++init_count_;

            if (init_time_s_ >= INIT_SECONDS && init_count_ > 0) {
                const Vector3f acc0 =
                    init_acc_sum_ / static_cast<float>(init_count_);

                // No yaw reference. Seed yaw = 0.
                initialized_ = filter_.initializeFromAccel(acc0, 0.0f);

                if (!initialized_) {
                    init_acc_sum_.setZero();
                    init_time_s_ = 0.0f;
                    init_count_ = 0;
                }
            }

            return;
        }

        Filter::Aux aux{};
        filter_.update(dt, gyr_meas_ned, acc_meas_ned, aux);
    }

    Snapshot snapshot() const override {
        Snapshot s;

        const Vector3f p_n = filter_.positionNED();
        const Vector3f v_n = filter_.velocityNED();
        const Vector3f fhat_n = filter_.specificForceNED();

        /*
          Filter state is NED, +Z down.
          W3D harness output is Z-up for displacement/velocity/acceleration.
        */
        const float a_z_down = fhat_n.z() + g_std;
        const float a_z_up = -a_z_down;

        s.disp_est_zu = Vector3f(0.0f, 0.0f, -p_n.z());
        s.vel_est_zu  = Vector3f(0.0f, 0.0f, -v_n.z());
        s.acc_est_zu  = Vector3f(0.0f, 0.0f, a_z_up);

        /*
          Mag=None: absolute yaw is unobservable.
          Report roll/pitch, force yaw to 0 so the harness does not treat
          free yaw drift as meaningful.
        */
        const Vector3f euler_rad = filter_.eulerRad();
        s.euler_nautical_deg = Vector3f(
            rad_to_deg(euler_rad.y()),
            rad_to_deg(euler_rad.x()),
            0.0f
        );

        s.acc_bias_est_ned = Vector3f::Zero();

        s.gyro_bias_est_ned = filter_.gyroBiasBody();
        s.mag_bias_est_ned_uT = Vector3f::Zero();

        s.tvg.k1 = filter_.gainK1();
        s.tvg.k2 = filter_.gainK2();
        s.tvg.kI = filter_.gainKI();
        s.tvg.vartheta = filter_.gainVartheta();
        s.tvg.p0z_hat = filter_.integratedVerticalPositionState();

        s.tvg.xi_n = filter_.xiNED();
        s.tvg.fhat_n = fhat_n;
        s.tvg.sigma_b = filter_.sigmaBody();
        s.tvg.gyro_bias_b = filter_.gyroBiasBody();

        s.tvg.xi_norm = s.tvg.xi_n.norm();
        s.tvg.fhat_norm = s.tvg.fhat_n.norm();
        s.tvg.sigma_norm = s.tvg.sigma_b.norm();
        s.tvg.gyro_bias_norm = s.tvg.gyro_bias_b.norm();

        s.direction.phase = NAN;
        s.direction.direction_deg = NAN;
        s.direction.direction_deg_generator_signed = NAN;
        s.direction.uncertainty_deg = NAN;
        s.direction.confidence = NAN;
        s.direction.amplitude = NAN;
        s.direction.direction_vec = Eigen::Vector2f::Zero();
        s.direction.filtered_signal = Eigen::Vector2f::Zero();
        s.direction.sign = UNCERTAIN;
        s.direction.sign_num = 0;

        return s;
    }

private:
    static Filter::Config make_config_(const Vector3f& sigma_a_init,
                                       const Vector3f& sigma_g)
    {
        (void)sigma_a_init;
        (void)sigma_g;

        Filter::Config cfg{};

        cfg.gravity_mps2 = g_std;

        /*
          No GNSS and no yaw reference:
            - vertical VVR channel remains useful
            - x/y position are not corrected
            - yaw and yaw gyro bias are unobservable

          Keep gains much more conservative than the paper's full
          GNSS+compass startup observer.
        */
        cfg.gyro_bias_limit_rad_s = 0.10f;
        cfg.max_specific_force_mps2 = 30.0f;

        cfg.use_time_varying_attitude_gains = true;
        cfg.attitude_gain_tau_s = 20.0f;
        cfg.attitude_gain_switch_s = 40.0f;

        // k2 is unused when Mag=None.
        cfg.k1_initial = 4.0f;
        cfg.k2_initial = 0.0f;
        cfg.kI_initial = 0.03f;

        cfg.k1_nominal = 0.70f;
        cfg.k2_nominal = 0.0f;
        cfg.kI_nominal = 0.004f;

        // Vertical VVR gains from the paper-like Riccati example.
        cfg.K_p0z_p0z = 5.4295f;
        cfg.K_pz_p0z  = 2.2396f;
        cfg.K_vz_p0z  = 0.4454f;
        cfg.K_xiz_p0z = 0.0354f;

        // GNSS horizontal correction is compiled out by WithGNSS=false.
        cfg.K_pp_scalar  = 0.0f;
        cfg.K_vp_scalar  = 0.0f;
        cfg.K_xip_scalar = 0.0f;

        cfg.theta = 1.0f;

        /*
          With no GNSS RMS input, do not let vartheta depend on PosRef quality.
          Use only fixed base + short startup transient.
        */
        cfg.use_time_varying_tmo_gain = true;
        cfg.vartheta0 = 0.65f;
        cfg.vartheta1_without_gnss = 0.0f;
        cfg.vartheta1_a = 2.0f;
        cfg.vartheta1_b = 0.0f;
        cfg.gnss_rms_lpf_tau_s = 125.0f;

        cfg.vartheta2_tau_s = 20.0f;
        cfg.vartheta2_switch_s = 30.0f;

        // Keep the paper's slow HP on integrated vertical innovation.
        cfg.p0z_highpass_tau_s = 600.0f;

        cfg.use_triad_style_force_injection = true;

        return cfg;
    }

private:
    static constexpr float INIT_SECONDS = 0.25f;

    bool initialized_ = false;
    float init_time_s_ = 0.0f;
    int init_count_ = 0;
    Vector3f init_acc_sum_ = Vector3f::Zero();

    Filter filter_;
};

static std::vector<float> finite_tail_values(const std::vector<float>& v,
                                             size_t start)
{
    std::vector<float> out;
    if (v.size() <= start) return out;

    out.reserve(v.size() - start);
    for (size_t i = start; i < v.size(); ++i) {
        if (std::isfinite(v[i])) {
            out.push_back(v[i]);
        }
    }

    return out;
}

static void print_tvg_nlo_vertical_summary(const TvgNloSimulationRunResult& result,
                                           float dt)
{
    constexpr float RMS_WINDOW_SEC = 60.0f;
    const int N_last = static_cast<int>(RMS_WINDOW_SEC / dt);
    if (result.errs_z.size() <= static_cast<size_t>(N_last)) return;

    const size_t start = result.errs_z.size() - N_last;

    RMSReport rms_z, rms_roll, rms_pitch, rms_yaw;
    for (size_t i = start; i < result.errs_z.size(); ++i) {
        rms_z.add(result.errs_z[i]);
        rms_roll.add(result.errs_roll[i]);
        rms_pitch.add(result.errs_pitch[i]);
        rms_yaw.add(result.errs_yaw[i]);
    }

    const float z_rms = rms_z.rms();
    const float z_pct = 100.0f * z_rms / result.wave_params.height;

    const auto& d = result.final_snapshot.tvg;

    std::cout << "=== Last 60 s TimeVarGain NLO VVR-only summary for "
              << result.output_name << " ===\n";

    std::cout << "Z RMS (m): " << z_rms << "\n";
    std::cout << "Z RMS (%Hs): " << z_pct
              << "% (Hs=" << result.wave_params.height << ")\n";

    std::cout << "Angles RMS (deg): Roll=" << rms_roll.rms()
              << " Pitch=" << rms_pitch.rms()
              << " Yaw(free/ignored)=" << rms_yaw.rms() << "\n";

    std::cout << "tvg_k1=" << d.k1
              << ", tvg_k2=" << d.k2
              << ", tvg_kI=" << d.kI
              << ", tvg_vartheta=" << d.vartheta << "\n";

    std::cout << "tvg_p0z_hat=" << d.p0z_hat << "\n";

    std::cout << "tvg_xi_n=("
              << d.xi_n.x() << ", "
              << d.xi_n.y() << ", "
              << d.xi_n.z() << "), norm=" << d.xi_norm << "\n";

    std::cout << "tvg_fhat_n=("
              << d.fhat_n.x() << ", "
              << d.fhat_n.y() << ", "
              << d.fhat_n.z() << "), norm=" << d.fhat_norm << "\n";

    std::cout << "tvg_sigma_b=("
              << d.sigma_b.x() << ", "
              << d.sigma_b.y() << ", "
              << d.sigma_b.z() << "), norm=" << d.sigma_norm << "\n";

    std::cout << "tvg_gyro_bias_b=("
              << d.gyro_bias_b.x() << ", "
              << d.gyro_bias_b.y() << ", "
              << d.gyro_bias_b.z() << "), norm=" << d.gyro_bias_norm << "\n";

    const std::vector<float> finite_freq =
        finite_tail_values(result.freq_hist, start);

    if (!finite_freq.empty()) {
        std::cout << "freq_hz: mean=" << mean_vec(finite_freq)
                  << " median=" << median_vec(finite_freq)
                  << " p05=" << percentile_vec(finite_freq, 0.05)
                  << " p95=" << percentile_vec(finite_freq, 0.95) << "\n";
    } else {
        std::cout << "freq_hz: n/a for this observer\n";
    }

    std::cout << "===============================================================\n\n";
}

static void fail_if_tvg_nlo_vertical_gates_breached(const TvgNloSimulationRunResult& result,
                                                    float dt)
{
    constexpr float RMS_WINDOW_SEC = 60.0f;
    const int N_last = static_cast<int>(RMS_WINDOW_SEC / dt);
    if (result.errs_z.size() <= static_cast<size_t>(N_last)) return;

    const size_t start = result.errs_z.size() - N_last;

    RMSReport rms_z;
    for (size_t i = start; i < result.errs_z.size(); ++i) {
        rms_z.add(result.errs_z[i]);
    }

    const float z_pct = 100.0f * rms_z.rms() / result.wave_params.height;
    const float z_limit = (result.wave_type == WaveType::JONSWAP)
        ? FAIL_LIMITS.err_limit_percent_z_jonswap
        : FAIL_LIMITS.err_limit_percent_z_pmstokes;

    if (z_pct > z_limit) {
        std::cerr << "ERROR: Z RMS above limit (" << z_pct << "% > "
                  << z_limit << "%). Failing.\n";
        std::exit(EXIT_FAILURE);
    }

    // Mag=None: yaw is unobservable, so no yaw gate.
}

static std::optional<TvgNloSimulationRunResult>
process_wave_file_for_tvg_nlo_nomag_nognss(const std::string& filename,
                                           float dt,
                                           bool add_noise)
{
    const float acc_sigma = 1.51e-3f * g_std;
    const float gyr_sigma = 0.00157f;

    const float acc_bias_range = 5e-3f * g_std;
    const float gyr_bias_range =
        0.05f * float(std::numbers::pi_v<float> / 180.0f);

    const float acc_bias_rw = 0.0005f;
    const float gyr_bias_rw = 0.00001f;

    SimulationNoiseModels noise_models;

    noise_models.accel_noise = make_imu_noise_model(
        acc_sigma,
        acc_bias_range,
        acc_bias_rw,
        1234
    );

    noise_models.gyro_noise = make_imu_noise_model(
        gyr_sigma,
        gyr_bias_range,
        gyr_bias_rw,
        5678
    );

    /*
      No magnetometer is used in this run. Leave mag_noise unset.
      TvgNloSimulationRunner also has options.with_mag=false below.
    */

    const Vector3f sigma_a_init(
        2.8f * acc_sigma,
        2.8f * acc_sigma,
        2.8f * acc_sigma
    );

    const Vector3f sigma_g(
        2.0f * gyr_sigma,
        2.0f * gyr_sigma,
        2.0f * gyr_sigma
    );

    FusionAdapterTimeVarGainNLO_NoGnssNoMag adapter(
        sigma_a_init,
        sigma_g
    );

    W3dSimulationOptions options;
    options.dt = dt;
    options.with_mag = false;
    options.add_noise = add_noise;
    options.mag_odr_hz = 0.0f;
    options.temperature_c = 35.0f;
    options.output_suffix_with_mag = "_tvg_nlo";
    options.output_suffix_no_mag = "_tvg_nlo_nomag_nognss";

    TvgNloSimulationRunner runner(options, std::move(noise_models), adapter);
    return runner.run(filename);
}

int main(int argc, char* argv[])
{
    const float dt = 1.0f / 200.0f;
    add_noise = true;

    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];

        if (arg == "--no-noise") {
            add_noise = false;
        } else if (arg == "--nomag") {
            // Accepted for compatibility. This simulator is always Mag=None.
        } else if (arg == "--with-mag") {
            std::cerr << "WARNING: --with-mag ignored. This simulator is compiled "
                         "as TimeVaryingGainNLO<false, NloMagType::None>.\n";
        } else if (arg == "--with-gnss") {
            std::cerr << "WARNING: --with-gnss ignored. This simulator is compiled "
                         "as TimeVaryingGainNLO<false, NloMagType::None>.\n";
        }
    }

    std::cout << "TimeVarGain NLO VVR-only simulation starting"
              << " with_gnss=false"
              << ", mag=None"
              << ", noise=" << (add_noise ? "true" : "false")
              << "\n";

    const auto files = collect_wave_data_files(".");

    for (const auto& fname : files) {
        auto result = process_wave_file_for_tvg_nlo_nomag_nognss(
            fname,
            dt,
            add_noise
        );

        if (!result) continue;

        print_tvg_nlo_vertical_summary(*result, dt);
        fail_if_tvg_nlo_vertical_gates_breached(*result, dt);
    }

    return 0;
}
