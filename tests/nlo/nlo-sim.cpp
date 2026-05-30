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
#include "nlo/TimeVarGainNLO_Adapter.h"

using Eigen::Vector3f;

bool add_noise = true;

static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 22.0f,
    .err_limit_percent_z_pmstokes  = 22.0f,
    .err_limit_yaw_deg             = 8.0f,
    .err_limit_percent_3d_jonswap  = 9999.0f,
    .err_limit_percent_3d_pmstokes = 9999.0f,
    .acc_z_bias_percent            = 9999.0f,
    .bias_3d_percent               = 9999.0f,
};

class FusionAdapterTimeVarGainNLO_NoGnssNoMag final
    : public IW3dFusionAdapterTyped<TvgNloFilterSnapshot> {
public:
    using Snapshot = TvgNloFilterSnapshot;
    using Core = TimeVarGainNloAdapter<false, NloMagType::None, float>;

    FusionAdapterTimeVarGainNLO_NoGnssNoMag(const Vector3f& sigma_a_init,
                                            const Vector3f& sigma_g)
        : core_(make_config_(sigma_a_init, sigma_g))
    {
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        (void)mag_body_ned;
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        (void)temperature_c;
        core_.update(dt, gyr_meas_ned, acc_meas_ned);
    }

    Snapshot snapshot() const override {
        const auto cs = core_.snapshot();

        Snapshot s;

        s.disp_est_zu = cs.disp_zu.template cast<float>();
        s.vel_est_zu  = cs.vel_zu.template cast<float>();
        s.acc_est_zu  = cs.acc_zu.template cast<float>();

        /*
          Existing harness comparison convention for this NLO path.
          Mag=None: yaw is unobservable, so force yaw to zero.
        */
        s.euler_nautical_deg = Vector3f(
            rad_to_deg(float(cs.euler_rad.y())),
            rad_to_deg(float(cs.euler_rad.x())),
            0.0f
        );

        s.acc_bias_est_ned = Vector3f::Zero();
        s.gyro_bias_est_ned = cs.tvg.gyro_bias_b.template cast<float>();
        s.mag_bias_est_ned_uT = Vector3f::Zero();

        s.tvg.k1 = float(cs.tvg.k1);
        s.tvg.k2 = float(cs.tvg.k2);
        s.tvg.kI = float(cs.tvg.kI);
        s.tvg.vartheta = float(cs.tvg.vartheta);
        s.tvg.p0z_hat = float(cs.tvg.p0z_hat);

        s.tvg.xi_n = cs.tvg.xi_n.template cast<float>();
        s.tvg.fhat_n = cs.tvg.fhat_n.template cast<float>();
        s.tvg.sigma_b = cs.tvg.sigma_b.template cast<float>();
        s.tvg.gyro_bias_b = cs.tvg.gyro_bias_b.template cast<float>();

        s.tvg.xi_norm = float(cs.tvg.xi_norm);
        s.tvg.fhat_norm = float(cs.tvg.fhat_norm);
        s.tvg.sigma_norm = float(cs.tvg.sigma_norm);
        s.tvg.gyro_bias_norm = float(cs.tvg.gyro_bias_norm);

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
    static Core::Config make_config_(const Vector3f& sigma_a_init,
                                     const Vector3f& sigma_g)
    {
        (void)sigma_a_init;
        (void)sigma_g;

        auto cfg = Core::makeDefaultConfig();

        cfg.gravity_mps2 = g_std;
        cfg.filter.gravity_mps2 = g_std;

        cfg.init_required_good_time_s = 2.0f;
        cfg.init_max_wait_s = 4.0f;
        cfg.init_gyro_max_rad_s = 0.08f;
        cfg.init_acc_norm_tol_frac = 0.18f;
        cfg.yaw_seed_rad = 0.0f;

        /*
          Critical for this no-GNSS/no-mag wave sim:
          do not run accelerometer tilt trim during waves.
          It chases horizontal wave acceleration and worsens roll/pitch.
        */
        cfg.tilt_trim_enabled = false;
        cfg.tilt_trim_duration_s = 0.0f;
        cfg.tilt_trim_tau_s = 12.0f;
        cfg.tilt_trim_acc_lpf_tau_s = 0.75f;
        cfg.tilt_trim_max_rate_rad_s = 0.025f;
        cfg.tilt_trim_gyro_max_rad_s = 0.08f;
        cfg.tilt_trim_acc_norm_tol_frac = 0.10f;

        cfg.tilt_bias_enabled = false;
        cfg.tilt_bias_ki = 0.0f;
        cfg.tilt_bias_limit_rad_s = 0.0f;

        cfg.run_filter_before_initialized = false;

        cfg.filter.gyro_bias_limit_rad_s = 0.10f;
        cfg.filter.max_specific_force_mps2 = 30.0f;

        cfg.filter.use_time_varying_attitude_gains = true;
        cfg.filter.attitude_gain_tau_s = 20.0f;
        cfg.filter.attitude_gain_switch_s = 40.0f;

        cfg.filter.k1_initial = 4.0f;
        cfg.filter.k2_initial = 0.0f;
        cfg.filter.kI_initial = 0.03f;

        cfg.filter.k1_nominal = 0.70f;
        cfg.filter.k2_nominal = 0.0f;
        cfg.filter.kI_nominal = 0.004f;

        cfg.filter.K_p0z_p0z = 5.4295f;
        cfg.filter.K_pz_p0z  = 2.2396f;
        cfg.filter.K_vz_p0z  = 0.4454f;
        cfg.filter.K_xiz_p0z = 0.0354f;

        cfg.filter.K_pp_scalar  = 0.0f;
        cfg.filter.K_vp_scalar  = 0.0f;
        cfg.filter.K_xip_scalar = 0.0f;

        cfg.filter.theta = 1.0f;

        cfg.filter.use_time_varying_tmo_gain = true;
        cfg.filter.vartheta0 = 0.65f;
        cfg.filter.vartheta1_without_gnss = 0.0f;
        cfg.filter.vartheta1_a = 2.0f;
        cfg.filter.vartheta1_b = 0.0f;
        cfg.filter.gnss_rms_lpf_tau_s = 125.0f;

        cfg.filter.vartheta2_tau_s = 20.0f;
        cfg.filter.vartheta2_switch_s = 30.0f;

        cfg.filter.p0z_highpass_tau_s = 600.0f;
        cfg.filter.use_triad_style_force_injection = true;

        return cfg;
    }

private:
    Core core_;
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
    float mean_z_err = 0.0f;
    size_t n_z = 0;

    for (size_t i = start; i < result.errs_z.size(); ++i) {
        rms_z.add(result.errs_z[i]);
        rms_roll.add(result.errs_roll[i]);
        rms_pitch.add(result.errs_pitch[i]);
        rms_yaw.add(result.errs_yaw[i]);

        if (std::isfinite(result.errs_z[i])) {
            mean_z_err += result.errs_z[i];
            ++n_z;
        }
    }

    if (n_z > 0) {
        mean_z_err /= static_cast<float>(n_z);
    }

    RMSReport rms_z_demeaned;
    for (size_t i = start; i < result.errs_z.size(); ++i) {
        if (std::isfinite(result.errs_z[i])) {
            rms_z_demeaned.add(result.errs_z[i] - mean_z_err);
        }
    }

    const float z_rms = rms_z.rms();
    const float z_pct = 100.0f * z_rms / result.wave_params.height;

    const float z_rms_demeaned = rms_z_demeaned.rms();
    const float z_pct_demeaned = 100.0f * z_rms_demeaned / result.wave_params.height;

    const auto& d = result.final_snapshot.tvg;

    std::cout << "=== Last 60 s TimeVarGain NLO VVR-only summary for "
              << result.output_name << " ===\n";

    std::cout << "Z RMS raw (m): " << z_rms << "\n";
    std::cout << "Z RMS raw (%Hs): " << z_pct
              << "% (Hs=" << result.wave_params.height << ")\n";
    std::cout << "Z mean error (m): " << mean_z_err << "\n";
    std::cout << "Z RMS de-meaned (m): " << z_rms_demeaned << "\n";
    std::cout << "Z RMS de-meaned (%Hs): " << z_pct_demeaned << "%\n";

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

    float mean_z_err = 0.0f;
    size_t n_z = 0;

    for (size_t i = start; i < result.errs_z.size(); ++i) {
        if (std::isfinite(result.errs_z[i])) {
            mean_z_err += result.errs_z[i];
            ++n_z;
        }
    }

    if (n_z == 0) return;

    mean_z_err /= static_cast<float>(n_z);

    RMSReport rms_z_demeaned;
    for (size_t i = start; i < result.errs_z.size(); ++i) {
        if (std::isfinite(result.errs_z[i])) {
            rms_z_demeaned.add(result.errs_z[i] - mean_z_err);
        }
    }

    const float z_pct =
        100.0f * rms_z_demeaned.rms() / result.wave_params.height;

    const float z_limit = (result.wave_type == WaveType::JONSWAP)
        ? FAIL_LIMITS.err_limit_percent_z_jonswap
        : FAIL_LIMITS.err_limit_percent_z_pmstokes;

    if (z_pct > z_limit) {
        std::cerr << "ERROR: de-meaned Z RMS above limit (" << z_pct << "% > "
                  << z_limit << "%). Mean Z error was " << mean_z_err
                  << " m. Failing.\n";
        std::exit(EXIT_FAILURE);
    }
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
                         "as TimeVarGainNloAdapter<false, NloMagType::None>.\n";
        } else if (arg == "--with-gnss") {
            std::cerr << "WARNING: --with-gnss ignored. This simulator is compiled "
                         "as TimeVarGainNloAdapter<false, NloMagType::None>.\n";
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
