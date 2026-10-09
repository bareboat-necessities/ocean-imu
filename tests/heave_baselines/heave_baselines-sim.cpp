#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

/*
  Copyright 2026, Mikhail Grushinskiy

  Literature heave baselines on the OU-II/OU-III replay.

  Every method sees the same corrupted IMU stream as the OU-III simulator:
  the record, the noise template (process_wave_file_for_tracker) and the
  seeds are shared, and the scored window is the same trailing 900 s.
  The three baselines share one Mahony attitude front end, the tuned,
  sea-state-adaptive one of the PII observer, and differ only in how they turn
  its levelled vertical acceleration into heave:

    pii      AdaptiveVerticalPII (this repository's non-Kalman observer)
    godhavn  Godhavn 1998 adaptive fourth-order heave filter
    kuchler  Kuchler et al. 2011 FFT-identified harmonic-mode EKF
*/

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include "heave_baselines/GodhavnHeaveFilter.h"
#include "heave_baselines/KuchlerHeaveEKF.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

constexpr float RMS_WINDOW_SEC = 900.0f;

enum class Method { PII, Godhavn, Kuchler };

const char* method_name(Method m) {
    switch (m) {
        case Method::PII: return "pii";
        case Method::Godhavn: return "godhavn";
        case Method::Kuchler: return "kuchler";
    }
    return "?";
}

bool env_double(const char* name, double& out) {
    if (const char* s = std::getenv(name)) {
        out = std::strtod(s, nullptr);
        return true;
    }
    return false;
}

// Regression sentinels for the deterministic single-realization protocol, not
// targets: the worst Z RMS (%Hs) of each method across the scored records plus
// about half a percent.
struct HeaveLimits { float jonswap; float pmstokes; };
constexpr HeaveLimits PII_LIMITS{7.74f, 8.99f};       // worst 7.70 / 8.94 (H8.5)
constexpr HeaveLimits GODHAVN_LIMITS{14.80f, 18.03f}; // worst 14.72 / 17.94 (H8.5)
constexpr HeaveLimits KUCHLER_LIMITS{16.25f, 18.52f}; // worst 16.17 / 18.42 (H8.5)

template <Method M>
class HeaveBaselineAdapter final : public IW3dFusionAdapter {
public:
    using Frontend = marine_obs::AdaptiveVerticalPIIMahony<float, true, TrackerType::PLL>;

    HeaveBaselineAdapter(bool with_mag,
                         const Vector3f& sigma_a_init,
                         const Vector3f& sigma_g,
                         const Vector3f& sigma_m)
        : with_mag_(with_mag), frontend_(make_frontend_config_(with_mag))
    {
        (void)sigma_a_init;
        (void)sigma_g;
        (void)sigma_m;

        heave_baselines::GodhavnHeaveFilter<double>::Config g{};
        double wc = 0.0;
        if (env_double("HB_GODHAVN_FIXED_WC", wc)) {
            g.adaptive = false;
            g.wc_init = wc;
        }
        godhavn_ = heave_baselines::GodhavnHeaveFilter<double>(g);

        heave_baselines::KuchlerHeaveEKF<double>::Config k{};
        // Datasheet white noise of the simulated accelerometer.
        k.accel_noise_std = 1.51e-3 * double(g_std);
        env_double("HB_KUCHLER_ZETA_Q", k.zeta_q);
        env_double("HB_KUCHLER_OMEGA_RW", k.omega_rw);
        kuchler_ = heave_baselines::KuchlerHeaveEKF<double>(k);
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        last_mag_body_ned_ = mag_body_ned;
        have_mag_ = true;
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        (void)temperature_c;
        const Vector3f g = ned_to_mahony_body_(gyr_meas_ned);
        const Vector3f a = ned_to_mahony_body_(acc_meas_ned);
        if (with_mag_ && have_mag_) {
            const Vector3f m = ned_to_mahony_mag_(last_mag_body_ned_);
            frontend_.updateIMUMag(g.x(), g.y(), g.z(), a.x(), a.y(), a.z(),
                                   m.x(), m.y(), m.z(), dt);
        } else {
            frontend_.updateIMU(g.x(), g.y(), g.z(), a.x(), a.y(), a.z(), dt);
        }

        const double a_up = frontend_.verticalWorldAccelUp();
        if constexpr (M == Method::PII) {
            heave_ = frontend_.displacement();
            vel_ = frontend_.velocity();
        } else if constexpr (M == Method::Godhavn) {
            heave_ = godhavn_.update(a_up, dt);
            vel_ = godhavn_.heaveRate();
        } else {
            heave_ = kuchler_.update(a_up, dt);
            vel_ = kuchler_.heaveRate();
        }
    }

    FilterSnapshot snapshot() const override {
        FilterSnapshot s;
        s.disp_est_zu = Vector3f(0.0f, 0.0f, float(heave_));
        s.vel_est_zu = Vector3f(0.0f, 0.0f, float(vel_));

        const auto q = frontend_.quaternionWorldToBody();
        float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
        quat_wb_zu_to_euler_nautical(
            Quaternionf(float(q.w), float(q.x), float(q.y), float(q.z)), roll, pitch, yaw);
        // Same magnetic-frame yaw convention as the PII simulator.
        yaw = with_mag_ ? wrapDeg(yaw + 2.0f * MagSim_WMM::default_declination_deg) : 0.0f;
        s.euler_nautical_deg = Vector3f(roll, pitch, yaw);

        if constexpr (M == Method::Godhavn) {
            s.tuning_applied = float(godhavn_.cutoff());
            s.tuning_target = float(godhavn_.cutoffTarget());
            const double wp = godhavn_.waveFrequencyEstimate();
            s.freq_hz = float(wp / (2.0 * M_PI));
        } else if constexpr (M == Method::Kuchler) {
            s.tuning_applied = float(kuchler_.modeCount());
            const double wp = kuchler_.dominantFrequency();
            s.freq_hz = float(wp / (2.0 * M_PI));
        }
        return s;
    }

private:
    static Vector3f ned_to_mahony_body_(const Vector3f& v_ned) {
        const Vector3f v_zu = ned_to_zu(v_ned);
        return Vector3f(-v_zu.x(), -v_zu.y(), v_zu.z());
    }

    static Vector3f ned_to_mahony_mag_(const Vector3f& v_ned) {
        return Vector3f(v_ned.x(), -v_ned.y(), -v_ned.z());
    }

    static Frontend::Config make_frontend_config_(bool with_mag) {
        Frontend::Config cfg{};
        cfg.gravity_mps2 = g_std;
        cfg.use_mag = with_mag;
        return cfg;
    }

    bool with_mag_ = true;
    bool have_mag_ = false;
    Vector3f last_mag_body_ned_ = Vector3f::Zero();
    Frontend frontend_;
    heave_baselines::GodhavnHeaveFilter<double> godhavn_;
    heave_baselines::KuchlerHeaveEKF<double> kuchler_;
    double heave_ = 0.0;
    double vel_ = 0.0;
};

bool any_gate_failed = false;

void summarize(Method method, const W3dSimulationRunResult& r, float dt) {
    const size_t n_last = static_cast<size_t>(RMS_WINDOW_SEC / dt);
    if (r.errs_z.size() <= n_last) {
        std::cout << "QUALITY_GATE: SKIPPED REASON=record_shorter_than_900s_window RECORD="
                  << r.output_name << "\n";
        return;
    }
    const size_t start = r.errs_z.size() - n_last;
    RMSReport rz, rref;
    for (size_t i = start; i < r.errs_z.size(); ++i) {
        rz.add(r.errs_z[i]);
        rref.add(r.ref_z[i]);
    }
    const float z_pct = 100.0f * rz.rms() / r.wave_params.height;
    std::vector<float> f(r.freq_hist.begin() + long(start), r.freq_hist.end());
    std::cout << "HEAVE_BASELINE method=" << method_name(method)
              << " record=" << r.input_name
              << " Hs=" << r.wave_params.height
              << " z_rms_m=" << rz.rms()
              << " z_pct_hs=" << z_pct
              << " heave_ref_rms_m=" << rref.rms()
              << " freq_hz_median=" << median_vec(f)
              << " final_tuning=" << r.final_tuning_applied << "\n";

    const HeaveLimits lim = method == Method::PII ? PII_LIMITS
        : method == Method::Godhavn ? GODHAVN_LIMITS : KUCHLER_LIMITS;
    const float limit = r.wave_type == WaveType::JONSWAP ? lim.jonswap : lim.pmstokes;
    if (!(z_pct <= limit)) {
        std::cerr << "ERROR: " << method_name(method) << " Z RMS above limit ("
                  << z_pct << "% > " << limit << "%) for " << r.input_name << "\n";
        any_gate_failed = true;
    }
}

template <Method M>
std::optional<W3dSimulationRunResult> run_file(const std::string& file, float dt,
                                               bool with_mag, bool add_noise,
                                               const W3dRandomSeeds& seeds,
                                               bool write_timeseries)
{
    const std::string suffix = std::string("_heave_") + method_name(M);
    return process_wave_file_for_tracker<HeaveBaselineAdapter<M>>(
        file, dt, with_mag, add_noise, 25.0f,
        suffix, suffix + "_nomag", seeds, write_timeseries);
}

} // namespace

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    bool add_noise = true;
    std::vector<Method> methods{Method::PII, Method::Godhavn, Method::Kuchler};
    std::vector<std::string> files;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        } else if (arg == "--input" && i + 1 < argc) {
            files.emplace_back(argv[++i]);
        } else if (arg == "--method" && i + 1 < argc) {
            const std::string m = argv[++i];
            if (m == "pii") methods = {Method::PII};
            else if (m == "godhavn") methods = {Method::Godhavn};
            else if (m == "kuchler") methods = {Method::Kuchler};
            else if (m != "all") {
                std::cerr << "ERROR: unknown method " << m << "\n";
                return 2;
            }
        } else {
            std::cerr << "Usage: " << argv[0]
                      << " [--nomag] [--no-noise] [--method pii|godhavn|kuchler|all]"
                         " [--input PATH]...\n";
            return 2;
        }
    }

    W3dRandomSeeds seeds;
    try {
        seeds = w3d_random_seeds_from_env();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    bool write_timeseries = false;
    if (const char* v = std::getenv("W3D_WRITE_TIMESERIES")) {
        write_timeseries = std::string(v) != "0";
    }
    if (files.empty()) files = collect_wave_data_files(".");

    for (const auto& file : files) {
        for (Method m : methods) {
            std::optional<W3dSimulationRunResult> r;
            switch (m) {
                case Method::PII:
                    r = run_file<Method::PII>(file, dt, with_mag, add_noise, seeds, write_timeseries);
                    break;
                case Method::Godhavn:
                    r = run_file<Method::Godhavn>(file, dt, with_mag, add_noise, seeds, write_timeseries);
                    break;
                case Method::Kuchler:
                    r = run_file<Method::Kuchler>(file, dt, with_mag, add_noise, seeds, write_timeseries);
                    break;
            }
            if (r) summarize(m, *r, dt);
        }
    }
    return any_gate_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
