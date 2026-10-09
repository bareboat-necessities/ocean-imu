#include <algorithm>
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

    pii          AdaptiveVerticalPII (this repository's non-Kalman observer)
    godhavn      Godhavn 1998 adaptive fourth-order heave filter, cutoff law
                 of Richter et al. 2014 (Eqs. 11-12) as published
    godhavn_ext  the same filter with the bias-instability term and input
                 mean subtraction (not in the papers)
    richter_zd   Richter et al. 2014 zero-displacement filter (Section 3.2)
                 with the same two extensions
    kuchler      Kuchler et al. 2011 FFT-identified harmonic-mode EKF
*/

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include "heave_baselines/GodhavnHeaveFilter.h"
#include "heave_baselines/KuchlerHeaveEKF.h"
#include "kalman_common/SeaStateFusionDefaults.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

constexpr float RMS_WINDOW_SEC = 900.0f;

enum class Method { PII, Godhavn, GodhavnExt, RichterZd, Kuchler };

const char* method_name(Method m) {
    switch (m) {
        case Method::PII: return "pii";
        case Method::Godhavn: return "godhavn";
        case Method::GodhavnExt: return "godhavn_ext";
        case Method::RichterZd: return "richter_zd";
        case Method::Kuchler: return "kuchler";
    }
    return "?";
}

constexpr double kDt = 1.0 / 200.0;
namespace sfd = seastate::common::defaults;

// Vertical reference of the baselines:
//   mahony       the PII observer's AHRS as shipped: sea-state-adaptive gains
//                (time constant ~1 s) and magnetometer correction;
//   proxy        the same Mahony core, IMU-only, at the gains of the private
//                Mahony observer OU-II, OU-III and TFG level their vertical
//                channel with (defaults::STARTUP_PROXY_TWO_KP/KI): the
//                correction corner sits an order of magnitude below the wave
//                band and the magnetometer cannot steer tilt;
//   truth        the measured accelerometer levelled with the record's true
//                attitude, an ideal VRU that isolates the heave filter.
// The record's attitude is loaded in every mode: it scores the tilt error and
// is consumed one sample per IMU step in the runner's order.
enum class Levelling { Mahony, Proxy, Truth };
Levelling g_frontend = Levelling::Mahony;
std::vector<Quaternionf> g_truth_q;
std::vector<float> g_tilt_deg;

const char* frontend_name(Levelling f) {
    switch (f) {
        case Levelling::Mahony: return "mahony";
        case Levelling::Proxy: return "proxy";
        case Levelling::Truth: return "truth";
    }
    return "?";
}

void load_truth_attitude(const std::string& file) {
    g_truth_q.clear();
    g_tilt_deg.clear();
    WaveDataCSVReader reader(file);
    reader.for_each_record([](const Wave_Data_Sample& rec) {
        Quaternionf q(rec.imu.q_wb_zu_w, rec.imu.q_wb_zu_x, rec.imu.q_wb_zu_y, rec.imu.q_wb_zu_z);
        if (!(q.norm() > 1e-6f) || !std::isfinite(q.norm())) q = Quaternionf::Identity();
        // The record's q_wb_zu maps world vectors into the body frame (the
        // runner builds the body magnetometer with it); keep body -> world.
        g_truth_q.push_back(q.normalized().conjugate());
    });
}

bool env_flag(const char* name) {
    const char* s = std::getenv(name);
    return s != nullptr && std::string(s) != "0";
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
constexpr HeaveLimits PII_LIMITS{7.74f, 8.99f};           // worst 7.70 / 8.94 (H8.5)
constexpr HeaveLimits GODHAVN_LIMITS{431.5f, 553.8f};     // worst 429.3 / 551.0 (H8.5)
constexpr HeaveLimits GODHAVN_EXT_LIMITS{21.20f, 26.66f}; // worst 21.09 / 26.52 (H8.5)
constexpr HeaveLimits KUCHLER_LIMITS{18.00f, 20.58f};     // worst 17.91 / 20.47 (H8.5)
constexpr HeaveLimits RICHTER_ZD_LIMITS{17.61f, 22.96f};    // worst 17.52 / 22.84 (H8.5)
// Shared proxy Mahony (--frontend proxy).
constexpr HeaveLimits PII_PROXY_LIMITS{6.02f, 5.90f};          // worst 5.99 / 5.87
constexpr HeaveLimits GODHAVN_PROXY_LIMITS{388.8f, 509.1f};     // worst 386.8 / 506.5
constexpr HeaveLimits GODHAVN_EXT_PROXY_LIMITS{9.82f, 10.70f};  // worst 9.77 / 10.64
constexpr HeaveLimits KUCHLER_PROXY_LIMITS{17.00f, 18.02f};     // worst 16.91 / 17.93
constexpr HeaveLimits RICHTER_ZD_PROXY_LIMITS{5.34f, 5.78f};    // worst 5.31 / 5.75
// Ideal vertical reference (--frontend truth).
constexpr HeaveLimits GODHAVN_TRUTH_LIMITS{386.0f, 495.0f};    // worst 384.0 / 492.4 (H8.5)
constexpr HeaveLimits GODHAVN_EXT_TRUTH_LIMITS{9.03f, 9.70f};  // worst 8.98 / 9.65 (H8.5)
constexpr HeaveLimits KUCHLER_TRUTH_LIMITS{16.62f, 20.23f};    // worst 16.53 / 20.12 (H8.5)
constexpr HeaveLimits RICHTER_ZD_TRUTH_LIMITS{5.21f, 5.36f};   // worst 5.18 / 5.33

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

        // Datasheet white noise of the simulated accelerometer, the same
        // figure both papers take from the sensor data sheet.
        const double accel_noise_std = 1.51e-3 * double(g_std);

        heave_baselines::GodhavnHeaveFilter<double>::Config g{};
        g.noise_density = accel_noise_std * accel_noise_std * kDt;
        double wc = 0.0;
        if (env_double("HB_GODHAVN_FIXED_WC", wc)) {
            g.adaptive = false;
            g.wc_init = wc;
        }
        // godhavn: Richter's Eqs. 11-12 as published.  godhavn_ext adds the
        // bias-instability variance and the input-mean subtraction.
        // richter_zd is Richter's zero-displacement filter with the same two
        // extensions.
        constexpr bool ext = M == Method::GodhavnExt || M == Method::RichterZd;
        g.bias_instability_term = ext || env_flag("HB_GODHAVN_BIAS_TERM");
        g.subtract_input_mean = ext || env_flag("HB_GODHAVN_SUBTRACT_MEAN");
        g.zero_displacement = M == Method::RichterZd || env_flag("HB_GODHAVN_ZD");
        godhavn_ = heave_baselines::GodhavnHeaveFilter<double>(g);

        heave_baselines::KuchlerHeaveEKF<double>::Config k{};
        k.accel_noise_std = accel_noise_std;
        env_double("HB_KUCHLER_ZETA_Q", k.zeta_q);
        env_double("HB_KUCHLER_OMEGA_RW", k.omega_rw);
        k.offset_from_fit = env_flag("HB_KUCHLER_OFFSET_FIT");
        kuchler_ = heave_baselines::KuchlerHeaveEKF<double>(k);
        // The paper feeds the body-z accelerometer and neglects roll/pitch.
        kuchler_body_z_ = env_flag("HB_KUCHLER_BODY_Z");
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
        if (with_mag_ && have_mag_ && g_frontend == Levelling::Mahony) {
            const Vector3f m = ned_to_mahony_mag_(last_mag_body_ned_);
            frontend_.updateIMUMag(g.x(), g.y(), g.z(), a.x(), a.y(), a.z(),
                                   m.x(), m.y(), m.z(), dt);
        } else {
            frontend_.updateIMU(g.x(), g.y(), g.z(), a.x(), a.y(), a.z(), dt);
        }

        if (sample_ >= g_truth_q.size())
            throw std::runtime_error("truth attitude shorter than the IMU stream");
        const Quaternionf& q_true = g_truth_q[sample_];
        double a_up = frontend_.verticalWorldAccelUp();
        if (g_frontend == Levelling::Truth) {
            const Vector3f f_world = q_true * ned_to_zu(acc_meas_ned);
            a_up = double(f_world.z()) - double(g_std);
            g_tilt_deg.push_back(0.0f);
        } else {
            g_tilt_deg.push_back(tilt_error_deg_(q_true));
        }
        ++sample_;
        input_ = a_up;
        if constexpr (M == Method::PII) {
            heave_ = frontend_.displacement();
            vel_ = frontend_.velocity();
        } else if constexpr (M == Method::Godhavn || M == Method::GodhavnExt ||
                             M == Method::RichterZd) {
            heave_ = godhavn_.update(a_up, dt);
            vel_ = godhavn_.heaveRate();
        } else {
            const double a_in = kuchler_body_z_ ? double(a.z()) - double(g_std) : a_up;
            heave_ = kuchler_.update(a_in, dt);
            vel_ = kuchler_.heaveRate();
        }
    }

    FilterSnapshot snapshot() const override {
        FilterSnapshot s;
        s.disp_est_zu = Vector3f(0.0f, 0.0f, float(heave_));
        s.vel_est_zu = Vector3f(0.0f, 0.0f, float(vel_));
        // The levelled vertical acceleration every baseline consumes.
        s.acc_est_zu = Vector3f(0.0f, 0.0f, float(input_));

        const auto q = frontend_.quaternionWorldToBody();
        float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
        quat_wb_zu_to_euler_nautical(
            Quaternionf(float(q.w), float(q.x), float(q.y), float(q.z)), roll, pitch, yaw);
        // Same magnetic-frame yaw convention as the PII simulator.
        yaw = with_mag_ ? wrapDeg(yaw + 2.0f * MagSim_WMM::default_declination_deg) : 0.0f;
        s.euler_nautical_deg = Vector3f(roll, pitch, yaw);

        if constexpr (M == Method::Godhavn || M == Method::GodhavnExt ||
                      M == Method::RichterZd) {
            s.tuning_applied = float(godhavn_.cutoff());
            s.tuning_target = float(godhavn_.cutoffTarget());
            const double wp = godhavn_.dominantFrequency();
            s.freq_hz = float(wp / (2.0 * M_PI));
        } else if constexpr (M == Method::Kuchler) {
            s.tuning_applied = float(kuchler_.modeCount());
            const double wp = kuchler_.dominantFrequency();
            s.freq_hz = float(wp / (2.0 * M_PI));
        }
        return s;
    }

private:
    // Angle between the estimated and true gravity directions in the body
    // frame: the levelling error that leaks horizontal acceleration into the
    // vertical channel, independent of yaw.
    float tilt_error_deg_(const Quaternionf& q_body_to_world) const {
        const auto q = frontend_.quaternionWorldToBody();
        const Vector3f up_est(2.0f * float(q.x * q.z - q.w * q.y),
                              2.0f * float(q.w * q.x + q.y * q.z),
                              float(q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z));
        const Vector3f up_zu = q_body_to_world.conjugate() * Vector3f::UnitZ();
        const Vector3f up_true(-up_zu.x(), -up_zu.y(), up_zu.z());
        const float c = std::clamp(up_est.normalized().dot(up_true), -1.0f, 1.0f);
        return float(std::acos(c) * 180.0 / M_PI);
    }

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
        cfg.use_mag = with_mag && g_frontend == Levelling::Mahony;
        // Fixed Mahony gains: the shared proxy gains, or a sweep override.
        if (g_frontend == Levelling::Proxy) {
            cfg.adapt_mahony_gains = false;
            cfg.mahony_twoKp = sfd::STARTUP_PROXY_TWO_KP;
            cfg.mahony_twoKi = sfd::STARTUP_PROXY_TWO_KI;
        }
        double kp = 0.0, ki = 0.0;
        if (env_double("HB_MAHONY_TWOKP", kp)) {
            cfg.adapt_mahony_gains = false;
            cfg.mahony_twoKp = float(kp);
            cfg.mahony_twoKi = env_double("HB_MAHONY_TWOKI", ki) ? float(ki) : 0.0f;
        }
        return cfg;
    }

    bool with_mag_ = true;
    bool have_mag_ = false;
    Vector3f last_mag_body_ned_ = Vector3f::Zero();
    Frontend frontend_;
    heave_baselines::GodhavnHeaveFilter<double> godhavn_;
    heave_baselines::KuchlerHeaveEKF<double> kuchler_;
    bool kuchler_body_z_ = false;
    double heave_ = 0.0;
    double vel_ = 0.0;
    double input_ = 0.0;
    size_t sample_ = 0;
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
    RMSReport rz, rref, rtilt;
    for (size_t i = start; i < r.errs_z.size(); ++i) {
        rz.add(r.errs_z[i]);
        rref.add(r.ref_z[i]);
        if (i < g_tilt_deg.size()) rtilt.add(g_tilt_deg[i]);
    }
    const float z_pct = 100.0f * rz.rms() / r.wave_params.height;
    std::vector<float> f(r.freq_hist.begin() + long(start), r.freq_hist.end());
    std::cout << "HEAVE_BASELINE method=" << method_name(method)
              << " frontend=" << frontend_name(g_frontend)
              << " record=" << r.input_name
              << " Hs=" << r.wave_params.height
              << " z_rms_m=" << rz.rms()
              << " z_pct_hs=" << z_pct
              << " heave_ref_rms_m=" << rref.rms()
              << " tilt_rms_deg=" << rtilt.rms()
              << " freq_hz_median=" << median_vec(f)
              << " final_tuning=" << r.final_tuning_applied << "\n";

    const HeaveLimits lim =
        g_frontend == Levelling::Truth
        ? (method == Method::Godhavn ? GODHAVN_TRUTH_LIMITS
           : method == Method::GodhavnExt ? GODHAVN_EXT_TRUTH_LIMITS
           : method == Method::RichterZd ? RICHTER_ZD_TRUTH_LIMITS : KUCHLER_TRUTH_LIMITS)
        : g_frontend == Levelling::Proxy
        ? (method == Method::PII ? PII_PROXY_LIMITS
           : method == Method::Godhavn ? GODHAVN_PROXY_LIMITS
           : method == Method::GodhavnExt ? GODHAVN_EXT_PROXY_LIMITS
           : method == Method::RichterZd ? RICHTER_ZD_PROXY_LIMITS : KUCHLER_PROXY_LIMITS)
        : (method == Method::PII ? PII_LIMITS
           : method == Method::Godhavn ? GODHAVN_LIMITS
           : method == Method::GodhavnExt ? GODHAVN_EXT_LIMITS
           : method == Method::RichterZd ? RICHTER_ZD_LIMITS : KUCHLER_LIMITS);
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
    const std::string suffix = std::string("_heave_") + method_name(M) +
        (g_frontend == Levelling::Mahony ? "" : std::string("_") + frontend_name(g_frontend));
    load_truth_attitude(file);
    return process_wave_file_for_tracker<HeaveBaselineAdapter<M>>(
        file, dt, with_mag, add_noise, 25.0f,
        suffix, suffix + "_nomag", seeds, write_timeseries);
}

} // namespace

int main(int argc, char* argv[]) {
    const float dt = float(kDt);
    bool with_mag = true;
    bool add_noise = true;
    std::vector<Method> methods{Method::PII, Method::Godhavn, Method::GodhavnExt,
                                Method::RichterZd, Method::Kuchler};
    std::vector<std::string> files;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        } else if (arg == "--input" && i + 1 < argc) {
            files.emplace_back(argv[++i]);
        } else if (arg == "--frontend" && i + 1 < argc) {
            const std::string f = argv[++i];
            if (f == "truth") g_frontend = Levelling::Truth;
            else if (f == "proxy") g_frontend = Levelling::Proxy;
            else if (f != "mahony") {
                std::cerr << "ERROR: unknown frontend " << f << "\n";
                return 2;
            }
        } else if (arg == "--method" && i + 1 < argc) {
            const std::string m = argv[++i];
            if (m == "pii") methods = {Method::PII};
            else if (m == "godhavn") methods = {Method::Godhavn};
            else if (m == "godhavn_ext") methods = {Method::GodhavnExt};
            else if (m == "richter_zd") methods = {Method::RichterZd};
            else if (m == "kuchler") methods = {Method::Kuchler};
            else if (m != "all") {
                std::cerr << "ERROR: unknown method " << m << "\n";
                return 2;
            }
        } else {
            std::cerr << "Usage: " << argv[0]
                      << " [--nomag] [--no-noise] [--frontend mahony|proxy|truth]"
                         " [--method pii|godhavn|godhavn_ext|richter_zd|kuchler|all]"
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
    // PII carries its own attitude estimator; it has no truth-levelled form.
    if (g_frontend == Levelling::Truth) std::erase(methods, Method::PII);

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
                case Method::GodhavnExt:
                    r = run_file<Method::GodhavnExt>(file, dt, with_mag, add_noise, seeds, write_timeseries);
                    break;
                case Method::RichterZd:
                    r = run_file<Method::RichterZd>(file, dt, with_mag, add_noise, seeds, write_timeseries);
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
