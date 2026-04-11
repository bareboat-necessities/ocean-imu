/*
   Copyright (c) 2025-2026  Mikhail Grushinskiy
*/

#include <iostream>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <cstdlib>

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const float g_std = 9.80665f;
const float MAG_DELAY_SEC = 5.0f;
const int   IMU_RATE_HZ = 200;
const int   MAG_RATE_HZ = 100;
const int   MAG_UPDATE_STRIDE = IMU_RATE_HZ / MAG_RATE_HZ;

#include "ahrs/KalmanQMEKF.h"
#include "util/WaveFilesSupport.h"
#include "ahrs/FrameConversions.h"

using Eigen::Vector3f;
using Eigen::Quaternionf;
using Eigen::Matrix3f;

// Noise model
bool add_noise = true;

struct NoiseModel {
    std::default_random_engine rng;
    std::normal_distribution<float> dist;
    Vector3f bias;
};

NoiseModel make_noise_model(float sigma, float bias_range, unsigned seed) {
    NoiseModel m{
        std::default_random_engine(seed),
        std::normal_distribution<float>(0.0f, sigma),
        Vector3f::Zero()
    };
    std::uniform_real_distribution<float> ub(-bias_range, bias_range);
    m.bias = Vector3f(ub(m.rng), ub(m.rng), ub(m.rng));
    return m;
}

// Measurement model: meas = truth + bias + white
Vector3f apply_noise(const Vector3f& v, NoiseModel& m) {
    return v + m.bias + Vector3f(m.dist(m.rng), m.dist(m.rng), m.dist(m.rng));
}

struct OutputRow {
    double t{};

    // Reference Euler (deg, nautical: world->body in ENU/Z-up)
    float roll_ref{}, pitch_ref{}, yaw_ref{};

    // Raw IMU inputs from file (nautical body frame)
    float acc_bx{}, acc_by{}, acc_bz{};
    float gyro_x{}, gyro_y{}, gyro_z{};

    // Kalman estimates (deg, nautical)
    float roll_est{}, pitch_est{}, yaw_est{};

    // Noisy IMU actually fed to filter before frame conversion
    float acc_noisy_x{}, acc_noisy_y{}, acc_noisy_z{};
    float gyro_noisy_x{}, gyro_noisy_y{}, gyro_noisy_z{};
};

struct AhrsFailureLimits {
    float err_limit_roll_deg = 10.2f;
    float err_limit_pitch_deg = 10.2f;
    float err_limit_yaw_deg = 30.0f;
    float min_processing_hz = 200.0f;
};

static float wrap_deg(float a) {
    a = std::fmod(a + 180.0f, 360.0f);
    if (a < 0.0f) a += 360.0f;
    return a - 180.0f;
}

static float diff_deg(float est_deg, float ref_deg) {
    return wrap_deg(est_deg - ref_deg);
}

static bool quat_is_placeholder_identity(const IMU_Sample& imu) {
    constexpr float eps = 1e-6f;
    return std::fabs(imu.q_wb_zu_w - 1.0f) < eps &&
           std::fabs(imu.q_wb_zu_x) < eps &&
           std::fabs(imu.q_wb_zu_y) < eps &&
           std::fabs(imu.q_wb_zu_z) < eps;
}

// Prefer the stored quaternion truth when present.
// Fall back to Euler-derived quaternion for legacy CSVs that lack quaternion columns.
static Quaternionf reference_quat_wb_zu(const IMU_Sample& imu) {
    const bool placeholder_q = quat_is_placeholder_identity(imu);
    const bool euler_nontrivial =
        std::fabs(imu.roll_deg) > 1e-4f ||
        std::fabs(imu.pitch_deg) > 1e-4f ||
        std::fabs(imu.yaw_deg) > 1e-4f;

    if (placeholder_q && euler_nontrivial) {
        return quat_from_euler(imu.roll_deg, imu.pitch_deg, imu.yaw_deg);
    }

    return quat_wb_zu_from_csv(
        imu.q_wb_zu_w,
        imu.q_wb_zu_x,
        imu.q_wb_zu_y,
        imu.q_wb_zu_z
    );
}

static void print_summary_and_fail_if_needed(const std::string& output_name,
                                             const std::vector<OutputRow>& rows,
                                             float dt,
                                             bool with_mag,
                                             double elapsed_sec,
                                             const AhrsFailureLimits& limits)
{
    constexpr float RMS_WINDOW_SEC = 60.0f;
    const int n_last = static_cast<int>(RMS_WINDOW_SEC / dt);
    if (rows.size() <= static_cast<size_t>(n_last)) return;

    const size_t start = rows.size() - static_cast<size_t>(n_last);
    float roll_ss = 0.0f, pitch_ss = 0.0f, yaw_ss = 0.0f;

    for (size_t i = start; i < rows.size(); ++i) {
        const float er = diff_deg(rows[i].roll_est, rows[i].roll_ref);
        const float ep = diff_deg(rows[i].pitch_est, rows[i].pitch_ref);
        const float ey = diff_deg(rows[i].yaw_est, rows[i].yaw_ref);
        roll_ss += er * er;
        pitch_ss += ep * ep;
        yaw_ss += ey * ey;
    }

    const float n = static_cast<float>(rows.size() - start);
    const float roll_rms = std::sqrt(roll_ss / n);
    const float pitch_rms = std::sqrt(pitch_ss / n);
    const float yaw_rms = std::sqrt(yaw_ss / n);

    const double processing_hz = (elapsed_sec > 0.0)
        ? static_cast<double>(rows.size()) / elapsed_sec
        : 0.0;
    const double realtime_x = (elapsed_sec > 0.0)
        ? (static_cast<double>(rows.size()) * dt) / elapsed_sec
        : 0.0;

    std::cout << "=== Last 60 s RMS summary for " << output_name << " ===\n";
    std::cout << "Angles RMS (deg): Roll=" << roll_rms
              << " Pitch=" << pitch_rms
              << " Yaw=" << yaw_rms << "\n";
    std::cout << "Performance: elapsed=" << elapsed_sec
              << " s, samples=" << rows.size()
              << ", throughput=" << processing_hz
              << " samples/s, realtime_x=" << realtime_x << "\n";

    if (roll_rms > limits.err_limit_roll_deg) {
        std::cerr << "ERROR: Roll RMS above limit (" << roll_rms << " deg > "
                  << limits.err_limit_roll_deg << " deg). Failing.\n";
        std::exit(EXIT_FAILURE);
    }
    if (pitch_rms > limits.err_limit_pitch_deg) {
        std::cerr << "ERROR: Pitch RMS above limit (" << pitch_rms << " deg > "
                  << limits.err_limit_pitch_deg << " deg). Failing.\n";
        std::exit(EXIT_FAILURE);
    }
    if (with_mag && yaw_rms > limits.err_limit_yaw_deg) {
        std::cerr << "ERROR: Yaw RMS above limit (" << yaw_rms << " deg > "
                  << limits.err_limit_yaw_deg << " deg). Failing.\n";
        std::exit(EXIT_FAILURE);
    }
    if (processing_hz < limits.min_processing_hz) {
        std::cerr << "ERROR: Processing throughput below limit (" << processing_hz
                  << " samples/s < " << limits.min_processing_hz
                  << " samples/s). Failing.\n";
        std::exit(EXIT_FAILURE);
    }
}

void process_wave_file(const std::string& filename, float dt, bool with_mag) {
    auto parsed = WaveFileNaming::parse_to_params(filename);
    if (!parsed) return;

    auto [kind, type, wp] = *parsed;
    if (kind != FileKind::Data) return;
    if (!(type == WaveType::JONSWAP || type == WaveType::PMSTOKES)) return;

    std::cout << "Processing " << filename
              << " (" << EnumTraits<WaveType>::to_string(type)
              << ") with_mag=" << (with_mag ? "true" : "false")
              << ", noise=" << (add_noise ? "true" : "false") << "\n";

    WaveDataCSVReader reader(filename);

    // accel update uses acc / g_std
    const Vector3f sigma_a(0.055f,   0.055f,   0.055f);
    const Vector3f sigma_g(0.00074f, 0.00074f, 0.00074f);
    const Vector3f sigma_m(0.015f,   0.015f,   0.015f);

    QuaternionMEKF<float, true> mekf(sigma_a, sigma_g, sigma_m, 0.0001f, 0.1f, 1e-06f);

    // Fixed world magnetic field in aerospace NED.
    const Vector3f B_world_ned = MagSim_WMM::mag_world_aero();
    mekf.set_mag_world_ref(B_world_ned);

    // Noise models
    static NoiseModel accel_noise = make_noise_model(0.04f,  0.05f,   1234);
    static NoiseModel gyro_noise  = make_noise_model(0.001f, 0.0004f, 5678);

    bool first = true;
    bool mag_enabled = false;
    int iter = 0;

    std::vector<OutputRow> rows;
    const auto t0 = std::chrono::steady_clock::now();

    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        ++iter;

        // Raw body-frame IMU in nautical ENU/Z-up from file
        Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        Vector3f gyr_b(rec.imu.gyro_x, rec.imu.gyro_y, rec.imu.gyro_z);

        // Apply optional additive noise/bias
        Vector3f acc_noisy = acc_b;
        Vector3f gyr_noisy = gyr_b;
        if (add_noise) {
            acc_noisy = apply_noise(acc_b, accel_noise);
            gyr_noisy = apply_noise(gyr_b, gyro_noise);
        }

        // Convert sensor feeds to filter body frame: NED
        const Vector3f acc_f = zu_to_ned(acc_noisy);
        const Vector3f gyr_f = zu_to_ned(gyr_noisy);

        // Truth attitude comes from the same source used to synthesize mag: q_wb_zu.
        const Quaternionf q_ref_wb_zu = reference_quat_wb_zu(rec.imu);

        float r_ref_n = 0.0f;
        float p_ref_n = 0.0f;
        float y_ref_n = 0.0f;
        quat_wb_zu_to_euler_nautical(q_ref_wb_zu, r_ref_n, p_ref_n, y_ref_n);

        // Simulated magnetometer from quaternion truth, then convert body ENU -> body NED.
        Vector3f mag_f = Vector3f::Zero();
        if (with_mag) {
            const Vector3f mag_b_enu = mag_body_from_quat_wb_zu(q_ref_wb_zu);
            mag_f = zu_to_ned(mag_b_enu);
        }

        // Initialize from accel only.
        if (first) {
            mekf.initialize_from_acc(acc_f / g_std);
            first = false;
        }

        // Propagation
        mekf.time_update(gyr_f, dt);

        // Updates
        if (with_mag && rec.time >= MAG_DELAY_SEC) {
            if (!mag_enabled) {
                mekf.measurement_update_mag_only(mag_f);
                mag_enabled = true;
            }

            if (iter % MAG_UPDATE_STRIDE == 0) {
                mekf.measurement_update_mag_only(mag_f);
            }

            mekf.measurement_update_acc_only(acc_f / g_std);
        } else {
            mekf.measurement_update_acc_only(acc_f / g_std);
        }

        // Filter quaternion is BODY->WORLD in NED.
        const auto coeffs = mekf.quaternion(); // [x,y,z,w]
        const Quaternionf q_bw_ned(coeffs(3), coeffs(0), coeffs(1), coeffs(2));

        float r_est = 0.0f;
        float p_est = 0.0f;
        float y_est = 0.0f;
        quat_to_euler_nautical(q_bw_ned, r_est, p_est, y_est);

        rows.push_back({
            rec.time,
            r_ref_n, p_ref_n, y_ref_n,
            rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz,
            rec.imu.gyro_x, rec.imu.gyro_y, rec.imu.gyro_z,
            r_est, p_est, y_est,
            acc_noisy.x(), acc_noisy.y(), acc_noisy.z(),
            gyr_noisy.x(), gyr_noisy.y(), gyr_noisy.z()
        });
    });

    // Output file
    std::string outname = filename;
    auto pos_prefix = outname.find("wave_data_");
    if (pos_prefix != std::string::npos) {
        outname.replace(pos_prefix, std::string("wave_data_").size(), "qmekf_");
    }
    auto pos_ext = outname.rfind(".csv");
    if (pos_ext != std::string::npos) {
        outname.insert(pos_ext, with_mag ? "_kalman" : "_kalman_nomag");
    } else {
        outname += with_mag ? "_kalman.csv" : "_kalman_nomag.csv";
    }

    std::ofstream ofs(outname);
    ofs << "time,"
        << "roll_ref,pitch_ref,yaw_ref,"
        << "acc_bx,acc_by,acc_bz,"
        << "gyro_x,gyro_y,gyro_z,"
        << "roll_est,pitch_est,yaw_est,"
        << "acc_noisy_x,acc_noisy_y,acc_noisy_z,"
        << "gyro_noisy_x,gyro_noisy_y,gyro_noisy_z\n";

    for (const auto& r : rows) {
        ofs << r.t << ","
            << r.roll_ref << "," << r.pitch_ref << "," << r.yaw_ref << ","
            << r.acc_bx   << "," << r.acc_by    << "," << r.acc_bz << ","
            << r.gyro_x   << "," << r.gyro_y    << "," << r.gyro_z << ","
            << r.roll_est << "," << r.pitch_est << "," << r.yaw_est << ","
            << r.acc_noisy_x << "," << r.acc_noisy_y << "," << r.acc_noisy_z << ","
            << r.gyro_noisy_x << "," << r.gyro_noisy_y << "," << r.gyro_noisy_z
            << "\n";
    }
    ofs.close();

    std::cout << "Wrote " << outname << "\n";

    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

    static constexpr AhrsFailureLimits FAIL_LIMITS{};
    print_summary_and_fail_if_needed(outname, rows, dt, with_mag, elapsed_sec, FAIL_LIMITS);
}

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;

    bool with_mag = true;
    add_noise = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        }
    }

    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", mag_delay=" << MAG_DELAY_SEC
              << " sec, noise=" << (add_noise ? "true" : "false") << "\n";

    std::vector<std::string> files;
    for (auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;

        std::string fname = entry.path().string();
        if (fname.find("wave_data_") == std::string::npos) continue;

        if (auto kind = WaveFileNaming::parse_kind_only(fname);
            kind && *kind == FileKind::Data) {
            files.push_back(fname);
        }
    }

    std::sort(files.begin(), files.end());

    for (const auto& fname : files) {
        process_wave_file(fname, dt, with_mag);
    }

    return 0;
}
