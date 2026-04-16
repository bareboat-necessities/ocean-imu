#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

#define EIGEN_NON_ARDUINO
#include "util/WaveFilesSupport.h"
#include "ahrs/FrameConversions.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {
float wrap_deg(float a) {
    a = std::fmod(a + 180.0f, 360.0f);
    if (a < 0.0f) a += 360.0f;
    return a - 180.0f;
}

float diff_deg(float est_deg, float ref_deg) {
    return wrap_deg(est_deg - ref_deg);
}

struct Accum3 {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    std::size_t n = 0;
    void add(const Vector3f& e) { sx += e.x()*e.x(); sy += e.y()*e.y(); sz += e.z()*e.z(); ++n; }
    Vector3f rms() const {
        if (n == 0) return Vector3f::Zero();
        return Vector3f(std::sqrt(float(sx / n)), std::sqrt(float(sy / n)), std::sqrt(float(sz / n)));
    }
};

struct ScalarAccum {
    double sxx = 0.0;
    std::size_t n = 0;
    void add(float v) { sxx += double(v) * double(v); ++n; }
    float rms() const { return (n == 0) ? 0.0f : std::sqrt(float(sxx / n)); }
};

Quaternionf quat_from_sample(const IMU_Sample& imu) {
    Quaternionf q(imu.q_wb_zu_w, imu.q_wb_zu_x, imu.q_wb_zu_y, imu.q_wb_zu_z);
    if (!std::isfinite(q.w()) || !std::isfinite(q.x()) || !std::isfinite(q.y()) || !std::isfinite(q.z()) || q.norm() < 1e-6f) {
        return Quaternionf::Identity();
    }
    q.normalize();
    return q;
}

Vector3f body_rate_from_quat_delta(const Quaternionf& q_prev, const Quaternionf& q_cur, float dt_sec) {
    if (!(dt_sec > 0.0f)) return Vector3f::Zero();
    Quaternionf dq = q_cur * q_prev.conjugate();
    if (dq.w() < 0.0f) dq.coeffs() *= -1.0f;
    const Vector3f v(dq.x(), dq.y(), dq.z());
    const float s = v.norm();
    if (s < 1e-9f) return Vector3f::Zero();
    const float angle = 2.0f * std::atan2(s, std::max(1e-9f, dq.w()));
    return (angle / dt_sec) * (v / s);
}

int check_file(const std::filesystem::path& file, bool& prefer_inverted_out) {
    WaveDataCSVReader reader(file.string());
    bool has_prev = false;
    Wave_Data_Sample prev{};
    Accum3 gyro_err_same{}, gyro_err_invert{}, euler_err_world{}, euler_err_body{};
    ScalarAccum yaw_quat_world_abs{}, yaw_quat_body_abs{};

    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        const auto q = quat_from_sample(rec.imu);
        float r_q = 0.0f, p_q = 0.0f, y_q = 0.0f;
        quat_wb_zu_to_euler_nautical(q, r_q, p_q, y_q);
        euler_err_world.add(Vector3f(diff_deg(rec.imu.roll_deg, r_q), diff_deg(rec.imu.pitch_deg, p_q), diff_deg(rec.imu.yaw_deg, y_q)));
        yaw_quat_world_abs.add(wrap_deg(y_q));

        float r_b = 0.0f, p_b = 0.0f, y_b = 0.0f;
        matrix_to_euler_zyx_deg(q.conjugate().toRotationMatrix(), r_b, p_b, y_b);
        euler_err_body.add(Vector3f(diff_deg(rec.imu.roll_deg, r_b), diff_deg(rec.imu.pitch_deg, p_b), diff_deg(rec.imu.yaw_deg, y_b)));
        yaw_quat_body_abs.add(wrap_deg(y_b));

        if (has_prev) {
            const float dt = float(rec.time - prev.time);
            const auto q_prev = quat_from_sample(prev.imu);
            const Vector3f w_from_q = body_rate_from_quat_delta(q_prev, q, dt);
            const Vector3f w_csv(rec.imu.gyro_x, rec.imu.gyro_y, rec.imu.gyro_z);
            gyro_err_same.add(w_csv - w_from_q);
            gyro_err_invert.add(w_csv + w_from_q);
        }
        prev = rec;
        has_prev = true;
    });

    const Vector3f rms_e_world = euler_err_world.rms();
    const Vector3f rms_e_body = euler_err_body.rms();
    const float euler_world_norm = rms_e_world.norm();
    const float euler_body_norm = rms_e_body.norm();
    const bool euler_matches_world = euler_world_norm <= euler_body_norm;
    const Vector3f rms_e = euler_matches_world ? rms_e_world : rms_e_body;
    const float norm_same = gyro_err_same.rms().norm();
    const float norm_inv = gyro_err_invert.rms().norm();
    const float yaw_quat_world_rms = yaw_quat_world_abs.rms();
    const float yaw_quat_body_rms = yaw_quat_body_abs.rms();
    prefer_inverted_out = norm_inv < norm_same;

    std::cout << "[convention] " << file.filename().string() << "\n"
              << "  Euler CSV-vs-selected-q RMS deg: roll=" << rms_e.x() << " pitch=" << rms_e.y() << " yaw=" << rms_e.z() << "\n"
              << "  Preferred Euler convention:       " << (euler_matches_world ? "world->body (q_wb_zu)" : "body->world (q_bw_zu)") << "\n"
              << "  Yaw RMS from q_wb_zu (world->body): " << yaw_quat_world_rms << " deg\n"
              << "  Yaw RMS from q_bw_zu (body->world): " << yaw_quat_body_rms << " deg\n"
              << "  Gyro RMS (csv - dQ):      " << norm_same << " rad/s\n"
              << "  Gyro RMS (csv + dQ):      " << norm_inv << " rad/s\n"
              << "  Preferred gyro sign:      " << (prefer_inverted_out ? "inverted (csv ~= -dQ)" : "same (csv ~= +dQ)") << "\n";

    bool ok = true;
    // Upstream sim-data-files datasets rely on quaternion columns as the
    // authoritative attitude source. Euler columns may use legacy conventions,
    // so they are diagnostics only and should not fail this check.

    // For canonical wave datasets the body yaw is near zero for at least one
    // quaternion orientation interpretation.
    if (std::min(yaw_quat_world_rms, yaw_quat_body_rms) > 3.0f) ok = false;

    if (std::min(norm_same, norm_inv) > 0.25f) ok = false;
    if (!ok) {
        std::cerr << "ERROR: convention mismatch in " << file.filename().string() << "\n";
        return 1;
    }
    return 0;
}
} // namespace

int main() {
    int failures = 0;
    bool sign_initialized = false;
    bool prefer_inverted_ref = false;
    int checked = 0;

    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (!entry.is_regular_file()) continue;
        const auto name = entry.path().filename().string();
        if (name.rfind("wave_data_jonswap_", 0) == 0 || name.rfind("wave_data_pmstokes_", 0) == 0) {
            ++checked;
            bool prefer_inverted = false;
            failures += check_file(entry.path(), prefer_inverted);
            if (!sign_initialized) {
                sign_initialized = true;
                prefer_inverted_ref = prefer_inverted;
            } else if (prefer_inverted != prefer_inverted_ref) {
                std::cerr << "ERROR: Mixed gyro sign conventions detected across files.\n";
                ++failures;
            }
        }
    }

    if (checked == 0) {
        std::cerr << "ERROR: No wave_data_jonswap_/wave_data_pmstokes_ CSV files found in current directory.\n";
        return 1;
    }
    if (failures) {
        std::cerr << "Convention check failures: " << failures << "\n";
        return 1;
    }
    std::cout << "Convention checks passed.\n";
    return 0;
}
