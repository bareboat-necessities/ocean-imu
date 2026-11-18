#pragma once
// Test-only, single-use intervention at the production handover.
#include "ahrs/FrameConversions.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace startup_oracle {

inline const std::string_view policy = [] {
    const char* name = std::getenv("STARTUP_ORACLE_POLICY");
    if (!name) throw std::invalid_argument("missing startup oracle policy");
    return std::string_view(name);
}();
inline bool previous() { return policy == "previous" || policy == "previous_exact_tilt"; }
inline bool used = false;
inline double epoch = 0;
inline int elapsed_samples = -1;
inline bool trace_sample = false;
inline Eigen::Quaternionf truth = Eigen::Quaternionf::Identity();

inline void observe(double t, const Eigen::Matrix3f& world_to_body_zu) {
    if (used) {
        elapsed_samples = static_cast<int>(std::lround((t - epoch) * 200.0));
        trace_sample = elapsed_samples == 1 || elapsed_samples == 10 || elapsed_samples == 200 ||
                       elapsed_samples == 1000 || elapsed_samples == 2000;
        if (!trace_sample) return;
    } else {
        epoch = t;
    }
    truth = Eigen::Quaternionf(rot_wb_zu_to_bw_ned(world_to_body_zu)).normalized();
}

inline float heading(const Eigen::Quaternionf& q) {
    const auto r = q.toRotationMatrix();
    return std::atan2(r(1, 0), r(0, 0));
}

inline float tilt_error(const Eigen::Quaternionf& q) {
    const Eigen::Vector3f a = q.conjugate() * Eigen::Vector3f::UnitZ();
    const Eigen::Vector3f b = truth.conjugate() * Eigen::Vector3f::UnitZ();
    return std::atan2(a.cross(b).norm(), a.dot(b)) * 180.0f / std::numbers::pi_v<float>;
}

inline void report_stage(const char* stage, const Eigen::Quaternionf& q) {
    std::cout << "GRAVITY_UPDATE sample=" << elapsed_samples << " stage=" << stage
              << " tilt_deg=" << tilt_error(q) << '\n';
}

template<class SetField>
void handover(Eigen::Quaternionf& seed, const Eigen::Vector3f& field, bool learned, SetField set_field) {
    if (used) throw std::logic_error("repeated oracle handover");
    used = true;
    const float before = tilt_error(seed);
    const float yaw_before = heading(seed);
    const bool exact_tilt = policy == "previous_exact_tilt" || policy == "exact_tilt" || policy == "exact_both";
    const bool exact_field = policy == "exact_field" || policy == "exact_both";
    if (exact_tilt) {
        // Left yaw rotation preserves body down. Keep the original heading, including its error.
        seed = (Eigen::AngleAxisf(yaw_before - heading(truth), Eigen::Vector3f::UnitZ()) * truth).normalized();
    }
    const Eigen::Vector3f exact = MagSim_WMM::mag_world_aero(0.0f);
    const float field_before = (field - exact).norm();
    if (exact_field) set_field(exact);
    const float yaw_delta = std::remainder(heading(seed) - yaw_before, 2.0f * std::numbers::pi_v<float>);
    std::cout << "GRAVITY_ORACLE epoch_s=" << epoch << " learned=" << learned
              << " initial_tilt_deg=" << before << " installed_tilt_deg=" << tilt_error(seed)
              << " yaw_change_rad=" << yaw_delta << " initial_field_error_uT=" << field_before
              << " installed_field_error_uT=" << (field - exact).norm() << '\n';
}

} // namespace startup_oracle
