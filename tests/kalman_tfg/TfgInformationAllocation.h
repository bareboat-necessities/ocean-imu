#pragma once

// Test-only passive readout of the already-retained accelerometer H, K and r.
// No estimator state, covariance, scheduling, noise or quality limit is changed.
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <type_traits>

class TfgInformationAllocation {
public:
    TfgInformationAllocation() {
        const char* value = std::getenv("TFG_ALLOCATION_TRACE");
        enabled_ = value != nullptr && std::string_view(value) == "1";
    }

    template<class Fusion>
    void observe(float dt, const Fusion& fusion) {
        if (!enabled_) return;
        time_ += static_cast<double>(dt);
        ++updates_;
        if (updates_ % 10U != 0U) return;  // 20 Hz at the pinned 200 Hz input.
        Row row{};
        row.time = time_;
        auto& v = row.values;
        v[1] = fusion.isLive() ? 1.0 : 0.0;
        v[2] = fusion.getTauApplied();
        v[3] = fusion.getSigmaApplied();
        v[4] = fusion.getRSTarget();
        v[5] = fusion.getRSApplied();
        v[6] = fusion.getRSTarget() <= 0.1500001f ? 1.0 : 0.0;
        const auto& diagnostic = fusion.mekf().lastAccDiag();
        if (fusion.isLive() && diagnostic.accepted) {
            using Kf = std::remove_cvref_t<decltype(fusion.mekf())>;
            const auto h = diagnostic.H.template cast<double>().eval();
            const auto k = diagnostic.K.template cast<double>().eval();
            const auto r = diagnostic.r.template cast<double>().eval();
            const auto gt = (h.template block<3, 3>(0, Kf::OFF_PHI) *
                             k.template block<3, 3>(Kf::OFF_PHI, 0)).eval();
            const auto ga = (h.template block<3, 3>(0, Kf::OFF_AW) *
                             k.template block<3, 3>(Kf::OFF_AW, 0)).eval();
            const auto gb = (h.template block<3, 3>(0, Kf::OFF_BA) *
                             k.template block<3, 3>(Kf::OFF_BA, 0)).eval();
            const auto ut = (gt * r).eval();
            const auto ua = (ga * r).eval();
            const auto ub = (gb * r).eval();
            v[0] = 1.0;
            v[7] = diagnostic.nis;
            v[8] = r.squaredNorm();
            v[9] = gt.trace(); v[10] = ga.trace(); v[11] = gb.trace();
            v[12] = gt.squaredNorm(); v[13] = ga.squaredNorm(); v[14] = gb.squaredNorm();
            v[15] = ut.squaredNorm(); v[16] = ua.squaredNorm(); v[17] = ub.squaredNorm();
            v[18] = ut.dot(ua); v[19] = ut.dot(ub); v[20] = ua.dot(ub);
            v[21] = (ut + ua + ub).squaredNorm();
            v[22] = r.dot(ut); v[23] = r.dot(ua); v[24] = r.dot(ub);
            v[25] = (gt + ga + gb - h * k).squaredNorm();
        }
        for (double value : v) {
            if (!std::isfinite(value)) { ++invalid_; return; }
        }
        rows_.push_back(row);
        for (std::size_t i = 0; i < names_.size(); ++i) sums_[i] += v[i];
        prune_();
    }

    void finish() {
        if (!enabled_ || finished_) return;
        finished_ = true;
        prune_();
        const auto precision = std::cout.precision();
        std::cout << std::setprecision(17)
                  << "TFG_ALLOCATION_SUMMARY samples=" << rows_.size()
                  << " invalid_samples=" << invalid_
                  << " end_sec=" << time_ << " window_sec=900 sample_stride=10";
        const double denominator = static_cast<double>(rows_.size());
        for (std::size_t i = 0; i < names_.size(); ++i) {
            std::cout << ' ' << names_[i] << '='
                      << (denominator > 0.0 ? static_cast<double>(sums_[i]) / denominator : 0.0);
        }
        std::cout << '\n';
        std::cout.precision(precision);
    }

private:
    static constexpr std::array<const char*, 26> names_{{
        "accepted", "live", "tau", "sigma_aw", "rs_target", "rs_applied", "rs_floor",
        "nis", "residual_energy", "theta_trace", "aw_trace", "ba_trace",
        "theta_fro2", "aw_fro2", "ba_fro2", "theta_energy", "aw_energy", "ba_energy",
        "theta_aw_cross", "theta_ba_cross", "aw_ba_cross", "net_energy",
        "theta_r_dot", "aw_r_dot", "ba_r_dot", "closure_sq"
    }};
    struct Row { double time{}; std::array<double, names_.size()> values{}; };
    void prune_() {
        while (!rows_.empty() && rows_.front().time <= time_ - 900.0) {
            for (std::size_t i = 0; i < names_.size(); ++i) sums_[i] -= rows_.front().values[i];
            rows_.pop_front();
        }
    }
    bool enabled_{false};
    bool finished_{false};
    double time_{0.0};
    std::uint64_t updates_{0};
    std::uint64_t invalid_{0};
    std::deque<Row> rows_{};
    std::array<long double, names_.size()> sums_{};
};
