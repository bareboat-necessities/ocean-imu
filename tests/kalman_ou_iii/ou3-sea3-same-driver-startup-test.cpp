// Canonical non-promoting complete-SEA3 same-driver startup bridge.
//
// Input is produced by ou3_sea3_validated_continuum_float_trace.py. Every
// binary32 accelerometer sample is uniquely implied by an outward enclosure of
// ONE exact legal continuum SEA3 member. This executable contains no sea
// generator, quadrature, replay, finite harmonics, or independent sample
// choices: it feeds those exact implementation inputs to unchanged shipping
// code and records the reachable Live covariance and applied schedule.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;

namespace {
constexpr float DT = 0.005f;
constexpr int PREHISTORY_SAMPLES = 12000;
constexpr int WORD_SAMPLES = 601;
constexpr float RS_X_FACTOR = 0.72f;
constexpr float RS_Y_FACTOR = 0.72f;
constexpr float RS_MIN = 0.15f;
constexpr float RS_MAX = 100.0f;

bool finite_positive(float x) { return std::isfinite(x) && x > 0.0f; }

void print_matrix_json(std::ostream& os, const auto& M) {
    os << '[';
    for (Eigen::Index i = 0; i < M.rows(); ++i) {
        if (i) os << ',';
        os << '[';
        for (Eigen::Index j = 0; j < M.cols(); ++j) {
            if (j) os << ',';
            os << std::setprecision(10) << M(i,j);
        }
        os << ']';
    }
    os << ']';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: ou3-sea3-same-driver-startup-test <validated-specific-force-z.txt>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open validated complete-SEA3 input trace\n";
        return 2;
    }
    std::vector<float> sfz;
    float sample = 0.0f;
    while (in >> sample) {
        if (!std::isfinite(sample)) {
            std::cerr << "nonfinite validated input\n";
            return 2;
        }
        sfz.push_back(sample);
    }
    if (sfz.size() != static_cast<std::size_t>(PREHISTORY_SAMPLES + WORD_SAMPLES)) {
        std::cerr << "validated input trace length mismatch: " << sfz.size() << "\n";
        return 2;
    }

    using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
    Filter f(false);
    f.setWithMag(false);
    f.setOnlineTuneWarmupSec(10.0f);
    f.initialize(Eigen::Vector3f::Constant(0.2f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.3f));
    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
    auto acc = [&](int k) {
        return Eigen::Vector3f(0.0f, 0.0f, sfz[static_cast<std::size_t>(k)]);
    };
    f.initialize_from_acc(acc(0));

    int first_ready = -1;
    float max_guard = 0.0f;
    for (int k = 0; k < PREHISTORY_SAMPLES; ++k) {
        f.updateFrontEnd(DT, gyro, acc(k));
        max_guard = std::max(max_guard, f.accelVibrationGuardEngagement());
        if (first_ready < 0 && f.isTunerReady()) first_ready = k;
    }
    if (first_ready < 0 || !f.isTunerReady() || !f.wavePeriodUsable() ||
        !f.startupProxyInitialized()) {
        std::cerr << "exact continuum SEA3 member did not reach shipping TunerReady\n";
        return 1;
    }
    if (max_guard != 0.0f) {
        std::cerr << "exact continuum member left dormant vibration-guard branch\n";
        return 1;
    }

    const float period_live = f.getWavePeriodSec();
    const float tau_live = f.getTauApplied();
    const float sigma_live = f.getSigmaApplied();
    const float rs_live = f.getRSApplied();
    if (!finite_positive(period_live) || !finite_positive(tau_live) ||
        !finite_positive(sigma_live) || !finite_positive(rs_live)) {
        std::cerr << "shipping schedule invalid at Live handoff\n";
        return 1;
    }
    if (f.getRSLaw() != RSAdaptationLaw::SpectralMSE) {
        std::cerr << "complete-SEA3 proof requires deployed SpectralMSE R_S law\n";
        return 1;
    }

    f.goLive(f.startupProxyQuat(), 0.035f, 1.5708f);
    if (!f.isAdaptiveLive()) {
        std::cerr << "shipping goLive failed\n";
        return 1;
    }
    const auto P0 = f.mekf().covariance_full();
    const auto P0_ldlt = P0.selfadjointView<Eigen::Lower>().ldlt();
    if (!P0.allFinite() || P0_ldlt.info() != Eigen::Success ||
        P0_ldlt.vectorD().minCoeff() <= 0.0f) {
        std::cerr << "reachable shipping Live covariance is not SPD\n";
        return 1;
    }

    struct Row {
        int k;
        float pseudo_period;
        float tau;
        float sigma;
        float rs_base;
        float rs_x;
        float rs_y;
        float rs_z;
    };
    std::vector<Row> word;
    word.reserve(WORD_SAMPLES);
    for (int k = 0; k < WORD_SAMPLES; ++k) {
        const int source_k = PREHISTORY_SAMPLES + k;
        f.updateTime(DT, gyro, acc(source_k), 35.0f);
        const float base = std::min(std::max(f.getRSApplied(), RS_MIN), RS_MAX);
        const float pseudo_period = f.getPseudoUpdatePeriodSec();
        if (!finite_positive(base) || !finite_positive(pseudo_period)) {
            std::cerr << "invalid actual-applied R_S/cadence inside word\n";
            return 1;
        }
        // SpectralMSE already contains realized T_S, so shipping applies no
        // additional cadence normalization. Thus these are the exact stds
        // passed by apply_RS_tune_() to set_RS_noise().
        word.push_back({k, pseudo_period, f.getTauApplied(), f.getSigmaApplied(),
                        f.getRSApplied(), base * RS_X_FACTOR,
                        base * RS_Y_FACTOR, base});
    }

    std::cout << std::setprecision(10) << '{';
    std::cout << "\"qualification\":\"COMPLETE_SEA3_SAME_DRIVER_SHIPPING_TRACE_V1\",";
    std::cout << "\"canonical_source\":\"COMPLETE_SEA3_NORMAL_LIVE_WORD\",";
    std::cout << "\"same_validated_continuum_member\":true,";
    std::cout << "\"trajectory_replay_used\":false,";
    std::cout << "\"finite_harmonic_source_used\":false,";
    std::cout << "\"first_tuner_ready_sample\":" << first_ready << ',';
    std::cout << "\"first_tuner_ready_time_s\":" << (first_ready + 1) * DT << ',';
    std::cout << "\"wave_period_at_live_s\":" << period_live << ',';
    std::cout << "\"tau_at_live_s\":" << tau_live << ',';
    std::cout << "\"sigma_at_live_mps2\":" << sigma_live << ',';
    std::cout << "\"RS_base_at_live\":" << rs_live << ',';
    std::cout << "\"vibration_guard_max\":" << max_guard << ',';
    std::cout << "\"live_covariance\":";
    print_matrix_json(std::cout, P0);
    std::cout << ",\"word_schedule\":[";
    for (std::size_t i = 0; i < word.size(); ++i) {
        if (i) std::cout << ',';
        const auto& r = word[i];
        std::cout << '{'
                  << "\"k\":" << r.k << ','
                  << "\"pseudo_period_s\":" << r.pseudo_period << ','
                  << "\"tau_applied_s\":" << r.tau << ','
                  << "\"sigma_applied_mps2\":" << r.sigma << ','
                  << "\"RS_base_applied\":" << r.rs_base << ','
                  << "\"RS_std_xyz\":[" << r.rs_x << ',' << r.rs_y << ',' << r.rs_z << ']'
                  << '}';
    }
    std::cout << "]}\n";
    return 0;
}
