// Time the plain AtomS3R TFG estimator's 200 Hz update and 25 Hz magnetic
// service. Input generation, I/O, and statistics are outside the timed region.
// In tests/kalman_tfg, use "make tfg_update_bench && ./tfg_update_bench".
// To count instructions over the first 512 fully live estimator updates:
//   make -B tfg_update_bench CPPFLAGS=-DTFG_BENCH_CALLGRIND
//   valgrind --tool=callgrind --instr-atstart=no --collect-atstart=no ./tfg_update_bench
#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_ImuUnits.h"
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"

#ifdef TFG_BENCH_CALLGRIND
#include <valgrind/callgrind.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

namespace {
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
using Vec = Eigen::Vector3f;
constexpr float dt = 0.005f;
constexpr float gravity = atoms3r_ical::ImuCalCfg::g_cal_local;
constexpr int warmup = 40000;  // Include delayed magnetic refinement and BA release.
constexpr int samples = 40000;

void report(const char* group, std::vector<double>& times) {
    std::sort(times.begin(), times.end());
    const auto n = times.size();
    const double mean = std::accumulate(times.begin(), times.end(), 0.0) / static_cast<double>(n);
    std::printf("%s steps=%zu mean_us=%.3f p50_us=%.3f p95_us=%.3f p99_us=%.3f\n",
                group, n, mean, times[n/2], times[n*95/100], times[n*99/100]);
}
} // namespace

int main() {
    Fusion fusion;
    Fusion::Config cfg;
    cfg.with_mag = true;
    cfg.sigma_a = Vec::Constant(0.12f);
    cfg.sigma_m = Vec::Constant(0.80f);
    cfg.gyro_noise_density = 0.00135f;
    cfg.gravity_magnitude = gravity;
    cfg.online_tune_warmup_sec = 10.0f;
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    cfg.mag_screen_installation = seastate::common::defaults::SCREEN_SAVED_MAG_INSTALLATION;
    cfg.mag_continuous_hard_iron = cfg.mag_screen_installation;
    fusion.begin(cfg);
    fusion.enableTuner(true);
    fusion.setTauScaledPseudoCadence(true);

    using Clock = std::chrono::steady_clock;
    std::vector<double> ordinary, magnetic, all;
    ordinary.reserve(samples);
    magnetic.reserve(samples/8);
    all.reserve(samples);
    int accepted_acc = 0, accepted_mag = 0;
    for (int i = 0; i < warmup + samples; ++i) {
        const float t = static_cast<float>(i) * dt;
        constexpr float omega = 1.06814f;  // 0.17 Hz roll and heave
        const float roll = 0.12f*std::sin(omega*t);
        const auto Rt = Eigen::AngleAxisf(-roll, Vec::UnitX()).toRotationMatrix();
        const Vec acc = Rt * Vec(0.07f*std::sin(0.83f*t), 0, 0.20f*std::sin(omega*t) - gravity);
        const Vec gyro(0.12f*omega*std::cos(omega*t), 0, 0);
        const Vec mag = Rt * Vec(20, 0, 43);
        const bool mag_due = i % 8 == 0;

        const auto start = Clock::now();
#ifdef TFG_BENCH_CALLGRIND
        constexpr int profiled_steps = 512;
        const bool profile = i >= warmup && i < warmup + profiled_steps;
        if (i == warmup) CALLGRIND_START_INSTRUMENTATION;
        if (profile) CALLGRIND_TOGGLE_COLLECT;
#endif
        fusion.update(dt, gyro, acc, 35.0f);
        if (mag_due) fusion.updateMag(mag);
#ifdef TFG_BENCH_CALLGRIND
        if (profile) CALLGRIND_TOGGLE_COLLECT;
        if (i + 1 == warmup + profiled_steps) CALLGRIND_STOP_INSTRUMENTATION;
#endif
        const auto end = Clock::now();

        if (i < warmup) continue;
        const auto& core = fusion.mekf();
        if (!fusion.isLive() || !core.acc_bias_updates_enabled() || !core.has_magnetic_reference()) {
            std::fprintf(stderr, "sample=%d live=%d acc_bias_enabled=%d magnetic_reference=%d\n",
                         i, fusion.isLive(), core.acc_bias_updates_enabled(), core.has_magnetic_reference());
            return 1;
        }
        const double us = std::chrono::duration<double, std::micro>(end-start).count();
        all.push_back(us);
        (mag_due ? magnetic : ordinary).push_back(us);
        accepted_acc += core.lastAccDiag().accepted;
        accepted_mag += mag_due && core.lastMagDiag().accepted;
    }

    report("all", all);
    report("imu", ordinary);
    report("imu_mag", magnetic);
    const auto& core = fusion.mekf();
    std::printf("accepted_acc=%d accepted_mag=%d final_heave_m=%.8f covariance_trace=%.8f\n",
                accepted_acc, accepted_mag, -core.get_position().z(), core.covariance_full().trace());
    return accepted_acc == samples && accepted_mag == samples/8 &&
           core.get_position().allFinite() && core.covariance_full().allFinite() ? 0 : 1;
}
