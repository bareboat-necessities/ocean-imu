/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    kalman_carma_shared_modes-timing -- host timing of the QMEKF and TFG cores (no
    wrapper, no front end): the CARMA3 baseline A41 (21 states) and the
    27-state shared-modes core (KalmanCarmaSharedModes.h), on a synthetic 200 Hz stream with
    model commits every 0.1 s (and without, --no-commits) and the
    magnetometer at 25 Hz.  --ftz sets x86 flush-to-zero / denormals-are-zero
    first: the float series of the discretizations underflow into denormals,
    which x86 handles in microcode (tens of times slower), so a host time
    without --ftz overstates the commit cost relative to an FPU without
    denormal support.  Ranking only; not an ESP32-S3 time
    (docs/kalman-carma-shared-modes.md).

    Build on demand (not part of `make build`):
        make -C tests/kalman_carma_shared_modes kalman_carma_shared_modes-timing
    Four repetitions rotate each core through every execution position.
    Output CSV includes rep/order so host-order sensitivity can be inspected.
*/

#define EIGEN_NON_ARDUINO

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <vector>
#include <algorithm>
#include <array>

#if defined(__SSE__)
#include <pmmintrin.h>
#include <xmmintrin.h>
#endif

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"

namespace csm = ocean_imu::carma_shared_modes;
using V3 = Eigen::Vector3f;

template<class F>
static void run(const char* name, bool commits, bool ftz, int seconds, bool jitter,
                size_t repetition, size_t position) {
    const V3 sa = V3::Constant(0.0294f), sg = V3::Constant(3.1e-5f), sm = V3::Constant(7.7f);
    auto f = std::make_unique<F>(sa, sg, sm);
    std::mt19937 rng(1);
    std::normal_distribution<float> n(0.0f, 1.0f);
    f->initialize_from_acc_mag(V3(0.0f, 0.0f, -9.8f), V3(20.0f, 0.0f, 45.0f));
    f->set_vertical_carma(1.0f, 0.6f);
    f->set_horizontal_carma(1.0f, 0.5f, 0.3f, 0.0f);
    f->reset_aw_covariance_to_stationary();
    const int N = 200 * seconds;
    double t_tu = 0.0, t_acc = 0.0, t_com = 0.0;
    float pz = 0.0f;
    std::vector<double> times;times.reserve(static_cast<size_t>(N));
    for (int k = 0; k < N; ++k) {
        const float t = static_cast<float>(k) * 0.005f;
        const V3 acc(0.3f * std::sin(1.9f * t) + 0.01f * n(rng), 0.2f * std::cos(1.9f * t) + 0.01f * n(rng),
                     -9.8f + 0.4f * std::cos(1.9f * t) + 0.01f * n(rng));
        const V3 gyr = 0.001f * V3(n(rng), n(rng), n(rng));
        const auto a = std::chrono::steady_clock::now();
        if (commits && k % 20 == 0) {
            const float kap = 1.0f + 0.1f * std::sin(0.01f * t) + 1e-4f * n(rng);
            f->set_vertical_carma(kap, 0.6f * (1.0f + 0.01f * n(rng)));
            f->set_horizontal_carma(kap, 0.5f, 0.3f, 0.0f);
        }
        const auto b = std::chrono::steady_clock::now();
        f->time_update(gyr, 0.005f * (jitter ? 1.0f + 0.25f * std::sin(0.13f*k) : 1.0f));
        const auto c = std::chrono::steady_clock::now();
        f->measurement_update_acc_only(acc, 35.0f);
        if (k % 8 == 0) f->measurement_update_mag_only(V3(20.0f, 0.0f, 45.0f) + 0.5f * V3(n(rng), n(rng), n(rng)));
        const auto d = std::chrono::steady_clock::now();
        t_com += std::chrono::duration<double, std::micro>(b - a).count();
        t_tu += std::chrono::duration<double, std::micro>(c - b).count();
        t_acc += std::chrono::duration<double, std::micro>(d - c).count();
        times.push_back(std::chrono::duration<double,std::micro>(d-a).count());
        pz += f->get_position().z();
    }
    const double nn = static_cast<double>(N);
    std::sort(times.begin(),times.end());
    std::printf("%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%zu,%zu,%zu\n", name, F::NX, commits ? 1 : 0, ftz ? 1 : 0, t_tu / nn, t_acc / nn,
                t_com / nn, (t_tu + t_acc + t_com) / nn,jitter?1:0,times[static_cast<size_t>(.99*(N-1))],times.back(),sizeof(F),repetition,position);
    if (!std::isfinite(pz)) std::printf("# non-finite position\n");
}

int main(int argc, char** argv) {
    bool ftz = false, commits = true, jitter = false;
    int seconds = 120;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--jitter")) jitter = true;
        else if (!std::strcmp(argv[i], "--ftz")) ftz = true;
        else if (!std::strcmp(argv[i], "--no-commits")) commits = false;
        else if (!std::strcmp(argv[i], "--seconds") && i + 1 < argc) seconds = static_cast<int>(std::strtol(argv[++i], nullptr, 10));
    }
    if (seconds <= 0) {
        std::fprintf(stderr, "--seconds must be positive\n");
        return 2;
    }
#if defined(__SSE__)
    if (ftz) {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    }
#else
    ftz = false;
#endif
    using A41 = ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false, true, true>;
    using CSM = csm::KalmanCarmaSharedModes<float>;
    using Runner = void (*)(const char*, bool, bool, int, bool, size_t, size_t);
    const std::array<Runner,4> runners = {run<A41>, run<CSM>,
        run<ocean_imu::kalman::Kalman3D_Wave_TFG_CARMA3<float,4>>,
        run<csm::KalmanTfgCarmaSharedModes<float>>};
    const std::array<const char*,4> names = {"A41", "CSM", "T41", "TCSM"};
    std::printf("core,states,commits,ftz,time_update_us,acc_mag_us,commit_us,total_us,jitter,p99_us,max_us,core_bytes,rep,order\n");
    for (size_t rep = 0; rep < runners.size(); ++rep) {
        for (size_t position = 0; position < runners.size(); ++position) {
            const size_t index = (rep + position) % runners.size();
            runners[index](names[index], commits, ftz, seconds, jitter, rep, position);
        }
    }
    return 0;
}
