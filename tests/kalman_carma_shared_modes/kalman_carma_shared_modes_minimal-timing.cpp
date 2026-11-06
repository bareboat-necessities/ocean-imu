/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    kalman_carma_shared_modes_minimal-timing -- host timing of the QMEKF cores
    (no wrapper, no front end): the CARMA3 baseline A41 (21 states), the
    33-state shared-modes core and the 27-state minimal realizations
    (KalmanCarmaSharedModesMinimal.h), on a synthetic 200 Hz stream with
    model commits every 0.1 s (and without, --no-commits) and the
    magnetometer at 25 Hz.  --ftz sets x86 flush-to-zero / denormals-are-zero
    first: the float series of the discretizations underflow into denormals,
    which x86 handles in microcode (tens of times slower), so a host time
    without --ftz overstates the commit cost relative to an FPU without
    denormal support.  Ranking only; not an ESP32-S3 time
    (docs/kalman-carma-shared-modes-minimal.md section 6).

    Build on demand (not part of `make build`):
        make -C tests/kalman_carma_shared_modes kalman_carma_shared_modes_minimal-timing
    Output CSV: core,states,commits,ftz,time_update_us,acc_mag_us,commit_us,total_us
*/

#define EIGEN_NON_ARDUINO

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>

#if defined(__SSE__)
#include <pmmintrin.h>
#include <xmmintrin.h>
#endif

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModesMinimal.h"

namespace csm = ocean_imu::carma_shared_modes;
using V3 = Eigen::Vector3f;

template<class F>
static void run(const char* name, bool commits, bool ftz, int seconds) {
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
        f->time_update(gyr, 0.005f);
        const auto c = std::chrono::steady_clock::now();
        f->measurement_update_acc_only(acc, 35.0f);
        if (k % 8 == 0) f->measurement_update_mag_only(V3(20.0f, 0.0f, 45.0f) + 0.5f * V3(n(rng), n(rng), n(rng)));
        const auto d = std::chrono::steady_clock::now();
        t_com += std::chrono::duration<double, std::micro>(b - a).count();
        t_tu += std::chrono::duration<double, std::micro>(c - b).count();
        t_acc += std::chrono::duration<double, std::micro>(d - c).count();
        pz += f->get_position().z();
    }
    const double nn = static_cast<double>(N);
    std::printf("%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f\n", name, F::NX, commits ? 1 : 0, ftz ? 1 : 0, t_tu / nn, t_acc / nn,
                t_com / nn, (t_tu + t_acc + t_com) / nn);
    if (!std::isfinite(pz)) std::printf("# non-finite position\n");
}

int main(int argc, char** argv) {
    bool ftz = false, commits = true;
    int seconds = 120;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--ftz")) ftz = true;
        else if (!std::strcmp(argv[i], "--no-commits")) commits = false;
        else if (!std::strcmp(argv[i], "--seconds") && i + 1 < argc) seconds = std::atoi(argv[++i]);
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
    using C33 = csm::KalmanCarmaSharedModes<float>;
    using C27 = csm::KalmanCarmaSharedModesMinimal<float, csm::MinimalRealization::Cascade>;
    using Q27 = csm::KalmanCarmaSharedModesMinimal<float, csm::MinimalRealization::Quotient>;
    std::printf("core,states,commits,ftz,time_update_us,acc_mag_us,commit_us,total_us\n");
    for (int rep = 0; rep < 3; ++rep) {   // interleaved repeats: the report takes the median
        run<A41>("A41", commits, ftz, seconds);
        run<C33>("CSM", commits, ftz, seconds);
        run<C27>("CSM27", commits, ftz, seconds);
        run<Q27>("CSM27Q", commits, ftz, seconds);
    }
    return 0;
}
