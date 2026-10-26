/*
  carma_compile_check -- ESP32-S3 (AtomS3R) build and timing probe for the
  EXPERIMENTAL CARMA heave filter.  Not a deployed sketch: it lives under
  tests/ so the shipped-sketch CI jobs do not pick it up.

  CARMA_CHECK_VARIANT selects the estimator compiled in:
    0  SeaStateFusion_OU_III (shipping baseline, for the size/time delta)
    1  SeaStateFusion_CARMA<KALMANF, 5>        CARMA(5,2), 22 states
    2  SeaStateFusion_CARMA<KALMANF, 4>        CARMA(4,1), 21 states
    3  SeaStateFusion_CARMA<KALMANF, 5, true>  CARMA(5,2) + S_z, 23 states

  On hardware it feeds a synthetic 0.3 Hz heave at 200 Hz with a jittered
  sample interval and prints the mean and maximum micros() per update.
*/

// Arduino-only translation unit: host compilers (the repository's warnings and
// clang-tidy gates compile every tests/*.cpp) see it empty.
#if __has_include(<Arduino.h>)

#include <Arduino.h>
#include <cmath>

#include <ArduinoEigenDense.h>
#include <ArduinoOceanImu.h>

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA.h"

#ifndef CARMA_CHECK_VARIANT
#define CARMA_CHECK_VARIANT 1
#endif

#if CARMA_CHECK_VARIANT == 0
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
#elif CARMA_CHECK_VARIANT == 1
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA<TrackerType::KALMANF, 5>;
#elif CARMA_CHECK_VARIANT == 2
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA<TrackerType::KALMANF, 4>;
#else
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA<TrackerType::KALMANF, 5, true>;
#endif

static Fusion* fusion = nullptr;
static uint32_t n_updates = 0;
static uint64_t sum_us = 0;
static uint32_t max_us = 0;
static float t_sim = 0.0f;
static uint32_t lcg = 12345u;

static float jitter() {
  lcg = lcg * 1664525u + 1013904223u;
  return 1.0f + 0.2f * (float(lcg >> 8) / 16777216.0f - 0.5f);
}

void setup() {
  Serial.begin(115200);
  fusion = new Fusion();
  Fusion::Config cfg;
  cfg.sigma_a = Eigen::Vector3f::Constant(0.0294f);
  cfg.sigma_g = Eigen::Vector3f::Constant(3.1e-5f);
  cfg.sigma_m = Eigen::Vector3f::Constant(7.7f);
  fusion->begin(cfg);
}

void loop() {
  const float dt = 0.005f * jitter();
  t_sim += dt;
  const float w = 2.0f * float(M_PI) * 0.3f;
  const Eigen::Vector3f acc(0.0f, 0.0f, 0.8f * w * w * std::sin(w * t_sim) - 9.80665f);
  const Eigen::Vector3f gyr(0.0f, 0.0f, 0.0f);
  const uint32_t t0 = micros();
  fusion->update(dt, gyr, acc);
  if ((n_updates & 7u) == 0u) fusion->updateMag(Eigen::Vector3f(20.0f, -4.0f, 45.0f));
  const uint32_t us = micros() - t0;
  if (fusion->isLive()) {
    sum_us += us;
    if (us > max_us) max_us = us;
    ++n_updates;
  }
  if (n_updates == 4000u) {
    Serial.printf("variant=%d updates=%u mean_us=%.1f max_us=%u heave=%.3f\n",
                  CARMA_CHECK_VARIANT, n_updates, double(sum_us) / n_updates, max_us,
                  double(fusion->raw().mekf().get_position().z()));
    n_updates = 0; sum_us = 0; max_us = 0;
  }
}

#endif  // __has_include(<Arduino.h>)
