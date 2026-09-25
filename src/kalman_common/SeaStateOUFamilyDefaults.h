#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Global-scope names of the shared defaults exported by
// SeaStateFusionFilter_OU_II.h and SeaStateFusionFilter_OU_III.h.  The two OU
// headers are never included in the same translation unit (their estimator
// constants and TuneState share names), so each includes this once; the TFG
// orchestrator reads seastate::common::defaults directly and keeps the global
// namespace clean.
//
// Values and their provenance live in SeaStateFusionDefaults.h.  Constants
// that differ between the two OU families (MAX_TUNE_FREQ_HZ, MAX_SIGMA_A,
// PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT, the regularizer bounds, ...) are defined
// in each family's own header.

#include "kalman_common/SeaStateFusionDefaults.h"

// Shared constants
extern const float g_std;

constexpr float ACC_NOISE_FLOOR_SIGMA_DEFAULT = seastate::common::defaults::ACC_NOISE_FLOOR_SIGMA;

constexpr float MIN_FREQ_HZ = seastate::common::defaults::ACCEL_TRACKER_MIN_FREQ_HZ;
constexpr float MAX_FREQ_HZ = seastate::common::defaults::ACCEL_TRACKER_MAX_FREQ_HZ;

constexpr float MIN_TUNE_FREQ_HZ   = seastate::common::defaults::MIN_TUNE_FREQ_HZ;
constexpr float TUNE_FREQ_PRIOR_HZ = seastate::common::defaults::TUNE_FREQ_PRIOR_HZ;

constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = seastate::common::defaults::SIGMA_BAND_LOW_RATIO;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = seastate::common::defaults::SIGMA_BAND_HIGH_RATIO;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = seastate::common::defaults::SIGMA_BAND_MIN_HZ;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = seastate::common::defaults::SIGMA_BAND_MAX_HZ;

constexpr float ADAPT_TAU_SEC          = seastate::common::defaults::ADAPT_TAU_SEC;
constexpr float ADAPT_TAU_SEA_PERIODS  = seastate::common::defaults::ADAPT_TAU_SEA_PERIODS;
constexpr float ADAPT_EVERY_SECS       = seastate::common::defaults::ADAPT_EVERY_SECS;
// Compatibility default of the inner filters for callers that drive them
// directly; the proxy-startup wrappers configure
// defaults::STARTUP_ONLINE_TUNE_WARMUP_SEC.
constexpr float ONLINE_TUNE_WARMUP_SEC = 5.0f;
constexpr float MAG_DELAY_SEC          = seastate::common::defaults::MAG_DELAY_SEC;

constexpr float STARTUP_PROXY_TWO_KP_DEFAULT = seastate::common::defaults::STARTUP_PROXY_TWO_KP;
constexpr float STARTUP_PROXY_TWO_KI_DEFAULT = seastate::common::defaults::STARTUP_PROXY_TWO_KI;

constexpr float ACC_VIBRATION_GUARD_HZ_DEFAULT    = seastate::common::defaults::ACC_VIBRATION_GUARD_HZ;
constexpr int   ACC_VIBRATION_GUARD_POLES_DEFAULT = seastate::common::defaults::ACC_VIBRATION_GUARD_POLES;
constexpr float ACC_VIBRATION_RACC_GAIN_DEFAULT   = seastate::common::defaults::ACC_VIBRATION_RACC_GAIN;

// Frequency smoother dt (the OU wrappers are designed for 200 Hz).
constexpr float FREQ_SMOOTHER_DT = seastate::common::defaults::NOMINAL_IMU_DT_S;

constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S     = seastate::common::defaults::PSEUDO_UPDATE_PERIOD_NOMINAL_S;
constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S        = seastate::common::defaults::PSEUDO_UPDATE_TAU_NOMINAL_S;
constexpr float PSEUDO_UPDATE_TAU_RATIO_DEFAULT    = seastate::common::defaults::PSEUDO_UPDATE_TAU_RATIO;
constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = seastate::common::defaults::PSEUDO_UPDATE_PERIOD_MIN_S;
