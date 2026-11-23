#pragma once

// The Arduino toolchain normally optimizes for flash. Let GCC unroll the
// few fixed-size numerical kernels on the sample path. Keep the caller's
// floating-point, target and debug settings.
#if defined(__GNUC__) && !defined(__clang__) && defined(__OPTIMIZE_SIZE__)
#define OCEAN_IMU_HOT_LOOP __attribute__((optimize("O2", "unroll-loops")))
#else
#define OCEAN_IMU_HOT_LOOP
#endif

// The AtomS3R release has a 16 KiB instruction cache. These larger loops run
// together on every sample and otherwise evict each other from flash cache.
#if defined(OCEAN_IMU_CARMA_SAMPLE_IRAM) && defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
#include <esp_attr.h>
#define OCEAN_IMU_SAMPLE_IRAM IRAM_ATTR
#else
#define OCEAN_IMU_SAMPLE_IRAM
#endif
