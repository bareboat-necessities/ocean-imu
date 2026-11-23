#pragma once

// The Arduino toolchain normally optimizes for flash. Let GCC unroll the
// few fixed-size numerical kernels on the sample path. Keep the caller's
// floating-point, target and debug settings.
#if defined(__GNUC__) && !defined(__clang__) && defined(__OPTIMIZE_SIZE__)
#define OCEAN_IMU_HOT_LOOP __attribute__((optimize("O2", "unroll-loops")))
#else
#define OCEAN_IMU_HOT_LOOP
#endif
