#pragma once

// Eigen templates take their optimization settings where their headers are
// first included. Limit speed optimization to the estimator and its kernels.
#if defined(__XTENSA__) && defined(__GNUC__) && !defined(__clang__) && defined(__OPTIMIZE_SIZE__)
#define ATOMS3R_BEGIN_FILTER_MATH _Pragma("GCC push_options") _Pragma("GCC optimize (\"O3\", \"unroll-loops\")")
#define ATOMS3R_END_FILTER_MATH _Pragma("GCC pop_options")
#else
#define ATOMS3R_BEGIN_FILTER_MATH
#define ATOMS3R_END_FILTER_MATH
#endif
