#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Production-default regression limits measured on main a053a1d2 (v2.7.6).
// Eight pinned v1.2.3 vessel-RAO records, default seeds, trailing 900 s.
// Worst score + 0.5%, rounded upward to 1e-5, capped by an existing
// stricter production limit. Historical --regression bars stay separate.
// Evidence and reproduction: docs/filter-simulator-gates.md.
#include "util/W3dSimCommon.h"
namespace ocean_imu::sim::current_limits {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 6.66210f, // measured 6.62894869
    .err_limit_percent_z_pmstokes  = 6.46219f, // measured 6.43003035
    .err_limit_yaw_deg             = 0.72934f, // measured 0.725703537
    .err_limit_roll_deg            = 0.22667f, // measured 0.225542024
    .err_limit_pitch_deg           = 0.21278f, // measured 0.211719394
    .err_limit_percent_3d_jonswap  = 12.20058f, // measured 12.1398802
    .err_limit_percent_3d_pmstokes = 12.02862f, // measured 11.968771
    .acc_z_bias_percent            = 4.58751f, // measured 4.56468
    .bias_3d_percent               = 61.99132f, // measured 61.6829
    .gyro_bias_3d_percent          = 12.64713f, // measured 12.5842
};

inline constexpr W3dFailureLimits ou3{
    .err_limit_percent_z_jonswap   = 4.27053f, // measured 4.24928093
    .err_limit_percent_z_pmstokes  = 4.23388f, // measured 4.21280766
    .err_limit_yaw_deg             = 0.82829f, // measured 0.824165642
    .err_limit_roll_deg            = 0.32284f, // measured 0.321227342
    .err_limit_pitch_deg           = 0.13815f, // measured 0.137457415
    .err_limit_percent_3d_jonswap  = 9.73471f, // measured 9.68627834
    .err_limit_percent_3d_pmstokes = 9.69867f, // measured 9.65040874
    .acc_z_bias_percent            = 4.35523f, // measured 4.33356
    .bias_3d_percent               = 69.10592f, // measured 68.7621
    .gyro_bias_3d_percent          = 12.25508f, // measured 12.1941
};

inline constexpr W3dFailureLimits tfg{
    .err_limit_percent_z_jonswap   = 4.40490f, // measured 4.38297987
    .err_limit_percent_z_pmstokes  = 4.35577f, // measured 4.33409595
    .err_limit_yaw_deg             = 0.95465f, // measured 0.949897051
    .err_limit_roll_deg            = 0.35239f, // measured 0.350634784
    .err_limit_pitch_deg           = 0.09390f, // measured 0.0934322327
    .err_limit_percent_3d_jonswap  = 9.82421f, // measured 9.77532673
    .err_limit_percent_3d_pmstokes = 9.92000f, // measured 9.88502121
    .acc_z_bias_percent            = 4.34785f, // measured 4.32621
    .bias_3d_percent               = 68.35106f, // measured 68.011
    .gyro_bias_3d_percent          = 10.40186f, // measured 10.3501
};

inline constexpr W3dFailureLimits carma{
    .err_limit_percent_z_jonswap   = 3.07238f, // measured 3.05709362
    .err_limit_percent_z_pmstokes  = 3.23081f, // measured 3.21473551
    .err_limit_yaw_deg             = 0.76470f, // measured 0.760889888
    .err_limit_roll_deg            = 0.35048f, // measured 0.348729193
    .err_limit_pitch_deg           = 0.20325f, // measured 0.202237442
    .err_limit_percent_3d_jonswap  = 4.86855f, // measured 4.84432125
    .err_limit_percent_3d_pmstokes = 5.25706f, // measured 5.23089981
    .acc_z_bias_percent            = 4.06168f, // measured 4.04147
    .bias_3d_percent               = 81.73414f, // measured 81.3275
    .gyro_bias_3d_percent          = 12.78000f, // measured 12.7244
};

inline constexpr W3dFailureLimits csm{
    .err_limit_percent_z_jonswap   = 2.46700f, // measured 2.45491838
    .err_limit_percent_z_pmstokes  = 2.41541f, // measured 2.40339065
    .err_limit_yaw_deg             = 0.75808f, // measured 0.754300177
    .err_limit_roll_deg            = 0.35291f, // measured 0.351149976
    .err_limit_pitch_deg           = 0.21288f, // measured 0.211819157
    .err_limit_percent_3d_jonswap  = 4.41992f, // measured 4.39792442
    .err_limit_percent_3d_pmstokes = 4.88513f, // measured 4.86081743
    .acc_z_bias_percent            = 4.03272f, // measured 4.01265
    .bias_3d_percent               = 83.27692f, // measured 82.8626
    .gyro_bias_3d_percent          = 12.71000f, // measured 12.6637
};

inline constexpr W3dFailureLimits tcarma{
    .err_limit_percent_z_jonswap   = 3.07200f, // measured 3.05673885
    .err_limit_percent_z_pmstokes  = 3.23061f, // measured 3.21453714
    .err_limit_yaw_deg             = 0.77351f, // measured 0.769657433
    .err_limit_roll_deg            = 0.35073f, // measured 0.348984271
    .err_limit_pitch_deg           = 0.20364f, // measured 0.202624053
    .err_limit_percent_3d_jonswap  = 4.87124f, // measured 4.84699917
    .err_limit_percent_3d_pmstokes = 5.25855f, // measured 5.23238039
    .acc_z_bias_percent            = 4.06068f, // measured 4.04047
    .bias_3d_percent               = 81.82580f, // measured 81.4187
    .gyro_bias_3d_percent          = 12.40000f, // measured 12.3488
};

inline constexpr W3dFailureLimits tcsm{
    .err_limit_percent_z_jonswap   = 2.46690f, // measured 2.45477271
    .err_limit_percent_z_pmstokes  = 2.41520f, // measured 2.40329385
    .err_limit_yaw_deg             = 0.76168f, // measured 0.757887661
    .err_limit_roll_deg            = 0.35224f, // measured 0.350478113
    .err_limit_pitch_deg           = 0.20862f, // measured 0.207578912
    .err_limit_percent_3d_jonswap  = 4.42280f, // measured 4.40079069
    .err_limit_percent_3d_pmstokes = 4.88809f, // measured 4.86376667
    .acc_z_bias_percent            = 4.02852f, // measured 4.00847
    .bias_3d_percent               = 82.68527f, // measured 82.2739
    .gyro_bias_3d_percent          = 12.34600f, // measured 12.3014
};

inline constexpr W3dFailureLimits tfg_legacy{
    .err_limit_percent_z_jonswap   = 4.34379f, // measured 4.32217884
    .err_limit_percent_z_pmstokes  = 4.31684f, // measured 4.29535675
    .err_limit_yaw_deg             = 1.06248f, // measured 1.05719304
    .err_limit_roll_deg            = 0.43509f, // measured 0.432924926
    .err_limit_pitch_deg           = 0.49514f, // measured 0.492675722
    .err_limit_percent_3d_jonswap  = 18.65808f, // measured 18.5652485
    .err_limit_percent_3d_pmstokes = 19.03321f, // measured 18.938509
    .acc_z_bias_percent            = 4.30249f, // measured 4.28108
    .bias_3d_percent               = 120.86733f, // measured 120.266
    .gyro_bias_3d_percent          = 20.99466f, // measured 20.8902
};

} // namespace ocean_imu::sim::current_limits
