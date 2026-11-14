#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Production-default regression limits measured on main a053a1d2 (v2.7.6);
// OU-II, CARMA and TFG-CARMA measured without stillness variance attenuation
// (docs/moving-vessel-regression-diagnosis.md).
// Eight pinned v1.2.3 vessel-RAO records, default seeds, trailing 900 s.
// Worst score + 0.5%, rounded upward to 1e-5, capped by an existing
// stricter production limit. Historical --regression bars stay separate.
// Evidence and reproduction: docs/filter-simulator-gates.md.
#include "util/W3dSimCommon.h"
namespace ocean_imu::sim::current_limits {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 6.67400f, // measured 6.64418745
    .err_limit_percent_z_pmstokes  = 6.47142f, // measured 6.43922281
    .err_limit_yaw_deg             = 0.74609f, // measured 0.742370069
    .err_limit_roll_deg            = 0.22529f, // measured 0.224164888
    .err_limit_pitch_deg           = 0.21279f, // measured 0.211723939
    .err_limit_percent_3d_jonswap  = 12.20210f, // measured 12.1413851
    .err_limit_percent_3d_pmstokes = 12.03067f, // measured 11.970809
    .acc_z_bias_percent            = 4.58776f, // measured 4.56493
    .bias_3d_percent               = 62.00278f, // measured 61.6943
    .gyro_bias_3d_percent          = 12.64763f, // measured 12.5847
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
    .err_limit_percent_z_jonswap   = 2.47299f, // measured 2.4606843
    .err_limit_percent_z_pmstokes  = 2.45054f, // measured 2.43834066
    .err_limit_yaw_deg             = 0.76470f, // measured 0.763234377
    .err_limit_roll_deg            = 0.35093f, // measured 0.349181592
    .err_limit_pitch_deg           = 0.20365f, // measured 0.202634886
    .err_limit_percent_3d_jonswap  = 4.88800f, // measured 4.88411903
    .err_limit_percent_3d_pmstokes = 5.27600f, // measured 5.25382423
    .acc_z_bias_percent            = 4.04775f, // measured 4.02761
    .bias_3d_percent               = 81.84731f, // measured 81.4401
    .gyro_bias_3d_percent          = 12.78000f, // measured 12.7238
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
    .err_limit_percent_z_jonswap   = 2.47298f, // measured 2.46067643
    .err_limit_percent_z_pmstokes  = 2.45038f, // measured 2.43818688
    .err_limit_yaw_deg             = 0.77351f, // measured 0.770234048
    .err_limit_roll_deg            = 0.35118f, // measured 0.349422902
    .err_limit_pitch_deg           = 0.20400f, // measured 0.202977106
    .err_limit_percent_3d_jonswap  = 4.89000f, // measured 4.88696814
    .err_limit_percent_3d_pmstokes = 5.27800f, // measured 5.2571454
    .acc_z_bias_percent            = 4.04668f, // measured 4.02654
    .bias_3d_percent               = 81.93012f, // measured 81.5225
    .gyro_bias_3d_percent          = 12.40000f, // measured 12.3564
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
