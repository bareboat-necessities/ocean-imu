#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// v1.3.1 28 ft, 5 kn, +/-3 deg heading / 120 s; exact 200 Hz clocks.
// Default seeds; trailing 900 s. Each profile: worst + 0.5%, rounded up
// to 1e-5. Evidence and reproduction: docs/filter-simulator-gates.md.
#include "util/W3dSimCommon.h"
namespace ocean_imu::sim::current_limits {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 20.94565f,
    .err_limit_percent_z_pmstokes  = 11.37117f,
    .err_limit_yaw_deg             = 0.89796f,
    .err_limit_roll_deg            = 0.58935f,
    .err_limit_pitch_deg           = 0.50524f,
    .err_limit_percent_3d_jonswap  = 31.14368f,
    .err_limit_percent_3d_pmstokes = 19.12513f,
    .acc_z_bias_percent            = 7.32224f,
    .bias_3d_percent               = 157.12572f,
    .gyro_bias_3d_percent          = 12.56401f,
};

inline constexpr W3dFailureLimits ou3{
    .err_limit_percent_z_jonswap   = 8.27723f,
    .err_limit_percent_z_pmstokes  = 7.87831f,
    .err_limit_yaw_deg             = 0.87554f,
    .err_limit_roll_deg            = 0.34667f,
    .err_limit_pitch_deg           = 0.15584f,
    .err_limit_percent_3d_jonswap  = 19.43353f,
    .err_limit_percent_3d_pmstokes = 19.15242f,
    .acc_z_bias_percent            = 4.45285f,
    .bias_3d_percent               = 75.64113f,
    .gyro_bias_3d_percent          = 12.69245f,
};

inline constexpr W3dFailureLimits tfg{
    .err_limit_percent_z_jonswap   = 7.85835f,
    .err_limit_percent_z_pmstokes  = 7.51301f,
    .err_limit_yaw_deg             = 1.64146f,
    .err_limit_roll_deg            = 0.67022f,
    .err_limit_pitch_deg           = 0.65965f,
    .err_limit_percent_3d_jonswap  = 17.75326f,
    .err_limit_percent_3d_pmstokes = 17.33916f,
    .acc_z_bias_percent            = 5.01148f,
    .bias_3d_percent               = 186.73905f,
    .gyro_bias_3d_percent          = 10.99611f,
};

inline constexpr W3dFailureLimits carma{
    .err_limit_percent_z_jonswap   = 4.21295f,
    .err_limit_percent_z_pmstokes  = 4.69083f,
    .err_limit_yaw_deg             = 0.94167f,
    .err_limit_roll_deg            = 0.31108f,
    .err_limit_pitch_deg           = 0.19556f,
    .err_limit_percent_3d_jonswap  = 6.89035f,
    .err_limit_percent_3d_pmstokes = 8.07388f,
    .acc_z_bias_percent            = 4.32003f,
    .bias_3d_percent               = 72.10192f,
    .gyro_bias_3d_percent          = 12.41417f,
};

inline constexpr W3dFailureLimits csm{
    .err_limit_percent_z_jonswap   = 3.56798f,
    .err_limit_percent_z_pmstokes  = 3.72017f,
    .err_limit_yaw_deg             = 0.99056f,
    .err_limit_roll_deg            = 0.21343f,
    .err_limit_pitch_deg           = 0.21032f,
    .err_limit_percent_3d_jonswap  = 6.11027f,
    .err_limit_percent_3d_pmstokes = 6.60492f,
    .acc_z_bias_percent            = 3.94444f,
    .bias_3d_percent               = 56.55326f,
    .gyro_bias_3d_percent          = 12.38804f,
};

inline constexpr W3dFailureLimits tcarma{
    .err_limit_percent_z_jonswap   = 4.21316f,
    .err_limit_percent_z_pmstokes  = 4.69209f,
    .err_limit_yaw_deg             = 0.93962f,
    .err_limit_roll_deg            = 0.31016f,
    .err_limit_pitch_deg           = 0.19382f,
    .err_limit_percent_3d_jonswap  = 6.89010f,
    .err_limit_percent_3d_pmstokes = 8.07473f,
    .acc_z_bias_percent            = 4.31899f,
    .bias_3d_percent               = 72.29649f,
    .gyro_bias_3d_percent          = 12.60914f,
};

inline constexpr W3dFailureLimits tcsm{
    .err_limit_percent_z_jonswap   = 3.56832f,
    .err_limit_percent_z_pmstokes  = 3.72067f,
    .err_limit_yaw_deg             = 0.95680f,
    .err_limit_roll_deg            = 0.21302f,
    .err_limit_pitch_deg           = 0.20167f,
    .err_limit_percent_3d_jonswap  = 6.11056f,
    .err_limit_percent_3d_pmstokes = 6.60548f,
    .acc_z_bias_percent            = 3.94022f,
    .bias_3d_percent               = 56.70200f,
    .gyro_bias_3d_percent          = 12.52693f,
};

inline constexpr W3dFailureLimits tfg_legacy{
    .err_limit_percent_z_jonswap   = 8.28795f,
    .err_limit_percent_z_pmstokes  = 8.05605f,
    .err_limit_yaw_deg             = 3.06465f,
    .err_limit_roll_deg            = 0.40280f,
    .err_limit_pitch_deg           = 1.16598f,
    .err_limit_percent_3d_jonswap  = 31.15484f,
    .err_limit_percent_3d_pmstokes = 30.90172f,
    .acc_z_bias_percent            = 5.51543f,
    .bias_3d_percent               = 250.68519f,
    .gyro_bias_3d_percent          = 27.59379f,
};

namespace interval {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 20.94551f,
    .err_limit_percent_z_pmstokes  = 11.37150f,
    .err_limit_yaw_deg             = 0.89936f,
    .err_limit_roll_deg            = 0.58666f,
    .err_limit_pitch_deg           = 0.50723f,
    .err_limit_percent_3d_jonswap  = 31.16022f,
    .err_limit_percent_3d_pmstokes = 19.15704f,
    .acc_z_bias_percent            = 7.32191f,
    .bias_3d_percent               = 156.98100f,
    .gyro_bias_3d_percent          = 12.55627f,
};

inline constexpr W3dFailureLimits ou3{
    .err_limit_percent_z_jonswap   = 8.27822f,
    .err_limit_percent_z_pmstokes  = 7.87925f,
    .err_limit_yaw_deg             = 0.86774f,
    .err_limit_roll_deg            = 0.34562f,
    .err_limit_pitch_deg           = 0.15798f,
    .err_limit_percent_3d_jonswap  = 19.35139f,
    .err_limit_percent_3d_pmstokes = 19.08362f,
    .acc_z_bias_percent            = 4.48809f,
    .bias_3d_percent               = 75.55932f,
    .gyro_bias_3d_percent          = 12.68089f,
};

inline constexpr W3dFailureLimits tfg{
    .err_limit_percent_z_jonswap   = 7.85941f,
    .err_limit_percent_z_pmstokes  = 7.51397f,
    .err_limit_yaw_deg             = 1.62600f,
    .err_limit_roll_deg            = 0.66877f,
    .err_limit_pitch_deg           = 0.65359f,
    .err_limit_percent_3d_jonswap  = 17.69845f,
    .err_limit_percent_3d_pmstokes = 17.28851f,
    .acc_z_bias_percent            = 5.13306f,
    .bias_3d_percent               = 185.60441f,
    .gyro_bias_3d_percent          = 10.99079f,
};

inline constexpr W3dFailureLimits carma{
    .err_limit_percent_z_jonswap   = 4.21299f,
    .err_limit_percent_z_pmstokes  = 4.69123f,
    .err_limit_yaw_deg             = 0.94383f,
    .err_limit_roll_deg            = 0.31009f,
    .err_limit_pitch_deg           = 0.19758f,
    .err_limit_percent_3d_jonswap  = 6.93610f,
    .err_limit_percent_3d_pmstokes = 8.12158f,
    .acc_z_bias_percent            = 4.35243f,
    .bias_3d_percent               = 72.03609f,
    .gyro_bias_3d_percent          = 12.40271f,
};

inline constexpr W3dFailureLimits csm{
    .err_limit_percent_z_jonswap   = 3.56817f,
    .err_limit_percent_z_pmstokes  = 3.72045f,
    .err_limit_yaw_deg             = 0.99092f,
    .err_limit_roll_deg            = 0.21307f,
    .err_limit_pitch_deg           = 0.21094f,
    .err_limit_percent_3d_jonswap  = 6.14746f,
    .err_limit_percent_3d_pmstokes = 6.64710f,
    .acc_z_bias_percent            = 3.98894f,
    .bias_3d_percent               = 56.53819f,
    .gyro_bias_3d_percent          = 12.37638f,
};

inline constexpr W3dFailureLimits tcarma{
    .err_limit_percent_z_jonswap   = 4.21322f,
    .err_limit_percent_z_pmstokes  = 4.69250f,
    .err_limit_yaw_deg             = 0.94157f,
    .err_limit_roll_deg            = 0.30920f,
    .err_limit_pitch_deg           = 0.19587f,
    .err_limit_percent_3d_jonswap  = 6.93581f,
    .err_limit_percent_3d_pmstokes = 8.12238f,
    .acc_z_bias_percent            = 4.35075f,
    .bias_3d_percent               = 72.23870f,
    .gyro_bias_3d_percent          = 12.59627f,
};

inline constexpr W3dFailureLimits tcsm{
    .err_limit_percent_z_jonswap   = 3.56850f,
    .err_limit_percent_z_pmstokes  = 3.72095f,
    .err_limit_yaw_deg             = 0.95754f,
    .err_limit_roll_deg            = 0.21266f,
    .err_limit_pitch_deg           = 0.20247f,
    .err_limit_percent_3d_jonswap  = 6.14765f,
    .err_limit_percent_3d_pmstokes = 6.64760f,
    .acc_z_bias_percent            = 3.98335f,
    .bias_3d_percent               = 56.69075f,
    .gyro_bias_3d_percent          = 12.51638f,
};
} // namespace interval
} // namespace ocean_imu::sim::current_limits
