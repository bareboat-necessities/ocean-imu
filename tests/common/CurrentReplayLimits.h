#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// v1.3.1 28 ft, 5 kn, +/-3 deg heading / 120 s; exact 200 Hz clocks.
// Default seeds; trailing 900 s. Each profile: worst + 0.5%, rounded up
// to 1e-5. Evidence and reproduction: docs/filter-simulator-gates.md.
#include "util/W3dSimCommon.h"
namespace ocean_imu::sim::current_limits {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 8.52521f,
    .err_limit_percent_z_pmstokes  = 8.06839f,
    .err_limit_yaw_deg             = 0.90353f,
    .err_limit_roll_deg            = 0.21405f,
    .err_limit_pitch_deg           = 0.18845f,
    .err_limit_percent_3d_jonswap  = 16.40907f,
    .err_limit_percent_3d_pmstokes = 16.08768f,
    .acc_z_bias_percent            = 4.60839f,
    .bias_3d_percent               = 53.91886f,
    .gyro_bias_3d_percent          = 12.57034f,
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
    .err_limit_percent_z_jonswap   = 3.58976f,
    .err_limit_percent_z_pmstokes  = 3.74505f,
    .err_limit_yaw_deg             = 0.94337f,
    .err_limit_roll_deg            = 0.21201f,
    .err_limit_pitch_deg           = 0.19463f,
    .err_limit_percent_3d_jonswap  = 6.41330f,
    .err_limit_percent_3d_pmstokes = 6.96590f,
    .acc_z_bias_percent            = 3.95950f,
    .bias_3d_percent               = 56.62221f,
    .gyro_bias_3d_percent          = 12.43879f,
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
    .err_limit_percent_z_jonswap   = 3.58978f,
    .err_limit_percent_z_pmstokes  = 3.74520f,
    .err_limit_yaw_deg             = 0.93217f,
    .err_limit_roll_deg            = 0.21205f,
    .err_limit_pitch_deg           = 0.19230f,
    .err_limit_percent_3d_jonswap  = 6.41332f,
    .err_limit_percent_3d_pmstokes = 6.96604f,
    .acc_z_bias_percent            = 3.95801f,
    .bias_3d_percent               = 56.69095f,
    .gyro_bias_3d_percent          = 12.57688f,
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
    .err_limit_percent_z_jonswap   = 8.52594f,
    .err_limit_percent_z_pmstokes  = 8.06917f,
    .err_limit_yaw_deg             = 0.90485f,
    .err_limit_roll_deg            = 0.21381f,
    .err_limit_pitch_deg           = 0.18872f,
    .err_limit_percent_3d_jonswap  = 16.45571f,
    .err_limit_percent_3d_pmstokes = 16.13084f,
    .acc_z_bias_percent            = 4.63046f,
    .bias_3d_percent               = 53.89946f,
    .gyro_bias_3d_percent          = 12.56140f,
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
    .err_limit_percent_z_jonswap   = 3.58991f,
    .err_limit_percent_z_pmstokes  = 3.74527f,
    .err_limit_yaw_deg             = 0.94446f,
    .err_limit_roll_deg            = 0.21179f,
    .err_limit_pitch_deg           = 0.19655f,
    .err_limit_percent_3d_jonswap  = 6.45718f,
    .err_limit_percent_3d_pmstokes = 7.01268f,
    .acc_z_bias_percent            = 4.00342f,
    .bias_3d_percent               = 56.62572f,
    .gyro_bias_3d_percent          = 12.42613f,
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
    .err_limit_percent_z_jonswap   = 3.58993f,
    .err_limit_percent_z_pmstokes  = 3.74543f,
    .err_limit_yaw_deg             = 0.93332f,
    .err_limit_roll_deg            = 0.21184f,
    .err_limit_pitch_deg           = 0.19404f,
    .err_limit_percent_3d_jonswap  = 6.45713f,
    .err_limit_percent_3d_pmstokes = 7.01274f,
    .acc_z_bias_percent            = 4.00092f,
    .bias_3d_percent               = 56.69618f,
    .gyro_bias_3d_percent          = 12.56572f,
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
