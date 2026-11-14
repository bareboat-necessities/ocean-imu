#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// v1.3.0 28 ft, 5 kn, +/-3 deg heading / 120 s; eight pinned records.
// Default seeds; trailing 900 s. Each profile: worst + 0.5%, rounded up
// to 1e-5. This changed physical experiment has its own regression baseline.
// Evidence and reproduction: docs/filter-simulator-gates.md.
#include "util/W3dSimCommon.h"
namespace ocean_imu::sim::current_limits {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 20.96609f,
    .err_limit_percent_z_pmstokes  = 11.36901f,
    .err_limit_yaw_deg             = 0.89815f,
    .err_limit_roll_deg            = 0.58950f,
    .err_limit_pitch_deg           = 0.50556f,
    .err_limit_percent_3d_jonswap  = 31.23156f,
    .err_limit_percent_3d_pmstokes = 19.23342f,
    .acc_z_bias_percent            = 7.32554f,
    .bias_3d_percent               = 157.19004f,
    .gyro_bias_3d_percent          = 12.49557f,
};

inline constexpr W3dFailureLimits ou3{
    .err_limit_percent_z_jonswap   = 8.27699f,
    .err_limit_percent_z_pmstokes  = 7.87839f,
    .err_limit_yaw_deg             = 0.87291f,
    .err_limit_roll_deg            = 0.34668f,
    .err_limit_pitch_deg           = 0.15585f,
    .err_limit_percent_3d_jonswap  = 19.51916f,
    .err_limit_percent_3d_pmstokes = 19.24972f,
    .acc_z_bias_percent            = 4.45288f,
    .bias_3d_percent               = 75.64183f,
    .gyro_bias_3d_percent          = 12.69014f,
};

inline constexpr W3dFailureLimits tfg{
    .err_limit_percent_z_jonswap   = 7.85818f,
    .err_limit_percent_z_pmstokes  = 7.51302f,
    .err_limit_yaw_deg             = 1.64249f,
    .err_limit_roll_deg            = 0.67032f,
    .err_limit_pitch_deg           = 0.66006f,
    .err_limit_percent_3d_jonswap  = 17.81642f,
    .err_limit_percent_3d_pmstokes = 17.42364f,
    .acc_z_bias_percent            = 5.01222f,
    .bias_3d_percent               = 186.81141f,
    .gyro_bias_3d_percent          = 10.86828f,
};

inline constexpr W3dFailureLimits carma{
    .err_limit_percent_z_jonswap   = 4.21293f,
    .err_limit_percent_z_pmstokes  = 4.69086f,
    .err_limit_yaw_deg             = 0.94186f,
    .err_limit_roll_deg            = 0.31108f,
    .err_limit_pitch_deg           = 0.19557f,
    .err_limit_percent_3d_jonswap  = 7.23442f,
    .err_limit_percent_3d_pmstokes = 8.39239f,
    .acc_z_bias_percent            = 4.32012f,
    .bias_3d_percent               = 72.10242f,
    .gyro_bias_3d_percent          = 12.47085f,
};

inline constexpr W3dFailureLimits csm{
    .err_limit_percent_z_jonswap   = 3.56795f,
    .err_limit_percent_z_pmstokes  = 3.72019f,
    .err_limit_yaw_deg             = 0.99081f,
    .err_limit_roll_deg            = 0.21343f,
    .err_limit_pitch_deg           = 0.21034f,
    .err_limit_percent_3d_jonswap  = 6.45326f,
    .err_limit_percent_3d_pmstokes = 6.95017f,
    .acc_z_bias_percent            = 3.94449f,
    .bias_3d_percent               = 56.55357f,
    .gyro_bias_3d_percent          = 12.42713f,
};

inline constexpr W3dFailureLimits tcarma{
    .err_limit_percent_z_jonswap   = 4.21313f,
    .err_limit_percent_z_pmstokes  = 4.69212f,
    .err_limit_yaw_deg             = 0.93966f,
    .err_limit_roll_deg            = 0.31016f,
    .err_limit_pitch_deg           = 0.19379f,
    .err_limit_percent_3d_jonswap  = 7.23423f,
    .err_limit_percent_3d_pmstokes = 8.39326f,
    .acc_z_bias_percent            = 4.31900f,
    .bias_3d_percent               = 72.29629f,
    .gyro_bias_3d_percent          = 12.68109f,
};

inline constexpr W3dFailureLimits tcsm{
    .err_limit_percent_z_jonswap   = 3.56830f,
    .err_limit_percent_z_pmstokes  = 3.72069f,
    .err_limit_yaw_deg             = 0.95692f,
    .err_limit_roll_deg            = 0.21302f,
    .err_limit_pitch_deg           = 0.20170f,
    .err_limit_percent_3d_jonswap  = 6.45375f,
    .err_limit_percent_3d_pmstokes = 6.95079f,
    .acc_z_bias_percent            = 3.94028f,
    .bias_3d_percent               = 56.70190f,
    .gyro_bias_3d_percent          = 12.62180f,
};

inline constexpr W3dFailureLimits tfg_legacy{
    .err_limit_percent_z_jonswap   = 8.28776f,
    .err_limit_percent_z_pmstokes  = 8.05588f,
    .err_limit_yaw_deg             = 3.06507f,
    .err_limit_roll_deg            = 0.40277f,
    .err_limit_pitch_deg           = 1.16615f,
    .err_limit_percent_3d_jonswap  = 31.25140f,
    .err_limit_percent_3d_pmstokes = 31.01273f,
    .acc_z_bias_percent            = 5.51599f,
    .bias_3d_percent               = 250.71836f,
    .gyro_bias_3d_percent          = 27.67821f,
};

namespace interval {

inline constexpr W3dFailureLimits ou2{
    .err_limit_percent_z_jonswap   = 20.96593f,
    .err_limit_percent_z_pmstokes  = 11.36931f,
    .err_limit_yaw_deg             = 0.89959f,
    .err_limit_roll_deg            = 0.58681f,
    .err_limit_pitch_deg           = 0.50755f,
    .err_limit_percent_3d_jonswap  = 31.24807f,
    .err_limit_percent_3d_pmstokes = 19.26528f,
    .acc_z_bias_percent            = 7.32516f,
    .bias_3d_percent               = 157.04733f,
    .gyro_bias_3d_percent          = 12.48783f,
};

inline constexpr W3dFailureLimits ou3{
    .err_limit_percent_z_jonswap   = 8.27797f,
    .err_limit_percent_z_pmstokes  = 7.87933f,
    .err_limit_yaw_deg             = 0.86510f,
    .err_limit_roll_deg            = 0.34562f,
    .err_limit_pitch_deg           = 0.15798f,
    .err_limit_percent_3d_jonswap  = 19.43768f,
    .err_limit_percent_3d_pmstokes = 19.18152f,
    .acc_z_bias_percent            = 4.48806f,
    .bias_3d_percent               = 75.56033f,
    .gyro_bias_3d_percent          = 12.67607f,
};

inline constexpr W3dFailureLimits tfg{
    .err_limit_percent_z_jonswap   = 7.85920f,
    .err_limit_percent_z_pmstokes  = 7.51395f,
    .err_limit_yaw_deg             = 1.62863f,
    .err_limit_roll_deg            = 0.66907f,
    .err_limit_pitch_deg           = 0.65465f,
    .err_limit_percent_3d_jonswap  = 17.76207f,
    .err_limit_percent_3d_pmstokes = 17.37332f,
    .acc_z_bias_percent            = 5.13534f,
    .bias_3d_percent               = 185.80139f,
    .gyro_bias_3d_percent          = 10.86345f,
};

inline constexpr W3dFailureLimits carma{
    .err_limit_percent_z_jonswap   = 4.21296f,
    .err_limit_percent_z_pmstokes  = 4.69127f,
    .err_limit_yaw_deg             = 0.94402f,
    .err_limit_roll_deg            = 0.31009f,
    .err_limit_pitch_deg           = 0.19759f,
    .err_limit_percent_3d_jonswap  = 7.27803f,
    .err_limit_percent_3d_pmstokes = 8.43838f,
    .acc_z_bias_percent            = 4.35245f,
    .bias_3d_percent               = 72.03660f,
    .gyro_bias_3d_percent          = 12.45909f,
};

inline constexpr W3dFailureLimits csm{
    .err_limit_percent_z_jonswap   = 3.56814f,
    .err_limit_percent_z_pmstokes  = 3.72048f,
    .err_limit_yaw_deg             = 0.99115f,
    .err_limit_roll_deg            = 0.21308f,
    .err_limit_pitch_deg           = 0.21095f,
    .err_limit_percent_3d_jonswap  = 6.48880f,
    .err_limit_percent_3d_pmstokes = 6.99057f,
    .acc_z_bias_percent            = 3.98894f,
    .bias_3d_percent               = 56.53879f,
    .gyro_bias_3d_percent          = 12.41588f,
};

inline constexpr W3dFailureLimits tcarma{
    .err_limit_percent_z_jonswap   = 4.21319f,
    .err_limit_percent_z_pmstokes  = 4.69255f,
    .err_limit_yaw_deg             = 0.94162f,
    .err_limit_roll_deg            = 0.30920f,
    .err_limit_pitch_deg           = 0.19583f,
    .err_limit_percent_3d_jonswap  = 7.27780f,
    .err_limit_percent_3d_pmstokes = 8.43922f,
    .acc_z_bias_percent            = 4.35069f,
    .bias_3d_percent               = 72.24021f,
    .gyro_bias_3d_percent          = 12.66793f,
};

inline constexpr W3dFailureLimits tcsm{
    .err_limit_percent_z_jonswap   = 3.56849f,
    .err_limit_percent_z_pmstokes  = 3.72097f,
    .err_limit_yaw_deg             = 0.95776f,
    .err_limit_roll_deg            = 0.21266f,
    .err_limit_pitch_deg           = 0.20247f,
    .err_limit_percent_3d_jonswap  = 6.48918f,
    .err_limit_percent_3d_pmstokes = 6.99106f,
    .acc_z_bias_percent            = 3.98333f,
    .bias_3d_percent               = 56.69045f,
    .gyro_bias_3d_percent          = 12.60632f,
};
} // namespace interval
} // namespace ocean_imu::sim::current_limits
