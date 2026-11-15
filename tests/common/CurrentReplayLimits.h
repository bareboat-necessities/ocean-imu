#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// v1.3.1 28 ft, 5 kn, +/-3 deg heading / 120 s; exact 200 Hz clocks.
// Default seeds; trailing 900 s. Worst + 0.5%, rounded up to 1e-5,
// capped by the previous limit: tightening never widens a gate.
// Evidence and reproduction: docs/filter-simulator-gates.md.
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
    .err_limit_percent_z_jonswap   = 5.74481f,
    .err_limit_percent_z_pmstokes  = 5.48631f,
    .err_limit_yaw_deg             = 1.07490f,
    .err_limit_roll_deg            = 0.53127f,
    .err_limit_pitch_deg           = 0.24688f,
    .err_limit_percent_3d_jonswap  = 16.05850f,
    .err_limit_percent_3d_pmstokes = 15.63066f,
    .acc_z_bias_percent            = 4.33085f,
    .bias_3d_percent               = 109.39526f,
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
    .err_limit_percent_z_jonswap   = 3.47282f,
    .err_limit_percent_z_pmstokes  = 3.54933f,
    .err_limit_yaw_deg             = 0.98215f,
    .err_limit_roll_deg            = 0.20289f,
    .err_limit_pitch_deg           = 0.20729f,
    .err_limit_percent_3d_jonswap  = 5.98179f,
    .err_limit_percent_3d_pmstokes = 6.38260f,
    .acc_z_bias_percent            = 3.94444f,
    .bias_3d_percent               = 54.86607f,
    .gyro_bias_3d_percent          = 12.31115f,
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
    .err_limit_percent_z_jonswap   = 3.47367f,
    .err_limit_percent_z_pmstokes  = 3.55102f,
    .err_limit_yaw_deg             = 0.95173f,
    .err_limit_roll_deg            = 0.19997f,
    .err_limit_pitch_deg           = 0.20021f,
    .err_limit_percent_3d_jonswap  = 5.98221f,
    .err_limit_percent_3d_pmstokes = 6.38402f,
    .acc_z_bias_percent            = 3.94022f,
    .bias_3d_percent               = 54.32046f,
    .gyro_bias_3d_percent          = 12.45045f,
};

inline constexpr W3dFailureLimits tfg_legacy{
    .err_limit_percent_z_jonswap   = 4.71275f,
    .err_limit_percent_z_pmstokes  = 4.85658f,
    .err_limit_yaw_deg             = 3.06465f,
    .err_limit_roll_deg            = 0.40280f,
    .err_limit_pitch_deg           = 1.16598f,
    .err_limit_percent_3d_jonswap  = 29.62384f,
    .err_limit_percent_3d_pmstokes = 29.43656f,
    .acc_z_bias_percent            = 5.05118f,
    .bias_3d_percent               = 250.55153f,
    .gyro_bias_3d_percent          = 27.49510f,
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
    .err_limit_percent_z_jonswap   = 5.74575f,
    .err_limit_percent_z_pmstokes  = 5.48726f,
    .err_limit_yaw_deg             = 1.06609f,
    .err_limit_roll_deg            = 0.53090f,
    .err_limit_pitch_deg           = 0.24785f,
    .err_limit_percent_3d_jonswap  = 15.99806f,
    .err_limit_percent_3d_pmstokes = 15.57394f,
    .acc_z_bias_percent            = 4.40964f,
    .bias_3d_percent               = 108.88371f,
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
    .err_limit_percent_z_jonswap   = 3.47315f,
    .err_limit_percent_z_pmstokes  = 3.54971f,
    .err_limit_yaw_deg             = 0.98251f,
    .err_limit_roll_deg            = 0.20238f,
    .err_limit_pitch_deg           = 0.20790f,
    .err_limit_percent_3d_jonswap  = 6.01909f,
    .err_limit_percent_3d_pmstokes = 6.42610f,
    .acc_z_bias_percent            = 3.98894f,
    .bias_3d_percent               = 54.93562f,
    .gyro_bias_3d_percent          = 12.30060f,
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
    .err_limit_percent_z_jonswap   = 3.47397f,
    .err_limit_percent_z_pmstokes  = 3.55136f,
    .err_limit_yaw_deg             = 0.95246f,
    .err_limit_roll_deg            = 0.19958f,
    .err_limit_pitch_deg           = 0.20097f,
    .err_limit_percent_3d_jonswap  = 6.01924f,
    .err_limit_percent_3d_pmstokes = 6.42732f,
    .acc_z_bias_percent            = 3.98335f,
    .bias_3d_percent               = 54.31051f,
    .gyro_bias_3d_percent          = 12.43839f,
};
} // namespace interval
} // namespace ocean_imu::sim::current_limits
