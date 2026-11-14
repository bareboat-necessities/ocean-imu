#pragma once
// Single unchanged fixture for TFG and tangent replay.
constexpr W3dFailureLimits kTfgSharedRegressionBars{
    .err_limit_percent_z_jonswap   = 2.4669f, // measured 2.4545393
    .err_limit_percent_z_pmstokes  = 2.4152f, // measured 2.4031157
    .err_limit_yaw_deg             = 0.7631f, // measured 0.7592051
    .err_limit_roll_deg            = 0.3585f, // measured 0.35668024
    .err_limit_pitch_deg           = 0.2188f, // measured 0.21763369
    .err_limit_percent_3d_jonswap  = 4.4354f, // measured 4.4132538
    .err_limit_percent_3d_pmstokes = 4.9068f, // measured 4.8822913
    .acc_z_bias_percent            = 4.1787f, // measured 4.1579
    .bias_3d_percent               = 83.4106f, // measured 82.9956
    .gyro_bias_3d_percent          = 12.3460f, // measured 12.2845
};
