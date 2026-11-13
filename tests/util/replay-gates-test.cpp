// Exercise the shared numeric gates with controlled estimation failures.
#include <limits>
#include <string>
#include "util/W3dSimCommon.h"

int main(int argc, char* argv[]) {
    if (argc != 2) return 2;
    const std::string mode = argv[1];
    W3dSimulationRunResult r;
    r.wave_type = WaveType::JONSWAP;
    r.wave_params.height = 1.0f;
    r.with_mag = true;
    constexpr size_t count = 901;
    for (auto* values : {&r.errs_x, &r.errs_y, &r.errs_z,
                         &r.errs_roll, &r.errs_pitch, &r.errs_yaw,
                         &r.ref_x, &r.ref_y, &r.ref_z,
                         &r.accb_err_x, &r.accb_err_y, &r.accb_err_z,
                         &r.gyrb_err_x, &r.gyrb_err_y, &r.gyrb_err_z,
                         &r.magb_err_x, &r.magb_err_y, &r.magb_err_z,
                         &r.accb_true_x, &r.accb_true_y, &r.accb_true_z,
                         &r.gyrb_true_x, &r.gyrb_true_y, &r.gyrb_true_z,
                         &r.magb_true_x, &r.magb_true_y, &r.magb_true_z}) {
        values->assign(count, 0.0f);
    }
    r.ref_z.assign(count, 0.5f);
    r.accb_true_z.assign(count, 0.1f);
    r.gyrb_true_z.assign(count, 0.001f);
    if (mode == "nan") r.errs_z.back() = std::numeric_limits<float>::quiet_NaN();
    else if (mode == "inf") r.accb_err_x.back() = std::numeric_limits<float>::infinity();
    else if (mode == "truth_nan") r.accb_true_x.back() = std::numeric_limits<float>::quiet_NaN();
    else if (mode == "missing") r.gyrb_err_x.clear();
    else if (mode == "acc_bias") r.accb_err_x.assign(count, 0.03f);
    else if (mode == "overflow") r.accb_err_x.assign(count, std::numeric_limits<float>::max());
    else if (mode == "gyro_bias") r.gyrb_err_x.assign(count, 0.00015f);
    else if (mode == "roll") r.errs_roll.assign(count, 0.2f);
    else if (mode == "quiet") {
        r.ref_z.assign(count, 0.0f);
        r.accb_true_z.assign(count, 0.0f);
        r.gyrb_true_z.assign(count, 0.0f);
    }
    else if (mode != "valid") return 2;
    const W3dFailureLimits limits{
        .err_limit_percent_z_jonswap = 1.0f,
        .err_limit_percent_z_pmstokes = 1.0f,
        .err_limit_yaw_deg = 0.1f,
        .err_limit_roll_deg = 0.1f,
        .err_limit_pitch_deg = 0.1f,
        .err_limit_percent_3d_jonswap = 1.0f,
        .err_limit_percent_3d_pmstokes = 1.0f,
        .acc_z_bias_percent = 20.0f,
        .bias_3d_percent = 20.0f,
        .gyro_bias_3d_percent = 10.0f,
    };
    print_summary_and_fail_if_needed(r, 1.0f, limits, {});
    return w3d_any_quality_gate_failed() ? 1 : 0;
}
