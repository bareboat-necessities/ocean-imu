#!/usr/bin/env python3
"""Patch the Actions checkout with the derived separate TFG X/Y R_S law.

Experimental PR only. The production header on main is not changed. The
registered tfg-validation workflow builds this checkout and an untouched exact
base checkout, so its two simulator logs are a paired control/experiment.
"""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
path = root / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"
text = path.read_text()
old = "        mekf_.set_RS_noise(Vector3f(rs * R_S_xy_factor_, rs * R_S_xy_factor_, rs));"
new = r'''        // EXPERIMENT ONLY: separate TFG X/Y regularization from the
        // reduced directional 3-D MSE model. NED X is generator North and
        // NED Y is generator East, hence the swapped measured wave ratios.
        constexpr float sigma_x_over_z = 0.569f;
        constexpr float sigma_y_over_z = 0.819f;
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f)
                                  ? sigma_coeff_ : 1.0f;
        const float sigma_z = std::max(tune_.sigma_applied / c_sigma, 1.0e-6f);
        const float qz = std::max(2.0f * rs_accel_noise_density_, 1.0e-12f);

        // X/pitch channel: magnetometer-observed tilt contribution. Use the
        // filter's configured magnetic measurement standard deviation and the
        // simulator's 25 Hz magnetic cadence; use the learned field norm once
        // it exists, otherwise the simulator's 52 uT WMM magnitude.
        const float mag_std_uT = std::max(cfg_.sigma_m.maxCoeff(), 1.0e-6f);
        constexpr float mag_dt_sec = 1.0f / 25.0f;
        const float B_uT = mekf_.has_magnetic_reference()
            ? std::max(mekf_.magnetic_reference_world().norm(), 1.0e-3f)
            : 52.0f;
        const float g = cfg_.gravity_magnitude;
        const float qx = qz + 2.0f * g * g *
            (mag_std_uT * mag_std_uT * mag_dt_sec) / (B_uT * B_uT);

        // Y/roll channel: reduced low-frequency roll/OU ambiguity model.
        // qz is a floor so the horizontal model never claims less residual
        // low-frequency acceleration error than the directly observed Z axis.
        const float sigma_y = sigma_y_over_z * sigma_z;
        const float qy = std::max(qz,
            4.0f * sigma_y * sigma_y * std::max(tune_.tau_applied, 1.0e-3f));

        const float kx = std::pow(sigma_x_over_z, 6.0f / 7.0f) *
                         std::pow(qx / qz, 1.0f / 14.0f);
        const float ky = std::pow(sigma_y_over_z, 6.0f / 7.0f) *
                         std::pow(qy / qz, 1.0f / 14.0f);
        mekf_.set_RS_noise(Vector3f(rs * kx, rs * ky, rs));'''

count = text.count(old)
if count == 1:
    path.write_text(text.replace(old, new))
elif count == 0 and "constexpr float sigma_x_over_z = 0.569f;" in text:
    pass
else:
    raise SystemExit(f"expected one production R_S assignment, found {count}")

# The production executable intentionally exits at the first failed sentinel.
# This experiment needs all eight records, so force the existing collect-all
# mode in the experimental simulator only. The workflow already uses `|| true`.
sim = root / "tests/kalman_tfg/kalman_tfg-sim.cpp"
sim_text = sim.read_text()
needle = "int main(int argc, char* argv[]) {\n"
insert = (
    "int main(int argc, char* argv[]) {\n"
    "    // EXPERIMENT ONLY: report every record even when a production gate fails.\n"
    "    setenv(\"W3D_COLLECT_ALL_GATES\", \"1\", 1);\n"
)
if insert not in sim_text:
    if sim_text.count(needle) != 1:
        raise SystemExit("simulator main() insertion point not found")
    sim.write_text(sim_text.replace(needle, insert))
