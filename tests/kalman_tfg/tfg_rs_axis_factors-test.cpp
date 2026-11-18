#define EIGEN_NON_ARDUINO

#include "kalman_tfg/SeaStateFusionFilter_TFG.h"

#include <cmath>
#include <iostream>

namespace {

using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
using Vector3f = Eigen::Vector3f;
using Matrix3f = Eigen::Matrix3f;

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near(float a, float b, float tol = 1.0e-5f) {
    return std::fabs(a - b) <= tol * std::max(1.0f, std::fabs(b));
}

Matrix3f integral_R(Fusion& f) {
    Vector3f r;
    Eigen::Matrix<float, 3, Fusion::Mekf::NX> H;
    Matrix3f R;
    f.mekf().integral_residual(r, H, R);
    return R;
}

void test_default_horizontal_factors_without_changing_rs_law() {
    Fusion f;
    Fusion::Config cfg;
    f.begin(cfg);

    check(f.getRSLaw() == Fusion::RSLaw::SpectralMSE,
          "axis-factor API changed the deployed adaptation law");
    check(near(f.getRSXFactor(), 1.08f), "default X R_S factor is not 1.08");
    check(near(f.getRSYFactor(), 1.15f), "default Y R_S factor is not 1.15");

    check(f.setFixedTuning(2.0f, 0.5f, 3.0f), "fixed tuning rejected");
    const Matrix3f R = integral_R(f);
    const float x_var = (3.0f * 1.08f) * (3.0f * 1.08f);
    const float y_var = (3.0f * 1.15f) * (3.0f * 1.15f);
    check(near(R(0,0), x_var), "default X R_S is not rho_x^2 times Z");
    check(near(R(1,1), y_var), "default Y R_S is not rho_y^2 times Z");
    check(near(R(2,2), 9.0f), "default Z R_S changed");
    f.mekf().set_acc_bias_updates_enabled(true);
    const Matrix3f after = integral_R(f);
    check(near(after(0,0), x_var*0.625f) && near(after(1,1), y_var*0.625f) && near(after(2,2), 9),
          "live spectral law did not select the matched horizontal observation");
    f.setEmbeddedFriendlyLegacyRSLaw(true);
    check(integral_R(f).isApprox(R), "legacy law changed its horizontal weighting");
    cfg.wave_matched_horizontal_integral = false;
    f.begin(cfg);
    f.setRSLaw(Fusion::RSLaw::SpectralMSE);
    f.setFixedTuning(2.0f, 0.5f, 3.0f);
    f.mekf().set_acc_bias_updates_enabled(true);
    check(integral_R(f).isApprox(R), "vertical-only configuration changed horizontal weighting");
}

void test_x_y_factors_are_independent() {
    Fusion f;
    f.begin(Fusion::Config{});
    f.setRSXFactor(1.25f);
    f.setRSYFactor(0.80f);
    check(near(f.getRSXFactor(), 1.25f), "X factor setter did not take effect");
    check(near(f.getRSYFactor(), 0.80f), "Y factor setter did not take effect");

    check(f.setFixedTuning(2.0f, 0.5f, 3.0f), "fixed tuning rejected");
    const Matrix3f R = integral_R(f);
    check(near(R(0,0), (3.0f * 1.25f) * (3.0f * 1.25f)),
          "X factor did not independently scale R_S");
    check(near(R(1,1), (3.0f * 0.80f) * (3.0f * 0.80f)),
          "Y factor did not independently scale R_S");
    check(near(R(2,2), 9.0f), "X/Y factors changed Z R_S");

    f.setRSXFactor(0.0f);
    f.setRSYFactor(-1.0f);
    check(near(f.getRSXFactor(), 1.25f), "invalid X factor was accepted");
    check(near(f.getRSYFactor(), 0.80f), "invalid Y factor was accepted");
}

void test_integral_uses_staged_period_and_tau_coefficient() {
    Fusion f;
    f.begin(Fusion::Config{});
    check(f.waveMatchedVerticalIntegral(), "wave-matched integral not enabled by default");
    check(f.setFixedTuning(2.0f, 0.5f, 3.0f), "fixed tuning rejected");
    const float beta = 4.0f / 9.8696044010893586188f;
    check(near(f.mekf().integral_velocity_weight(), beta), "wrong fixed-period weight");
    f.setTauCoeff(2.0f);
    check(near(f.mekf().integral_velocity_weight(), beta/4.0f),
          "OU scaling coefficient leaked into inferred wave period");
    check(f.mekf().integral_attitude_protected(), "spectral law lost its constrained gain");
    f.setEmbeddedFriendlyLegacyRSLaw(true);
    check(!f.mekf().integral_attitude_protected(), "legacy law did not retain its coupled gain");
    check(near(f.mekf().integral_velocity_weight(), beta/4.0f),
          "law change altered the wave period");
    f.setRSLaw(Fusion::RSLaw::SpectralMSE);
    check(f.mekf().integral_attitude_protected(), "law switch left a stale gain policy");
    const Matrix3f noise = integral_R(f);
    f.setWaveMatchedVerticalIntegral(false);
    check(f.mekf().integral_velocity_weight() == 0.0f, "ablation did not restore S=0");
    check((noise-integral_R(f)).norm() == 0.0f, "ablation also changed pseudo noise");
}

}  // namespace

int main() {
    test_default_horizontal_factors_without_changing_rs_law();
    test_x_y_factors_are_independent();
    test_integral_uses_staged_period_and_tau_coefficient();
    if (failures != 0) {
        std::cerr << failures << " TFG R_S axis-factor test(s) failed\n";
        return 1;
    }
    std::cout << "TFG independent R_S axis-factor tests passed\n";
    return 0;
}
