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

void test_default_is_isotropic_without_changing_rs_law() {
    Fusion f;
    Fusion::Config cfg;
    f.begin(cfg);

    check(f.getRSLaw() == Fusion::RSLaw::SpectralMSE,
          "axis-factor API changed the deployed adaptation law");
    check(near(f.getRSXFactor(), 1.0f), "default X R_S factor is not 1.0");
    check(near(f.getRSYFactor(), 1.0f), "default Y R_S factor is not 1.0");

    check(f.setFixedTuning(2.0f, 0.5f, 3.0f), "fixed tuning rejected");
    const Matrix3f R = integral_R(f);
    check(near(R(0,0), 9.0f), "default X R_S is not isotropic");
    check(near(R(1,1), 9.0f), "default Y R_S is not isotropic");
    check(near(R(2,2), 9.0f), "default Z R_S changed");
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

}  // namespace

int main() {
    test_default_is_isotropic_without_changing_rs_law();
    test_x_y_factors_are_independent();
    if (failures != 0) {
        std::cerr << failures << " TFG R_S axis-factor test(s) failed\n";
        return 1;
    }
    std::cout << "TFG independent R_S axis-factor tests passed\n";
    return 0;
}
