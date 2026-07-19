#define EIGEN_NON_ARDUINO

#include "kalman_ou_common/KalmanOUCoreMath.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace detail = ocean_imu::kalman::ou_detail;

namespace {

using T = double;

bool check_close(T actual, T expected, T tolerance, const char* message) {
    if (std::abs(actual - expected) <= tolerance) return true;
    std::cerr << message << ": actual=" << actual
              << " expected=" << expected
              << " tolerance=" << tolerance << '\n';
    return false;
}

T phi_pa_exact(T x, T tau) {
    return tau * tau * (x + std::expm1(-x));
}

T phi_Sa_exact(T x, T tau) {
    return tau * tau * tau * (T(0.5) * x * x - x - std::expm1(-x));
}

bool test_phi_coefficients_at_branch_boundary() {
    constexpr T threshold = T(1e-2);
    constexpr T tau = T(2.3);
    constexpr T relative_offset = T(1e-9);

    const T x_below = threshold * (T(1) - relative_offset);
    const T x_at = threshold;
    const T x_above = threshold * (T(1) + relative_offset);

    const auto below = detail::safe_phi_A_coeffs(tau * x_below, tau);
    const auto at = detail::safe_phi_A_coeffs(tau * x_at, tau);
    const auto above = detail::safe_phi_A_coeffs(tau * x_above, tau);

    const T scale_pa = std::max(T(1), std::abs(phi_pa_exact(threshold, tau)));
    const T scale_Sa = std::max(T(1), std::abs(phi_Sa_exact(threshold, tau)));

    // Below the threshold the polynomial branch must approximate the exact
    // expressions through its retained order. At and above the threshold the
    // exact expm1 branch is used.
    if (!check_close(below.phi_pa, phi_pa_exact(x_below, tau),
                     T(1e-11) * scale_pa,
                     "phi_pa small-x branch mismatch")) return false;
    if (!check_close(below.phi_Sa, phi_Sa_exact(x_below, tau),
                     T(1e-13) * scale_Sa,
                     "phi_Sa small-x branch mismatch")) return false;
    if (!check_close(at.phi_pa, phi_pa_exact(x_at, tau),
                     T(1e-15) * scale_pa,
                     "phi_pa exact branch mismatch at threshold")) return false;
    if (!check_close(at.phi_Sa, phi_Sa_exact(x_at, tau),
                     T(1e-15) * scale_Sa,
                     "phi_Sa exact branch mismatch at threshold")) return false;
    if (!check_close(above.phi_pa, phi_pa_exact(x_above, tau),
                     T(1e-15) * scale_pa,
                     "phi_pa exact branch mismatch above threshold")) return false;
    if (!check_close(above.phi_Sa, phi_Sa_exact(x_above, tau),
                     T(1e-15) * scale_Sa,
                     "phi_Sa exact branch mismatch above threshold")) return false;

    // Guard against the manuscript typo that inserted a spurious factor 1/2
    // in tau^2 * (x + expm1(-x)).
    const T erroneous_half_factor = T(0.5) * phi_pa_exact(x_at, tau);
    if (std::abs(at.phi_pa - erroneous_half_factor) <= T(0.1) * std::abs(at.phi_pa)) {
        std::cerr << "phi_pa regressed to the erroneous half-factor formula\n";
        return false;
    }

    return true;
}

} // namespace

int main() {
    return test_phi_coefficients_at_branch_boundary() ? 0 : 1;
}
