// Verify the deployed OU-II PhysicalMSE regularizers stay clear of their floors.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;

struct OperatingPoint { const char* record; float tau; float sigma; };
static constexpr OperatingPoint CALIBRATED[] = {
    {"jonswap H0.27",1.28034f,0.339120f},{"jonswap H1.50",2.18073f,0.688910f},
    {"jonswap H4.00",3.59057f,1.06115f},{"jonswap H8.50",4.22184f,1.34176f},
    {"pmstokes H0.27",1.17259f,0.364507f},{"pmstokes H1.50",2.04482f,0.746216f},
    {"pmstokes H4.00",3.28991f,1.06387f},{"pmstokes H8.50",4.10027f,1.41327f},
};

static bool check(bool ok, const char* msg) {
    if (!ok) std::cerr << "FAIL: " << msg << '\n';
    return ok;
}

int main() {
    Filter f(false);
    if (f.getPseudoLaw() != PseudoAdaptationLaw::PhysicalMSE)
        return 1;

    float min_rp=1e30f, min_rv=1e30f, min_tau=0.0f, min_sigma=0.0f;
    for (const auto& op : CALIBRATED) {
        float rp=0.0f, rv=0.0f;
        f.pseudo_targets_from_law_(op.tau, op.sigma, rp, rv);
        if (rp < min_rp) { min_rp=rp; min_tau=op.tau; min_sigma=op.sigma; }
        if (rv < min_rv) min_rv=rv;
    }

    bool ok = true;
    ok &= check(min_rp > 2.0f*MIN_R_p0_std,
                "PhysicalMSE r_p floor is too close to calibrated envelope");
    ok &= check(min_rv > 2.0f*MIN_R_v0_std,
                "PhysicalMSE r_v floor is too close to calibrated envelope");

    constexpr float still_scale = 0.05f/0.27f;
    float still_rp=0.0f, still_rv=0.0f;
    f.pseudo_targets_from_law_(min_tau, min_sigma*still_scale, still_rp, still_rv);
    ok &= check(still_rp > MIN_R_p0_std,
                "PhysicalMSE r_p floor clips near-still case");
    ok &= check(still_rv > MIN_R_v0_std,
                "PhysicalMSE r_v floor clips near-still case");

    if (!ok) return 1;
    std::cout << "OU-II PhysicalMSE floor checks passed\n";
    return 0;
}
