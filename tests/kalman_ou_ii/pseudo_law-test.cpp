// Pins the deployed OU-II dual pseudo-measurement PhysicalMSE law only.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;

static bool rel_near(float a, float b, float rel) {
    return std::fabs(a - b) <= rel * std::max(1e-12f, std::max(std::fabs(a), std::fabs(b)));
}
static int fail(const char* msg) { std::cerr << "FAIL: " << msg << "\n"; return 1; }
static void init(Filter& f) {
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
}

int main() {
    constexpr float tau = 2.1815f;
    constexpr float sigma = 0.6805f;
    constexpr float d = 1.02f;

    Filter f(false);
    init(f);
    if (f.getPseudoLaw() != PseudoAdaptationLaw::PhysicalMSE)
        return fail("OU-II default is not PhysicalMSE");

    float rp0=0, rv0=0, rp1=0, rv1=0;
    f.pseudo_targets_from_law_(tau, sigma, rp0, rv0);
    f.pseudo_targets_from_law_(tau, sigma*d, rp1, rv1);
    const float ap = std::log(rp1/rp0)/std::log(d);
    const float av = std::log(rv1/rv0)/std::log(d);
    if (!rel_near(ap, 4.0f/5.0f, 2e-3f) || !rel_near(av, 4.0f/5.0f, 2e-3f))
        return fail("PhysicalMSE amplitude exponent is not 4/5");

    f.pseudo_targets_from_law_(tau*d, sigma, rp1, rv1);
    const float tp = std::log(rp1/rp0)/std::log(d);
    const float tv = std::log(rv1/rv0)/std::log(d);
    if (!rel_near(tp, 19.0f/10.0f, 2e-3f) || !rel_near(tv, 9.0f/10.0f, 2e-3f))
        return fail("PhysicalMSE cadence-normalized tau exponents moved");

    for (float t : {0.3f, 1.1f, tau, 4.2f, 11.0f}) {
        for (float s : {0.05f, sigma, 3.0f}) {
            float rp=0, rv=0;
            f.pseudo_targets_from_law_(t, s, rp, rv);
            if (!rel_near(rp/rv, f.getPseudoMseRatio()*t, 1e-5f))
                return fail("PhysicalMSE channel ratio is not (C_P/C_V) tau");
        }
    }

    Filter g(false);
    init(g);
    float rpa=0, rva=0, rpb=0, rvb=0;
    g.setSigmaCoeff(0.85f);
    g.pseudo_targets_from_law_(tau, sigma, rpa, rva);
    g.setSigmaCoeff(1.70f);
    g.pseudo_targets_from_law_(tau, 2.0f*sigma, rpb, rvb);
    if (!rel_near(rpa, rpb, 1e-5f) || !rel_near(rva, rvb, 1e-5f))
        return fail("PhysicalMSE changed at fixed physical acceleration RMS");

    if (!rel_near(f.getPseudoMseCoeff(), 0.1116f, 1e-5f) ||
        !rel_near(f.getPseudoMseRatio(), 0.4611f, 1e-5f))
        return fail("deployed PhysicalMSE constants moved");

    std::cout << "OU-II deployed PhysicalMSE law passed\n";
    return 0;
}
