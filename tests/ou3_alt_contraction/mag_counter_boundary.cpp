// Host boundary correspondence, not a reset-to-INT_MAX reachability replay.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

int main(int argc, char**) {
    static_assert(std::numeric_limits<int>::max() == 2147483647);
    Fusion f;
    f.begin(Fusion::Config{});
    const Eigen::Vector3f field(15.0f, 0.0f, 20.0f);
    int samples = 0;
    for (; samples < 40000 && !(f.isLive() && f.hasMagNorthLock()); ++samples) {
        f.update(.005f, Eigen::Vector3f::Zero(), Eigen::Vector3f(0, 0, -g_std));
        if (samples % 8 == 0) f.updateMag(field);
    }
    if (!(f.isLive() && f.hasMagNorthLock())) {
        std::cerr << "quiet public-API prefix did not reach gauged Live\n";
        return 2;
    }
    // ONLY the count is installed for a two-edge boundary check. The universal
    // counter-projection induction, not this installation, proves the burst.
    auto& inner = f.raw();
    inner.mag_updates_applied_ = std::numeric_limits<int>::max() - 1;
    f.updateMag(field);
    if (inner.mag_updates_applied_ != std::numeric_limits<int>::max()) return 3;
    std::cout << "gauged_Live_samples=" << samples
              << " last_defined_count=" << inner.mag_updates_applied_ << std::endl;
    if (argc > 1) {
        f.updateMag(field); // UBSan must stop at the unchanged shipping ++.
        return 4;
    }
    return 0;
}
