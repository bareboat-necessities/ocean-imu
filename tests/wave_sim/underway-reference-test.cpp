#include <cassert>
#include <cmath>
#include "util/UnderwayReference.h"

int main() {
    constexpr double speed = 5.0 * 1852.0 / 3600.0;
    constexpr double h = 0.001;
    const auto origin = UnderwayReference::at(0);
    assert(origin.displacement[0] == 0 && origin.displacement[1] == 0);
    assert(std::abs(origin.velocity[0] - speed) < 1e-14);
    // Independent finite differences test the primitive and its derivatives,
    // including late-record cancellation and both steering signs.
    for (double t : {0.0, 17.3, 30.0, 75.0, 120.0, 1199.995}) {
        const auto left = UnderwayReference::at(t - h);
        const auto center = UnderwayReference::at(t);
        const auto right = UnderwayReference::at(t + h);
        double norm2 = 0;
        for (int axis = 0; axis < 3; ++axis) {
            const auto i = static_cast<std::size_t>(axis);
            assert(std::abs((right.displacement[i] - left.displacement[i]) / (2*h) - center.velocity[i]) < 1e-8);
            assert(std::abs((right.velocity[i] - left.velocity[i]) / (2*h) - center.acceleration[i]) < 1e-10);
            norm2 += center.velocity[i] * center.velocity[i];
        }
        assert(std::abs(std::sqrt(norm2) - speed) < 1e-13);
    }
    const auto cycle = UnderwayReference::at(120);
    assert(std::abs(cycle.displacement[1]) < 1e-12);
    assert(cycle.displacement[0] > 308.0 && cycle.displacement[0] < 308.7);
}
