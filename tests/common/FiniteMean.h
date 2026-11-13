#pragma once
#include <cmath>
#include <cstddef>
#include <limits>

// Missing/rejected diagnostics do not contribute zero-valued observations.
struct FiniteMean {
    double sum = 0;
    std::size_t count = 0;
    void add(double value) {
        if (std::isfinite(value)) { sum += value; ++count; }
    }
    double mean() const {
        return count ? sum / static_cast<double>(count)
                     : std::numeric_limits<double>::quiet_NaN();
    }
};
