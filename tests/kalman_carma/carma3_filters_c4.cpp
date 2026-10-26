// Complete-INS adapters: three-axis CARMA with CARMA(3,0)/(4,1) horizontal axes.
#include <cstdlib>

#include "carma3/FullAdapters.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_c4(const std::string& n) {
    if (n == "A41") return std::make_unique<FullCARMA3<4, false, 4, false>>(n);
    if (n == "A41S") return std::make_unique<FullCARMA3<4, true, 4, false>>(n);
    if (n == "M45") return std::make_unique<FullCARMA3<4, false, 5, false>>(n);
    if (n == "H41") return std::make_unique<FullCARMA3<4, false, 4, false>>(n, &CARMA41_SWAY_SHAPE);
    if (n == "A30") return std::make_unique<FullCARMA3<3, false, 3, false>>(n);
    // A41C:<r> -- horizontal driving-noise correlation r (coupled x-y prior).
    if (n.rfind("A41C:", 0) == 0) {
        const float r = static_cast<float>(std::strtod(n.c_str() + 5, nullptr));
        return std::make_unique<FullCARMA3<4, false, 4, false>>(n, nullptr, r);
    }
    return nullptr;
}

}  // namespace carma3_bench
