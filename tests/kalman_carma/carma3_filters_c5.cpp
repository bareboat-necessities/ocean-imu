// Complete-INS adapters: three-axis CARMA with CARMA(5,2) horizontal axes.
#include "carma3/FullAdapters.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_c5(const std::string& n) {
    if (n == "A52") return std::make_unique<FullCARMA3<5, false, 5, false>>(n);
    if (n == "A52S") return std::make_unique<FullCARMA3<5, true, 5, false>>(n);
    if (n == "M54") return std::make_unique<FullCARMA3<5, false, 4, false>>(n);
    return nullptr;
}

}  // namespace carma3_bench
