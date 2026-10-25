// Complete-INS adapters: shipping OU-III, the PR #704 hybrids, and the
// three-axis class with the OU-III horizontal chain (equivalence control).
#include "carma3/FullAdapters.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_ou(const std::string& n) {
    if (n == "OU3") return std::make_unique<FullOU3>();
    if (n == "Z41") return std::make_unique<FullHybrid<4, false>>(n);
    if (n == "Z52") return std::make_unique<FullHybrid<5, false>>(n);
    if (n == "Q41") return std::make_unique<FullCARMA3<0, false, 4, false>>(n);
    if (n == "Q52") return std::make_unique<FullCARMA3<0, false, 5, false>>(n);
    return nullptr;
}

}  // namespace carma3_bench
