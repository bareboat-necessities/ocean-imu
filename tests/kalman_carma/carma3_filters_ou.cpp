// Complete-INS adapter: the shipping OU-III baseline.
#include "carma3/FullAdapters.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_ou(const std::string& n) {
    if (n == "OU3") return std::make_unique<FullOU3>();
    return nullptr;
}

}  // namespace carma3_bench
