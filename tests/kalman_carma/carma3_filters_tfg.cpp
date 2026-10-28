// Complete-INS adapters: the two-frame Lie-group CARMA filter and its ablations.
#include "carma3/FullAdaptersTFG.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_tfg(const std::string& n) {
    using M4 = typename FullTFGCARMA3<4>::M;
    using M5 = typename FullTFGCARMA3<5>::M;
    // T41    TFG-CARMA(4,1) x/y/z: full right-Jacobian reset, CARMA3-equivalent
    //        covariance writes (the filter)
    // T41R0  no reset transport
    // T41RA  attitude-only reset (MEKF-like transport)
    // T41NW  covariance writes applied in the invariant coordinates
    // T52    TFG-CARMA(5,2) x/y/z (follow-up experiment)
    if (n == "T41") return std::make_unique<FullTFGCARMA3<4>>(n, M4::ResetMode::Full, M4::CovarianceWriteMode::Physical);
    if (n == "T41R0") return std::make_unique<FullTFGCARMA3<4>>(n, M4::ResetMode::None, M4::CovarianceWriteMode::Physical);
    if (n == "T41RA")
        return std::make_unique<FullTFGCARMA3<4>>(n, M4::ResetMode::AttitudeOnly, M4::CovarianceWriteMode::Physical);
    if (n == "T41NW") return std::make_unique<FullTFGCARMA3<4>>(n, M4::ResetMode::Full, M4::CovarianceWriteMode::Native);
    if (n == "T52") return std::make_unique<FullTFGCARMA3<5>>(n, M5::ResetMode::Full, M5::CovarianceWriteMode::Physical);
    if (n == "T52R0") return std::make_unique<FullTFGCARMA3<5>>(n, M5::ResetMode::None, M5::CovarianceWriteMode::Physical);
    return nullptr;
}

}  // namespace carma3_bench
