// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#define EIGEN_NON_ARDUINO
#include "CsmBenchAdapter.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_TfgCarmaSharedModes.h"
namespace carma3_bench {
std::unique_ptr<IEstimator3> make_tfg_csm(const std::string& n) {
    if (n.rfind("TCSM", 0) != 0) return nullptr;
    CsmSpec sp;
    if (!parse_csm_spec("CSM" + n.substr(4), sp)) return nullptr;
    return std::make_unique<FullSharedModes<csm::SeaStateFusion_TfgCarmaSharedModes>>(n, sp);
}
}
