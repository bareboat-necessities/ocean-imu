// Complete-INS adapters: coupled three-axis CARMA priors (Experiment A and the modal form).
#include "carma3/FullAdaptersCoupled.h"

namespace carma3_bench {

namespace {
namespace kc = ocean_imu::kalman::carma;
template<class M> using CoupledMekf = ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false, true, true, M>;
}  // namespace

std::unique_ptr<IEstimator3> make_full_coupled(const std::string& n) {
    const size_t c = n.find(':');
    if (c == std::string::npos) return nullptr;
    const std::string model = n.substr(0, c);
    CoupledSpec sp;
    if (model != "COR" && model != "MOD") return nullptr;
    if (!parse_coupled_spec(n.substr(c + 1), sp)) return nullptr;
    if (model == "MOD") return std::make_unique<FullCoupled<CoupledMekf<kc::CarmaModalMotion<float>>>>(n, sp);
    return std::make_unique<FullCoupled<CoupledMekf<kc::CarmaCorrelatedMotion<float>>>>(n, sp);
}

}  // namespace carma3_bench
