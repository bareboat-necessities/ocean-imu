#define EIGEN_NON_ARDUINO
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"
#include "StudyAdapter.h"
namespace tangent_study {
std::unique_ptr<Filter> make_tfg(const Config& c) {return std::make_unique<Adapter<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>,true>>(c);}
#ifndef TANGENT_BASELINE
std::unique_ptr<Filter> make_tangent(const Config& c) {
    Config d=c;d.tangent_geometry=true;return make_tfg(d);
}
#endif
}
