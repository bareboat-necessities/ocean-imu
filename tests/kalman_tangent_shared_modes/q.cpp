#define EIGEN_NON_ARDUINO
#include "StudyAdapter.h"
namespace tangent_study {
std::unique_ptr<Filter> make_q(const Config& c) {return std::make_unique<Adapter<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>,false>>(c);}
}
