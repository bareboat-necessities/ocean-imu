#define EIGEN_NON_ARDUINO
#include "StudyAdapter.h"
namespace tangent_study {
std::unique_ptr<Filter> make_independent(const Config& c) {return std::make_unique<Adapter<ocean_imu::kalman::Kalman3D_Wave_CARMA3<float>,false,false>>(c);}
}
