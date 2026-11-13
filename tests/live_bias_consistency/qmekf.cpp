#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "Adapter.h"
namespace followup {std::unique_ptr<Experiment> make_qmekf(const Options& o,const V& b) {
return std::make_unique<Adapter<ocean_imu::carma_shared_modes::KalmanCarmaSharedModes<float>,false>>(o,b);}}
