#define EIGEN_NON_ARDUINO
#include "Adapter.h"
namespace followup {std::unique_ptr<Experiment> make_tfg(const Options& o,const V& b) {
return std::make_unique<Adapter<ocean_imu::carma_shared_modes::KalmanTfgCarmaSharedModes<float>,true>>(o,b);}}
