#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#include "DefaultSamplingChecks.h"
int main(){default_sampling_check<ocean_imu::tfg::SeaStateFusionFilter_TFG<>>("TFG");}
