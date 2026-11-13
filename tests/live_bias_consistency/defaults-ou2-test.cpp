#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#include "DefaultSamplingChecks.h"
int main(){default_sampling_check<ocean_imu::ou2::SeaStateFusion_OU_II<TrackerType::KALMANF>>("OU-II");}
