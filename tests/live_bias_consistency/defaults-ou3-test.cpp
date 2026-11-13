#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "DefaultSamplingChecks.h"
int main(){default_sampling_check<ocean_imu::ou3::SeaStateFusion_OU_III<TrackerType::KALMANF>>("OU-III");}
