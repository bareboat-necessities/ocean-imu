#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#include "MainGyroRollout.h"
int gyro_rollout::ou2(int s,double t,const std::string& p){return run<ocean_imu::ou2::SeaStateFusion_OU_II<TrackerType::KALMANF>,false>(s,t,p);}
