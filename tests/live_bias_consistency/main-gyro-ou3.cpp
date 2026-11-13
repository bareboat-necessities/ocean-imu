#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "MainGyroRollout.h"
int gyro_rollout::ou3(int s,double t,const std::string& p){return run<ocean_imu::ou3::SeaStateFusion_OU_III<TrackerType::KALMANF>,false>(s,t,p);}
