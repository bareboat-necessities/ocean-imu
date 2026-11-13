#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#include "MainGyroRollout.h"
int gyro_rollout::tfg(int s,double t,const std::string& p){return run<ocean_imu::tfg::SeaStateFusionFilter_TFG<>,true>(s,t,p);}
