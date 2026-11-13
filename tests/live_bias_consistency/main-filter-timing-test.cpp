#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include "StandardHeaders.h"
#include <Eigen/Dense>
#define private public
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#undef private
#include "GyroPredictionChecks.h"
struct TfgOuProbe : ocean_imu::kalman::Kalman3D_Wave_TFG<double> {
 using Base=ocean_imu::kalman::Kalman3D_Wave_TFG<double>;
 TfgOuProbe(const V& a,const V& g,const V& m):Base(g.x()){set_Racc_std(a);set_Rmag_std(m);}
 void initialize_from_attitude(const Q& q,double,double){state().R=q.toRotationMatrix();}
 Q quaternion_boat()const{return quaternion();}
 void set_gyro_noise_density_rad_sqrt_s(const V& g){Base::set_gyro_noise_density_rad_sqrt_s(g.x());}
};
int main(){
 tests<ocean_imu::kalman::Kalman3D_Wave_OU_II<double>>("OU-II");
 tests<ocean_imu::kalman::Kalman3D_Wave_OU_III<double>>("OU-III");
 tests<TfgOuProbe>("TFG-OU");
 std::cout<<"main-filter timing tests PASS\n";
}
