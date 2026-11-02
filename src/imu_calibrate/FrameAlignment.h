#pragma once
// Copyright 2026, Mikhail Grushinskiy
#include "imu_calibrate/CalibrateIMU.h"
namespace imu_cal {
// Wahba alignment of measured, known physical directions. Gravity alone does
// not determine yaw: two non-collinear directions are mandatory. This utility
// is shared by fixture/installation tooling, separate from intrinsic fitting.
template<typename T>
bool solveFrameAlignment(const Eigen::Matrix<T,3,1>* source,const Eigen::Matrix<T,3,1>* target,
                         int n,Eigen::Matrix<T,3,3>& rotation,T max_rms_rad=T(.02)) {
  using M=Eigen::Matrix<T,3,3>;using V=Eigen::Matrix<T,3,1>;
  if(!source || !target || n<2)return false;
  M B=M::Zero();
  for(int i=0;i<n;++i) {
    if(!source[i].allFinite() || !target[i].allFinite() || source[i].norm()<T(1e-6) ||
       target[i].norm()<T(1e-6))return false;
    B.noalias()+=target[i].normalized()*source[i].normalized().transpose();
  }
  Eigen::JacobiSVD<M> svd(B,Eigen::ComputeFullU|Eigen::ComputeFullV);
  if(svd.singularValues()[1]<T(.01)*svd.singularValues()[0])return false;
  M sign=M::Identity();sign(2,2)=(svd.matrixU()*svd.matrixV().transpose()).determinant();
  const M R=svd.matrixU()*sign*svd.matrixV().transpose();
  T s=0;
  for(int i=0;i<n;++i){const V p=R*source[i].normalized();
    const T a=std::acos(clamp(p.dot(target[i].normalized()),T(-1),T(1)));s+=a*a;}
  if(!R.allFinite() || std::sqrt(s/n)>max_rms_rad)return false;
  rotation=R;return true;
}
} // namespace imu_cal
