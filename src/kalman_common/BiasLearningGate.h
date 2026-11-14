#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Experimental finite-memory excitation veto. This CONDITIONAL information
// model assumes a constant world attitude error, known field and zero-mean
// wave acceleration in the window. It is NOT the full marine observability
// Gramian. It can remove bias gain but never contributes a measurement or
// adds information to the covariance. The caller must use full Joseph with
// the modified gain, including all cross-covariances.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#else
#include <ArduinoEigenDense.h>
#endif
#include "lie/SO3Jacobians.h"
#include <algorithm>
#include <cmath>

namespace ocean_imu::kalman::filter_detail {
template<typename T> struct BiasLearningConfig {
    bool enabled=false;
    T horizon_s=T(30), commit_s=T(.1), correlation_s=T(1);
    T tilt_scale_rad=T(.05), bias_scale_mps2=T(.05);
    T force_std_mps2=T(.15), field_std=T(8);
    T information_threshold=T(1);
};
template<typename T> class BiasLearningGate {
public:
    using V=Eigen::Matrix<T,3,1>;
    using M=Eigen::Matrix<T,3,3>;
    using I6=Eigen::Matrix<T,6,6>;
    void configure(const BiasLearningConfig<T>& c) { cfg_=c; reset(); }
    void reset() { info_.setZero(); scale_.setIdentity(); eigen_.setZero(); clock_=T(0); age_=T(0); }
    bool enabled() const { return cfg_.enabled; }
    const M& scale() const { return scale_; }
    const V& information_eigenvalues() const { return eigen_; }
    T age() const { return age_; }
    // Called once per accepted IMU sample. Information is accumulated at a
    // bounded rate with an explicit effective correlation time; high-rate
    // copies are not counted as fresh stationary poses. Motion power only
    // inflates the conditional force residual model.
    void observe(const M& R,const V& field,const V& acceleration,T dt,T gravity) {
        if(!cfg_.enabled || !(dt>T(0)) || !std::isfinite(dt) || !R.allFinite() || !field.allFinite() || !acceleration.allFinite()) return;
        clock_+=dt;
        if(clock_<std::max(T(.01),cfg_.commit_s)) return;
        const T h=clock_; clock_=T(0); age_+=h;
        const T fade=std::exp(-h/std::max(T(.1),cfg_.horizon_s));
        const T weight=(T(1)-fade)*cfg_.horizon_s/std::max(T(.01),cfg_.correlation_s);
        Eigen::Matrix<T,3,6> Ha,Hm;
        const T variance=std::max(T(1e-6),cfg_.force_std_mps2*cfg_.force_std_mps2+acceleration.squaredNorm()/T(3));
        Ha.template leftCols<3>()=lie::skew<T>(V(0,0,gravity))*cfg_.tilt_scale_rad/std::sqrt(variance);
        Ha.template rightCols<3>()=-R*cfg_.bias_scale_mps2/std::sqrt(variance);
        Hm.setZero();
        Hm.template leftCols<3>()=-lie::skew<T>(field)*cfg_.tilt_scale_rad/std::max(T(1e-3),cfg_.field_std);
        info_=fade*info_+weight*(Ha.transpose()*Ha+Hm.transpose()*Hm);
        const M A=info_.template topLeftCorner<3,3>();
        const M B=info_.template topRightCorner<3,3>();
        Eigen::SelfAdjointEigenSolver<M> eigA(A);
        if(eigA.info()!=Eigen::Success) { scale_.setZero(); return; }
        V ai=eigA.eigenvalues();
        const T tol=std::max(T(1e-10),ai.maxCoeff()*T(1e-6));
        for(int i=0;i<3;++i) ai(i)=ai(i)>tol?T(1)/ai(i):T(0);
        const M Ai=eigA.eigenvectors()*ai.asDiagonal()*eigA.eigenvectors().transpose();
        M C=info_.template bottomRightCorner<3,3>()-B.transpose()*Ai*B;
        C=(C+C.transpose()).eval()*T(.5);
        Eigen::SelfAdjointEigenSolver<M> es(C);
        if(es.info()!=Eigen::Success) { scale_.setZero(); return; }
        eigen_=es.eigenvalues().cwiseMax(T(0));
        const V gain=(eigen_.array()/(eigen_.array()+std::max(T(1e-6),cfg_.information_threshold))).matrix();
        scale_=es.eigenvectors()*gain.asDiagonal()*es.eigenvectors().transpose();
    }
private:
    BiasLearningConfig<T> cfg_{};
    I6 info_=I6::Zero();
    M scale_=M::Identity();
    V eigen_=V::Zero();
    T clock_=T(0),age_=T(0);
};
}
