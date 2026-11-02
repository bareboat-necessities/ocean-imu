#pragma once
// Copyright 2026, Mikhail Grushinskiy
// Extension of the existing GyroCalibration, not an alternative correction path.
// Tedaldi/Pretto/Menegatti: integrate measured rate between calibrated static
// gravity endpoints. Nine matrix parameters; bias comes from GyroCapture.
#include "imu_calibrate/CalibrateIMU.h"
#include <algorithm>
#include <limits>

namespace imu_cal {
enum class GyroRotationFail : uint8_t {
  OK, BAD_SAMPLE, GAP, CAPACITY, ENDPOINT, TOO_FEW, GEOMETRY, SOLVE, PHYSICAL, VALIDATION
};
inline const char* gyroRotationReason(GyroRotationFail r) {
  switch(r) {
    case GyroRotationFail::OK:return "OK";
    case GyroRotationFail::GAP:return "Rotation sample gap";
    case GyroRotationFail::CAPACITY:return "Rotation buffer full";
    case GyroRotationFail::ENDPOINT:return "Rotate to another pose";
    case GyroRotationFail::TOO_FEW:return "More rotation segments";
    case GyroRotationFail::GEOMETRY:return "Need varied tilt axes";
    case GyroRotationFail::SOLVE:return "Gyro solve failed";
    case GyroRotationFail::PHYSICAL:return "Gyro matrix implausible";
    case GyroRotationFail::VALIDATION:return "Gyro validation failed";
    default:return "Invalid rotation sample";
  }
}

struct GyroRotationCfg {
  double max_gap_s=.020, max_rate=8.0, max_temperature_span=2.0;
  double min_tilt_rad=.15, max_condition=1e5, min_information=.002;
  double huber_rad=.025, max_train_rms=.025, max_validation_rms=.025, max_validation_error=.06;
  int min_train=12, min_validation=4, max_iterations=25;
};

// Bound storage. Rejected/aborted segments release their samples. Capture raw
// rates and endpoint specific forces so the FINAL accel/bias fit can be used.
template<int MAX_SAMPLES=8192,int MAX_SEGMENTS=24>
class GyroRotationFit {
public:
  using V=Eigen::Vector3d; using M=Eigen::Matrix3d;
  using J=Eigen::Matrix<double,3,9>; using H9=Eigen::Matrix<double,9,9>;
  struct Sample {float w[3],temp,dt;};
  struct Segment {int first=0,count=0;Eigen::Vector3f a0=Eigen::Vector3f::Zero(),a1=a0;
    float t0=25,t1=25;bool validation=false;};
  GyroRotationCfg cfg{};
  GyroRotationFail reason=GyroRotationFail::OK;
  void clear() {n_=ns_=0;active_=false;reason=GyroRotationFail::OK;}
  int segments() const {return ns_;}
  int samples() const {return n_;}
  bool beginSegment(uint32_t us,const Eigen::Vector3f& a,float temp,bool validation) {
    if(active_ || ns_>=MAX_SEGMENTS || !a.allFinite() || !std::isfinite(temp) || a.norm()<1)
      return fail_(GyroRotationFail::BAD_SAMPLE);
    pending_={};pending_.first=n_;pending_.a0=a;pending_.t0=temp;pending_.validation=validation;
    last_us_=us;have_rate_=false;active_=true;reason=GyroRotationFail::OK;return true;
  }
  bool add(uint32_t us,const Eigen::Vector3f& w,float temp) {
    if(!active_)return false;
    if(!w.allFinite() || !std::isfinite(temp) || w.norm()>cfg.max_rate ||
       std::fabs(temp-pending_.t0)>cfg.max_temperature_span) return abort_(GyroRotationFail::BAD_SAMPLE);
    const uint32_t du=us-last_us_;
    if(!du || du>uint32_t(cfg.max_gap_s*1e6))return abort_(GyroRotationFail::GAP);
    if(!have_rate_) {last_rate_=w;last_temp_=temp;last_us_=us;have_rate_=true;return true;}
    if(n_>=MAX_SAMPLES)return abort_(GyroRotationFail::CAPACITY);
    Sample& s=data_[n_++];
    const Eigen::Vector3f midpoint=.5f*(last_rate_+w);
    for(int k=0;k<3;++k)s.w[k]=midpoint[k];
    s.temp=.5f*(last_temp_+temp);s.dt=du*1e-6f;
    last_rate_=w;last_temp_=temp;last_us_=us;return true;
  }
  // GyroCapture READY proves the final useful_seconds were quiet. Removing
  // just that tail saves RAM/integration work without changing the endpoint.
  bool finishSegment(const Eigen::Vector3f& a,float temp,double quiet_tail_s=0) {
    if(!active_ || !a.allFinite() || !std::isfinite(temp))return abort_(GyroRotationFail::BAD_SAMPLE);
    double removed=0;
    while(n_>pending_.first && removed+data_[n_-1].dt<=quiet_tail_s)removed+=data_[--n_].dt;
    pending_.count=n_-pending_.first;pending_.a1=a;pending_.t1=temp;
    if(pending_.count<10)return abort_(GyroRotationFail::ENDPOINT);
    seg_[ns_++]=pending_;active_=false;return true;
  }
  void cancelSegment() {if(active_) {n_=pending_.first;active_=false;}}

  // Public for finite-difference regression of the analytic SO(3) sensitivity.
  static M skew(const V& v) {M x;x<<0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0;return x;}
  static void step(V& u,J& D,const M& S,const V& rate_dt) {
    const V theta=-S*rate_dt;const double t2=theta.squaredNorm();
    const M X=skew(theta),X2=X*X;
    double a,b,c;
    if(t2<1e-8) {a=1-t2/6;b=.5-t2/24;c=1.0/6-t2/120;}
    else {const double t=std::sqrt(t2);a=std::sin(t)/t;b=(1-std::cos(t))/t2;c=(1-a)/t2;}
    const M R=M::Identity()+a*X+b*X2,Jl=M::Identity()+b*X+c*X2;
    u=(R*u).eval();D=(R*D).eval();
    const M B=skew(u)*Jl;
    for(int row=0;row<3;++row)for(int col=0;col<3;++col)
      D.col(3*row+col).noalias()+=B.col(row)*rate_dt[col];
  }
  template<typename T>
  bool fit(const AccelCalibration<T>& acc,GyroCalibration<T>& gyro,
           void (*idle)(void*)=nullptr,void* context=nullptr) {
    idle_=idle;idle_context_=context;
    gyro.matrix_ok=false;reason=GyroRotationFail::OK;
    if(active_ || !acc.ok || !gyro.ok)return fail_(GyroRotationFail::BAD_SAMPLE);
    int nt=0,nv=0;for(int i=0;i<ns_;++i)seg_[i].validation?++nv:++nt;
    if(nt<cfg.min_train || nv<cfg.min_validation)return fail_(GyroRotationFail::TOO_FEW);
    M S=M::Identity();H9 H;Eigen::Matrix<double,9,1> g;
    double cost=0;
    if(!normal_(S,acc,gyro,H,g,cost))return false;
    double lambda=1e-5;
    for(int it=0;it<cfg.max_iterations;++it) {
      H9 damped=H;damped.diagonal().array()+=lambda;
      const auto ldlt=damped.ldlt();
      if(ldlt.info()!=Eigen::Success || !ldlt.isPositive())return fail_(GyroRotationFail::SOLVE);
      const Eigen::Matrix<double,9,1> dx=ldlt.solve(-g);
      if(!dx.allFinite())return fail_(GyroRotationFail::SOLVE);
      M candidate=S;for(int r=0;r<3;++r)for(int c=0;c<3;++c)candidate(r,c)+=dx[3*r+c];
      H9 Hn;Eigen::Matrix<double,9,1> gn;double cn=0;
      if(!normal_(candidate,acc,gyro,Hn,gn,cn))return false;
      if(cn<=cost) {S=candidate;H=Hn;g=gn;cost=cn;lambda=std::max(1e-10,lambda*.2);
        if(dx.norm()<1e-8)break;}
      else lambda*=10;
    }
    Eigen::SelfAdjointEigenSolver<H9> es(H);
    if(es.info()!=Eigen::Success || es.eigenvalues()[0]<cfg.min_information ||
       es.eigenvalues()[8]/es.eigenvalues()[0]>cfg.max_condition)return fail_(GyroRotationFail::GEOMETRY);
    Eigen::JacobiSVD<M> svd(S);
    if(!S.allFinite() || S.determinant()<=0 || svd.singularValues().minCoeff()<.8 ||
       svd.singularValues().maxCoeff()>1.25 || (S-M::Identity()).cwiseAbs().maxCoeff()>.20)
      return fail_(GyroRotationFail::PHYSICAL);
    // Validate the exact stored float/T matrix, including independent segments.
    const M stored=S.template cast<T>().template cast<double>();
    double train=0,validation=0,vmax=0;
    for(int i=0;i<ns_;++i) {
      V e;J D;if(!residual_(seg_[i],stored,acc,gyro,e,D))return false;
      const double angle=2*std::asin(std::min(1.0,e.norm()*.5));
      if(seg_[i].validation){validation+=angle*angle;vmax=std::max(vmax,angle);}else train+=angle*angle;
    }
    train=std::sqrt(train/nt);validation=std::sqrt(validation/nv);
    if(train>cfg.max_train_rms || validation>cfg.max_validation_rms || vmax>cfg.max_validation_error)
      return fail_(GyroRotationFail::VALIDATION);
    gyro.S=stored.template cast<T>();gyro.matrix_ok=true;
    gyro.rotation_train_rms=T(train);gyro.rotation_validation_rms=T(validation);
    gyro.rotation_validation_max=T(vmax);gyro.rotation_condition=T(es.eigenvalues()[8]/es.eigenvalues()[0]);
    gyro.rotation_train_segments=nt;gyro.rotation_validation_segments=nv;return true;
  }
private:
  bool fail_(GyroRotationFail r) {reason=r;return false;}
  bool abort_(GyroRotationFail r) {cancelSegment();return fail_(r);}
  template<typename T>
  bool residual_(const Segment& s,const M& S,const AccelCalibration<T>& acc,
                 const GyroCalibration<T>& gyro,V& e,J& D) {
    V u=acc.apply(s.a0.template cast<T>(),T(s.t0)).template cast<double>();
    V target=acc.apply(s.a1.template cast<T>(),T(s.t1)).template cast<double>();
    if(!u.allFinite() || !target.allFinite() || std::fabs(u.norm()-acc.g)>.25 ||
       std::fabs(target.norm()-acc.g)>.25)return fail_(GyroRotationFail::ENDPOINT);
    u.normalize();target.normalize();
    if(std::acos(clamp(u.dot(target),-1.0,1.0))<cfg.min_tilt_rad)
      return fail_(GyroRotationFail::ENDPOINT);
    D.setZero();
    for(int k=s.first;k<s.first+s.count;++k) {
      // Embedded workers must periodically let the idle task run. The hook
      // has no mathematical effect and introduces no RTOS dependency here.
      if(idle_ && ((k-s.first)&127)==0)idle_(idle_context_);
      const Sample& x=data_[k];
      if(gyro.thermal==GyroThermal::UNLEARNED &&
         std::fabs(double(x.temp)-double(gyro.biasT.T0))>cfg.max_temperature_span)
        return fail_(GyroRotationFail::BAD_SAMPLE);
      const V raw(x.w[0],x.w[1],x.w[2]);
      const V rate=raw-gyro.biasT.bias(T(x.temp)).template cast<double>();
      step(u,D,S,rate*double(x.dt));
    }
    e=u-target;return e.allFinite() && D.allFinite();
  }
  template<typename T>
  bool normal_(const M& S,const AccelCalibration<T>& acc,const GyroCalibration<T>& gyro,
               H9& H,Eigen::Matrix<double,9,1>& g,double& cost) {
    H.setZero();g.setZero();cost=0;
    for(int i=0;i<ns_;++i)if(!seg_[i].validation) {
      V e;J D;if(!residual_(seg_[i],S,acc,gyro,e,D))return false;
      const double n=e.norm(),w=n<=cfg.huber_rad?1:cfg.huber_rad/n;
      H.noalias()+=w*D.transpose()*D;g.noalias()+=w*D.transpose()*e;
      cost+=n<=cfg.huber_rad?.5*n*n:cfg.huber_rad*(n-.5*cfg.huber_rad);
    }
    return true;
  }
  Sample data_[MAX_SAMPLES];Segment seg_[MAX_SEGMENTS],pending_{};
  int n_=0,ns_=0;bool active_=false,have_rate_=false;uint32_t last_us_=0;
  Eigen::Vector3f last_rate_=Eigen::Vector3f::Zero();float last_temp_=25;
  void (*idle_)(void*)=nullptr;void* idle_context_=nullptr;
};
} // namespace imu_cal
