#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "imu_calibrate/FrameAlignment.h"
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <random>
using namespace atoms3r_ical;
struct Kv {
  std::map<std::string,std::vector<uint8_t>> records;
  size_t getBytesLength(const char* k){return records[k].size();}
  size_t getBytes(const char* k,void* p,size_t n){const auto& v=records[k];if(n!=v.size())return 0;memcpy(p,v.data(),n);return n;}
  size_t putBytes(const char* k,const void* p,size_t n){const auto* b=static_cast<const uint8_t*>(p);records[k]={b,b+n};return n;}
  bool remove(const char* k){records.erase(k);return true;}
};
using Store=ImuCalStoreT<Kv>;
static int errors=0;static void check(bool ok,const char* s){if(!ok){++errors;std::fprintf(stderr,"FAIL %s\n",s);}}
static ImuCalBlobV6 record() {
  ImuCalBlobV6 b;b.accel_ok=1;b.accel_S[0]=1.02f;b.accel_S[4]=.99f;b.accel_S[8]=1.01f;
  b.accel_b0[0]=.1f;b.accel_T_lo=b.accel_T_hi=b.accel_T0;
  b.accel_coeff_crc=accelCoeffCrc(b);
  imu_cal::GyroCalibration<float> g;g.ok=g.biasT.ok=g.matrix_ok=true;
  g.biasT.T_lo=g.biasT.T_hi=g.biasT.T0;
  g.S<<1.03f,.02f,0,0,.98f,.01f,-.01f,0,1.02f;
  g.rotation_condition=3;g.rotation_train_segments=8;g.rotation_validation_segments=4;
  g.biasT.b0<<.01f,-.02f,.005f;fillGyroFromFit(b,g);
  b.mag_ok=1;b.mag_A[0]=1.1f;b.mag_A[4]=.9f;b.mag_A[8]=1;b.mag_field_uT=45;b.mag_b[1]=2;
  b.mag_alignment_ok=b.housing_alignment_ok=b.vessel_alignment_ok=1;
  b.alignment_condition=3;b.alignment_train_poses=9;b.alignment_validation_poses=4;
  b.sensor_id_lo=12;b.sensor_id_hi=34;b.imu_type=5;
  return b;
}
int main(){
  ImuCalBlobV6 b=record();
  const Matrix3f rim=Eigen::AngleAxisf(.04f,Vector3f::UnitY()).toRotationMatrix();
  const Matrix3f rhi=Eigen::AngleAxisf(.2f,Vector3f::UnitX()).toRotationMatrix();
  const Matrix3f rvh=Eigen::AngleAxisf(-.3f,Vector3f::UnitZ()).toRotationMatrix();
  mat_to_rowmajor9_(rim,b.R_IM);mat_to_rowmajor9_(rhi,b.R_HI);mat_to_rowmajor9_(rvh,b.R_VH);
  Store store;ImuCalBlobV6 rb;check(store.saveVerified(b,rb),"new schema round-trip");
  RuntimeCals runtime;runtime.rebuildFromBlob(rb);check(!calibrationUnusableReason(true,runtime),"full calibration usable");
  const Vector3f x(1,2,-9);
  check(runtime.applyAccel(x,25).isApprox(rvh*rhi*runtime.acc.apply(x,25),1e-6f),"accel intrinsic then housing then vessel");
  check(runtime.applyGyro(x,25).isApprox(rvh*rhi*runtime.gyr.apply(x,25),1e-6f),"full gyro matrix applied once");
  check(runtime.applyMag(x).isApprox(rvh*rhi*rim*runtime.mag.apply(x),1e-6f),"mag alignment applied once in correct order");
  b=rb;b.gyro_train_segments=7;b.gyro_coeff_crc=gyroCoeffCrc(b);
  check(!validateBlob(Store::sealed_(b)),"too few short-plan training rotations rejected");
  b=rb;b.alignment_train_poses=8;
  check(!validateBlob(Store::sealed_(b)),"too few short-plan alignment poses rejected");
  b=rb;b.gyro_validation_segments=3;b.gyro_coeff_crc=gyroCoeffCrc(b);
  check(!validateBlob(Store::sealed_(b)),"short plan still requires four independent validation rotations");
  b=rb;b.version=5;b.crc=computeBlobCrc(b);
  check(!validateBlob(b) && !store.save(b),"unsupported header cannot be sealed as new");
  Store old;std::vector<uint8_t> v5(372,0);old.kv.putBytes("blob_m5v5",v5.data(),v5.size());
  check(!old.load(b),"old key never used");
  old.kv.putBytes(old.kKeyV6,&rb,sizeof(rb)-4);check(!old.load(b),"truncated new record rejected");
  b=rb;b.R_IM[0]+=.01f;check(!validateBlob(b),"rotation integrity protected");
  runtime.rebuildFromBlob(b);check(!runtime.record_valid && !runtime.acc.ok,"invalid record clears earlier runtime");
  b=rb;b.R_HI[0]=-b.R_HI[0];check(!validateBlob(Store::sealed_(b)),"reflection cannot become mounting rotation");
  b=rb;b.gyro_S[0]=2;b.gyro_coeff_crc=gyroCoeffCrc(b);
  check(!validateBlob(Store::sealed_(b)),"implausible gyro scale rejected");
  Vector3f src[3]={Vector3f::UnitX(),Vector3f::UnitY(),Vector3f::UnitZ()},dst[3];
  for(int i=0;i<3;++i)dst[i]=rvh*rhi*src[i];
  Matrix3f solved;
  check(imu_cal::solveFrameAlignment(src,dst,3,solved) && solved.isApprox(rvh*rhi,1e-5f),"physical fixture alignment");
  src[1]=src[2]=src[0];check(!imu_cal::solveFrameAlignment(src,dst,3,solved),"gravity alone cannot determine mounting yaw");
  // Actual wizard pose count and train/validation split, not hundreds of poses.
  atoms3r_magdiag::Pose train[9],validation[4];std::mt19937 rng(19);std::normal_distribution<float> n;
  for(int i=0;i<13;++i) {
    Eigen::Quaternionf q(n(rng),n(rng),n(rng),n(rng));q.normalize();
    atoms3r_magdiag::Pose pose;pose.a=q.conjugate()*Vector3f(0,0,-9.80256f);
    pose.m_raw=rim.transpose()*(q.conjugate()*Vector3f(20,0,43));
    if(i<9)train[i]=pose;else validation[i-9]=pose;
  }
  atoms3r_magdiag::MagModel model;
  const auto alignment=atoms3r_magdiag::solveAlignment(train,9,model);float rms=NAN;
  check(alignment.ok && alignment.R.isApprox(rim,1e-4f) &&
    atoms3r_magdiag::validateAlignment(alignment,validation,4,model,rms),"wizard-sized independent magnetic alignment");
  for(auto& p:train)p=train[0];
  check(!atoms3r_magdiag::solveAlignment(train,9,model).ok,"ridge cannot qualify degenerate alignment");
  std::printf("calibration frames: blob_bytes=%zu independent_alignment_rms=%.8g\n",sizeof(rb),double(rms));
  return errors?1:0;
}
