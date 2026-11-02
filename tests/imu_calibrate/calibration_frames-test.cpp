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
  // Explicit basic choices are usable after verified save/reload. A missing
  // full fit without operator consent must still stop the startup gate.
  for(int profile:{1,2,3}) {
    ImuCalBlobV6 basic=rb;
    if(profile&1) {
      imu_cal::GyroCalibration<float> bias;bias.ok=bias.biasT.ok=true;
      bias.biasT.T_lo=bias.biasT.T_hi=bias.biasT.T0;bias.biasT.b0=Vector3f(.012f,-.008f,.003f);
      fillGyroFromFit(basic,bias,true);
    }
    if(profile&2) {
      basic.mag_alignment_ok=0;basic.mag_alignment_skipped=1;
      mat_to_rowmajor9_(Matrix3f::Identity(),basic.R_IM);
      basic.alignment_train_sd=basic.alignment_validation_sd=basic.alignment_condition=0;
      basic.alignment_train_poses=basic.alignment_validation_poses=0;
    }
    ImuCalBlobV6 readback;
    check(store.saveVerified(basic,readback),"basic profile must reach verified persistence");
    RuntimeCals loaded;loaded.rebuildFromBlob(readback);
    check(!calibrationUnusableReason(true,loaded) && loaded.gyro_matrix_skipped==bool(profile&1) &&
      loaded.mag_alignment_skipped==bool(profile&2),"basic profile must survive reload without reentering the wizard");
    if(profile&1) {
      const Vector3f input(.3f,.2f,-.4f);
      check(!loaded.gyr.matrix_ok && loaded.applyGyro(input,25).isApprox(rvh*rhi*(input-loaded.gyr.biasT.b0),1e-6f),
        "basic mode applies measured bias and nominal scale/axes exactly once");
      auto unchosen=readback;unchosen.gyro_matrix_skipped=0;unchosen.gyro_coeff_crc=gyroCoeffCrc(unchosen);
      loaded.rebuildFromBlob(Store::sealed_(unchosen));
      check(calibrationUnusableReason(true,loaded)!=nullptr,"missing gyro fit cannot silently opt into basic mode");
      auto failed=readback;failed.gyro_S[0]=1.1f;failed.gyro_coeff_crc=gyroCoeffCrc(failed);
      check(!validateBlob(Store::sealed_(failed)),"explicit skip cannot apply an unqualified nonidentity gyro matrix");
      failed=readback;failed.gyro_train_segments=12;failed.gyro_coeff_crc=gyroCoeffCrc(failed);
      check(!validateBlob(Store::sealed_(failed)),"skipped rotation fit cannot claim training qualification");
      failed=readback;failed.gyro_ok=0;failed.gyro_coeff_crc=gyroCoeffCrc(failed);
      loaded.rebuildFromBlob(Store::sealed_(failed));
      check(calibrationUnusableReason(true,loaded)!=nullptr,"basic mode cannot bypass missing measured gyro bias");
    }
    if(profile&2) {
      auto unchosen=readback;unchosen.mag_alignment_skipped=0;
      loaded.rebuildFromBlob(Store::sealed_(unchosen));
      check(calibrationUnusableReason(true,loaded)!=nullptr,"missing mag alignment cannot silently select nominal axes");
      auto failed=readback;mat_to_rowmajor9_(rim,failed.R_IM);
      check(!validateBlob(Store::sealed_(failed)),"skipped alignment cannot retain a failed magnetic rotation");
      failed=readback;failed.mag_alignment_ok=1;
      check(!validateBlob(Store::sealed_(failed)),"alignment cannot be both qualified and skipped");
    }
    auto tampered=readback;tampered.gyro_matrix_skipped^=1;
    loaded.rebuildFromBlob(tampered);
    check(!loaded.record_valid && !loaded.gyro_matrix_skipped && !loaded.mag_alignment_skipped,
      "integrity failure clears profile permissions with the coefficients");
  }
  b=rb;b.gyro_scale_only=1;
  mat_to_rowmajor9_(Vector3f(1.04f,.975f,1.02f).asDiagonal().toDenseMatrix(),b.gyro_S);
  b.gyro_coeff_crc=gyroCoeffCrc(b);ImuCalBlobV6 scale_record;
  check(store.saveVerified(b,scale_record),"qualified scale-only model must survive verified persistence");
  RuntimeCals scale_runtime;scale_runtime.rebuildFromBlob(scale_record);
  check(!calibrationUnusableReason(true,scale_runtime) && scale_runtime.gyr.matrix_ok &&
    scale_runtime.gyr.scale_only && !scale_runtime.gyro_matrix_skipped &&
    scale_runtime.gyr.S.diagonal().isApprox(Vector3f(1.04f,.975f,1.02f)),
    "startup must retain measured scales and distinguish them from bias-only/full models");
  auto bad_scale=scale_record;bad_scale.gyro_S[1]=.01f;bad_scale.gyro_coeff_crc=gyroCoeffCrc(bad_scale);
  check(!validateBlob(Store::sealed_(bad_scale)),"scale-only flag cannot qualify unmeasured cross-axis terms");
  bad_scale=scale_record;bad_scale.gyro_matrix_ok=0;bad_scale.gyro_coeff_crc=gyroCoeffCrc(bad_scale);
  check(!validateBlob(Store::sealed_(bad_scale)),"scale-only model must have passed fit qualification");
  bad_scale=scale_record;bad_scale.gyro_scale_only=2;bad_scale.gyro_coeff_crc=gyroCoeffCrc(bad_scale);
  check(!validateBlob(Store::sealed_(bad_scale)),"unknown gyro model flags rejected");
  const Vector3f x(1,2,-9);
  check(runtime.applyAccel(x,25).isApprox(rvh*rhi*runtime.acc.apply(x,25),1e-6f),"accel intrinsic then housing then vessel");
  check(runtime.applyGyro(x,25).isApprox(rvh*rhi*runtime.gyr.apply(x,25),1e-6f),"full gyro matrix applied once");
  check(runtime.applyMag(x).isApprox(rvh*rhi*rim*runtime.mag.apply(x),1e-6f),"mag alignment applied once in correct order");
  struct NoiseConsumer {
    Vector3f density=Vector3f::Zero(),walk=Vector3f::Zero();unsigned mask=0,calls=0;
    void set_gyro_noise_density_rad_sqrt_s(const Vector3f& n){density=n;++calls;}
    void set_measured_gyro_bias_walk(const Vector3f& n,uint8_t m){walk=n;mask=m;++calls;}
  } consumer;
  check(!runtime.configureProcessNoise(consumer) && !consumer.calls,"absent noise cannot configure an estimator");
  b=rb;b.noise.valid=1;b.noise.white_mask=63;b.noise.walk_mask=1<<4;
  b.noise.duration_ms=2400000;b.noise.samples=480001;b.noise.period_s=.005f;
  for(int k=0;k<6;++k){b.noise.white_density[k]=k<3?.002f:.0002f;
    b.noise.sample_std[k]=b.noise.white_density[k]/std::sqrt(b.noise.period_s);}
  b.noise.bias_random_walk[4]=.00003f;
  ImuCalBlobV6 with_noise;
  check(store.saveVerified(b,with_noise),"measured noise round-trips in the unified record");
  runtime.rebuildFromBlob(with_noise);
  check(runtime.configureProcessNoise(consumer) && consumer.calls==2 && consumer.mask==2 &&
    consumer.density.isApprox(Vector3f::Constant(.0002f)) && consumer.walk==Vector3f(0,.00003f,0),
    "runtime passes continuous gyro density and only identifiable walk axes");
  check(runtime.noiseAccelStd().isApprox(Vector3f::Constant(.002f/std::sqrt(.005f))),
    "runtime retains measured accelerometer sample scatter in vessel coordinates");
  b=with_noise;b.noise.duration_ms=2399999;
  check(!validateBlob(Store::sealed_(b)),"a short acquisition cannot masquerade as measured noise");
  runtime.rebuildFromBlob(Store::sealed_(b));
  consumer.calls=0;
  check(!runtime.configureProcessNoise(consumer) && !consumer.calls,"invalid noise record clears earlier runtime parameters");
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
  // Actual wizard pose count and train/validation split, including full PCB rotations.
  atoms3r_magdiag::Pose train[9],validation[4];float rms=NAN;
  atoms3r_magdiag::MagModel model;
  for(float degrees:{0.f,2.f,30.f,90.f,140.f,179.f}) {
    std::mt19937 rng(19);std::normal_distribution<float> n;
    const Matrix3f expected=Eigen::AngleAxisf(degrees/atoms3r_magdiag::kRadToDeg,
      Vector3f(.3f,1,.2f).normalized()).toRotationMatrix();
    for(int i=0;i<13;++i) {
      Eigen::Quaternionf q(n(rng),n(rng),n(rng),n(rng));q.normalize();
      atoms3r_magdiag::Pose pose;pose.a=q.conjugate()*Vector3f(0,0,-9.80256f);
      pose.m_raw=expected.transpose()*(q.conjugate()*Vector3f(20,0,43));
      if(i<9)train[i]=pose;else validation[i-9]=pose;
    }
    const auto alignment=atoms3r_magdiag::solveAlignment(train,9,model,true);
    check(alignment.ok && alignment.R.isApprox(expected,1e-3f) && alignment.R.determinant()>.999f &&
      atoms3r_magdiag::validateAlignment(alignment,validation,4,model,rms),
      "wizard-sized independent magnetic alignment across full rotation range");
    // Validation stays outside the solve and cannot rescue a corrupted endpoint.
    validation[0].m_raw=-validation[0].m_raw;float corrupt_rms=NAN;
    check(!atoms3r_magdiag::validateAlignment(alignment,validation,4,model,corrupt_rms),
      "independent magnetic corruption remains rejected");
  }
  for(auto& p:train)p=train[0];
  const auto degenerate=atoms3r_magdiag::solveAlignment(train,9,model,true);
  check(!degenerate.ok && degenerate.reason==atoms3r_magdiag::AlignmentFail::GEOMETRY,
    "ridge or multiple seeds cannot qualify degenerate alignment");
  std::printf("calibration frames: blob_bytes=%zu independent_alignment_rms=%.8g\n",sizeof(rb),double(rms));
  return errors?1:0;
}
