#!/usr/bin/env python3
"""Run the production gyro-only save path with real capture, fit and storage."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()
repair = source[source.index("  bool repairGyroBias_("):source.index("  bool captureMounting_()")]
# Reuse the independent current-record fixture and memory-backed KV adapter.
fixtures = (root / "tests/imu_calibrate/calibration_frames-test.cpp").read_text()
fixtures = fixtures[:fixtures.index("int main(){")]
cpp = fixtures + r'''
#include "imu_calibrate/GyroCalCapture.h"
using ImuCalStoreNvs=Store;
class ImuCalWizard {
public:
  Store store_;
  imu_cal::GyroCalibrator<float,400,8> gyroCal_;
  imu_cal::GyroCalibration<float> gyr_out_;
  bool timing_saved_=true,done=false,cancel=false,fail_fit=false;
  unsigned captures=0,fits=0;
  enum class FitKind {GYRO};
  struct {imu_cal::FitFail reason=imu_cal::FitFail::OK;} fit_;
  struct UI {
    unsigned retries=0;
    void waitTapHere(const char*,const char*,const char*){}
    void fail(const char*,const char*){}
    bool retryMenuHere(const char*,const char*,const char*){++retries;return false;}
  } ui_;
  bool captureGyro_() {
    ++captures;if(cancel)return false;
    imu_cal::GyroCapture<float,400,8> capture(gyroCal_);capture.begin(0);
    for(uint32_t t=5;t<=12000;t+=5) {
      const float f=float(t)*.001f;
      const Vector3f a(.0001f*std::sin(13*f),0,-9.80256f);
      const Vector3f w=Vector3f(.012f,-.006f,.004f)+
        Vector3f(.0001f*std::sin(11*f),.0001f*std::cos(17*f),.0001f*std::sin(19*f));
      const Vector3f m(20+.01f*std::sin(float(t/100)),0,43);
      if(capture.update(t,&a,&w,26,&m)==imu_cal::GyroCaptureStatus::READY)return true;
    }
    return false;
  }
  bool runFitTask_(FitKind,const char*,bool) {
    ++fits;
    if(fail_fit){fit_.reason=imu_cal::FitFail::GYRO_NOT_STILL;return false;}
    return gyroCal_.fit(gyr_out_,&fit_.reason);
  }
  void showDone_(const ImuCalBlobV6&){done=true;}
REPAIR
};
int main(){
  auto previous=record();
  mat_to_rowmajor9_(Eigen::AngleAxisf(.04f,Vector3f::UnitY()).toRotationMatrix(),previous.R_IM);
  mat_to_rowmajor9_(Eigen::AngleAxisf(.2f,Vector3f::UnitX()).toRotationMatrix(),previous.R_HI);
  mat_to_rowmajor9_(Eigen::AngleAxisf(-.3f,Vector3f::UnitZ()).toRotationMatrix(),previous.R_VH);
  previous.noise.valid=1;previous.noise.white_mask=63;
  previous.noise.duration_ms=2400000;previous.noise.samples=480001;previous.noise.period_s=.005f;
  for(int k=0;k<6;++k){previous.noise.white_density[k]=k<3?.002f:.0002f;
    previous.noise.sample_std[k]=previous.noise.white_density[k]/std::sqrt(previous.noise.period_s);}
  ImuCalWizard wizard;ImuCalBlobV6 verified;
  check(wizard.store_.saveVerified(previous,verified),"repair fixture must be usable");previous=verified;
  ImuCalBlobV6 result;
  check(wizard.repairGyroBias_(previous,result) && wizard.done && !wizard.timing_saved_ &&
    wizard.captures==1 && wizard.fits==1,"gyro-only repair must capture, fit, verify and finish once");
  check(wizard.gyr_out_.biasT.b0.isApprox(Vector3f(.012f,-.006f,.004f),.001f),
    "repair must use fresh strict stationary measurements");
  auto expected=previous;fillGyroFromFit(expected,wizard.gyr_out_,true);
  expected.noise={};expected.mag_timing={};expected=Store::sealed_(expected);
  check(std::memcmp(&result,&expected,sizeof(result))==0,
    "repair must preserve every accel, mag, mounting, identity and unrelated metadata byte");
  RuntimeCals runtime;runtime.rebuildFromBlob(result);
  check(!calibrationUnusableReason(true,runtime) && runtime.gyro_matrix_skipped &&
    !runtime.gyr.matrix_ok && !runtime.gyr.scale_only && runtime.gyr.S.isIdentity() &&
    result.gyro_train_segments==0 && result.gyro_validation_segments==0 && !result.noise.valid,
    "repaired record must boot with measured bias and explicit nominal scale/axes");
  ImuCalBlobV6 loaded;check(wizard.store_.load(loaded) && std::memcmp(&loaded,&result,sizeof(result))==0,
    "repair must persist exactly the verified output");
  for(bool cancel:{true,false}) {
    ImuCalWizard stopped;ImuCalBlobV6 old;
    check(stopped.store_.saveVerified(previous,old),"cancellation fixture must save");
    stopped.cancel=cancel;stopped.fail_fit=!cancel;
    auto output=old;
    check(!stopped.repairGyroBias_(old,output) && !stopped.done &&
      std::memcmp(&output,&old,sizeof(old))==0 && stopped.store_.load(loaded) &&
      std::memcmp(&loaded,&old,sizeof(old))==0,
      "cancelled or failed gyro repair must not publish or overwrite any saved calibration");
  }
  std::printf("gyro repair: strict capture, preserved fits/mounts, verified save and cancellation PASS\n");
  return errors?1:0;
}
'''.replace("REPAIR", repair)
with tempfile.TemporaryDirectory() as directory:
    cpp_file = Path(directory) / "gyro-repair.cpp"
    cpp_file.write_text(cpp)
    eigen = root / "third_party/eigen"
    if not (eigen / "Eigen/Dense").exists():
        eigen = Path("/usr/include/eigen3")
    binary = Path(directory) / "gyro-repair"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-I" + str(root / "src"), "-isystem", str(eigen),
        str(cpp_file), "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True)
