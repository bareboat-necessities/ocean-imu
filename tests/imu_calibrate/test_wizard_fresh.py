#!/usr/bin/env python3
"""Exercise full-calibration admission with an existing current-schema record."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()
# Execute the production entry through begin(), stopping before physical
# acquisition. Hardware and menus are stubs; coefficient selection is real.
entry = source[source.index("bool runAndSave_("):source.index("      AccelIo io(*this);")]
entry += "return true; }}"
cpp = r'''
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/AccelCalFit.h"
#include <cstdlib>
#include <iostream>
using Vector3f=Eigen::Vector3f;
struct ImuCalBlobV6 {bool gyro_ok=true;float gyro_b0[3]={.03f,-.02f,.01f};unsigned token=47;};
struct RuntimeCals {void rebuildFromBlob(const ImuCalBlobV6&){};};
const char* calibrationUnusableReason(bool,const RuntimeCals&){return nullptr;}
imu_cal::AccelThermalPrior accelThermalPriorFrom(const ImuCalBlobV6&,unsigned,unsigned,unsigned){
  imu_cal::AccelThermalPrior p;p.valid=true;p.k[0]=.003;return p;
}
struct M5Ui {enum class StartAction {FULL,NOISE,CANCEL};};
struct {void println(const char*){}template<class... A>void printf(const char*,A...){} } Serial;
void clearM5UnifiedImuCalibration(){}
void require(bool b,const char* why){if(!b){std::cerr<<why<<'\n';std::exit(1);}}
class ImuCalWizard {
public:
  struct Buffer {void clear(){}} gyroCal_,magCal_;
  struct Store {bool load(ImuCalBlobV6& out){out=ImuCalBlobV6{};return true;}} store_;
  struct UI {
    M5Ui::StartAction action=M5Ui::StartAction::FULL;bool noise_available=false;
    M5Ui::StartAction startMenu(bool noise){noise_available=noise;return action;}
  } ui_;
  struct Config {double g=9.8025605;} accel_fcfg_;
  struct Accel {
    int begins=0;bool used_prior=true,used_reference=true;Vector3f reference;
    void begin(int,const Config&,const imu_cal::AccelThermalPrior& prior,
               const Vector3f& gyro,bool known){++begins;used_prior=prior.valid;reference=gyro;used_reference=known;}
  } accel_;
  int accel_ccfg_=0;
  unsigned sensor_id_lo_=12,sensor_id_hi_=34,imu_type_=5;
  int alignment_pose_count_=0,timing_fit_=0;
  bool mag_verified_=false,timing_saved_=false,noise_called=false;
  imu_cal::GyroCalibration<float> gyr_out_;
  imu_cal::MagCalibration<float> mag_out_;
  void configureCalibrators_(){}void sensorIdentity_(){}
  bool prepareGyroRotationStorage_(){return true;}
  bool captureNoise_(const ImuCalBlobV6& previous,ImuCalBlobV6&){noise_called=previous.token==47;return true;}
ENTRY
};
int main(){
  ImuCalWizard wizard;ImuCalBlobV6 out;
  require(wizard.runAndSave_(out) && wizard.accel_.begins==1,"full calibration did not reach acquisition");
  require(!wizard.accel_.used_prior && !wizard.accel_.used_reference && wizard.accel_.reference.isZero(),
          "fresh calibration must not substitute a previous slope or stationary gyro reference");
  require(wizard.ui_.noise_available,"a valid record must still admit separate noise acquisition");
  wizard.ui_.action=M5Ui::StartAction::NOISE;
  require(wizard.runAndSave_(out) && wizard.noise_called && wizard.accel_.begins==1,
          "noise acquisition must keep its calibrated frame without restarting intrinsic calibration");
  wizard.ui_.action=M5Ui::StartAction::CANCEL;
  require(!wizard.runAndSave_(out) && wizard.accel_.begins==1,"cancel must not begin acquisition");
  std::cout<<"wizard admission: fresh coefficients, separate calibrated noise capture and cancellation PASS\n";
}
'''.replace("ENTRY", entry)
with tempfile.TemporaryDirectory() as directory:
    cpp_file = Path(directory) / "fresh-calibration.cpp"
    cpp_file.write_text(cpp)
    eigen = root / "third_party/eigen"
    if not (eigen / "Eigen/Dense").exists():
        eigen = Path("/usr/include/eigen3")
    binary = Path(directory) / "fresh-calibration"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-I" + str(root / "src"), "-I" + str(eigen),
        str(cpp_file), "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True)
