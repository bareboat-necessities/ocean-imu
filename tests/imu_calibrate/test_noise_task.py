#!/usr/bin/env python3
"""Run the production noise worker with a sensor clock and real capture gates."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()
start = source.index("static void fitTaskNoise_(")
end = source.index("{", start)
depth = 1
while depth:
    end += 1
    depth += (source[end] == "{") - (source[end] == "}")
worker = source[start:end + 1]

cpp = r'''
#define EIGEN_NON_ARDUINO
#define ATOMS3R_NOISE_RAW_LOG 0
#include "imu_calibrate/SensorNoise.h"
#include "imu_calibrate/GyroCalCapture.h"
#include <atomic>
#include <memory>
#include <random>
#include <string>
#include <iostream>
using Vector3f=Eigen::Vector3f;
enum class Mode {COMPLETE,LOST,GAP,FROZEN,TEMPERATURE,CANCEL};
Mode mode;uint64_t clock_us;uint32_t calls=0,lost_us=0;bool stopped=false;
uint32_t micros(){return uint32_t(clock_us);}
uint32_t millis(){return uint32_t(clock_us/1000);}
void delayMicroseconds(uint32_t us){clock_us+=us;}
void vTaskDelay(unsigned){clock_us+=1000;}
void vTaskDelete(void*){}
struct {std::string log;void println(const char* s){log+=s;}
  template<class... A>void printf(const char*,A...){} } Serial;
struct {int Imu=0;} M5;
struct ImuSample {Vector3f a,w,m;float tempC=25;uint32_t sample_us=0;};
bool readImuMapped(int&,ImuSample&);
struct ImuCalCfg {static constexpr float g_cal_local=9.8025605f;};
struct Blob {unsigned token=42;};
struct RuntimeCals {
  void rebuildFromBlob(const Blob&){}
  Vector3f applyAccel(const Vector3f& a,float)const{return a;}
  Vector3f applyGyro(const Vector3f& w,float)const{return w;}
  Vector3f applyMag(const Vector3f& m)const{return m;}
};
const char* calibrationUnusableReason(bool,const RuntimeCals&){return nullptr;}
class ImuCalWizard {
public:
  imu_cal::GyroCalibrator<float,400,8> gyroCal_;
  std::unique_ptr<imu_cal::StationaryNoiseCapture<>> noise_capture_;
  std::atomic<bool> noise_cancel_{false};std::atomic<float> noise_progress_{0};
  struct FitCtx {ImuCalWizard* wiz;const Blob* blob;bool ok=false,done=false;
    imu_cal::FitFail reason=imu_cal::FitFail::OK;};
  static uint32_t hwmBytes_(){return 20000;}
WORKER
};
ImuCalWizard* active=nullptr;
std::mt19937 rng(318);std::normal_distribution<float> gaussian;
bool readImuMapped(int&,ImuSample& s){
  ++calls;s.sample_us=micros();
  const bool acquired=active->noise_capture_->samples()>=100;
  if(acquired && mode==Mode::LOST){
    if(!stopped){lost_us=micros();stopped=true;}
    // Bound a broken worker in the regression itself; cancellation must not
    // be needed to notice that the sensor stopped publishing samples.
    if(uint32_t(micros()-lost_us)>1000000)active->noise_cancel_=true;
    return false;
  }
  if(acquired && mode==Mode::GAP && !stopped){stopped=true;return false;}
  if(acquired && mode==Mode::CANCEL)active->noise_cancel_=true;
  if(acquired && mode==Mode::TEMPERATURE)s.tempC=28;
  s.a=Vector3f(0,0,-ImuCalCfg::g_cal_local);s.w.setZero();
  if(!(acquired && mode==Mode::FROZEN))for(int k=0;k<3;++k){
    s.a[k]+=.028284271f*gaussian(rng);s.w[k]=.0014142136f*gaussian(rng);
  }
  static Vector3f mag(23,1,41);
  if(calls%20==0)for(int k=0;k<3;++k)mag[k]=(k==0?23.f:k==1?1.f:41.f)+.01f*gaussian(rng);
  s.m=mag;return true;
}
void require(bool ok,const char* why){if(!ok){std::cerr<<why<<'\n';std::exit(1);}}
int main(){
  for(Mode test:{Mode::COMPLETE,Mode::LOST,Mode::GAP,Mode::FROZEN,Mode::TEMPERATURE,Mode::CANCEL}){
    mode=test;clock_us=0xffff0000ULL;calls=lost_us=0;stopped=false;Serial.log.clear();
    ImuCalWizard wizard;active=&wizard;
    wizard.noise_capture_=std::make_unique<imu_cal::StationaryNoiseCapture<>>();
    const Blob previous;ImuCalWizard::FitCtx ctx{&wizard,&previous};
    ImuCalWizard::fitTaskNoise_(&ctx);
    imu_cal::SensorNoiseModel model;
    require(ctx.done && previous.token==42,"worker must publish completion without changing calibration");
    if(test==Mode::COMPLETE){
      if(!ctx.ok || !wizard.noise_capture_->report(model))std::cerr<<"status="<<int(wizard.noise_capture_->status())
        <<" samples="<<wizard.noise_capture_->samples()<<" calls="<<calls<<" reason="<<int(ctx.reason)<<'\n';
      require(ctx.ok && wizard.noise_capture_->report(model),"production worker cannot complete forty minutes");
      require(model.duration_ms>=2400000 && model.samples>=480001 && wizard.noise_progress_==1,
              "capture duration, sensor sample count or progress is wrong");
      require(std::fabs(model.white_density[3]/.0001f-1)<.1f,"worker noise units differ from analyzer");
    } else {
      require(!ctx.ok && !wizard.noise_capture_->report(model),"interrupted acquisition must not yield a model");
      if(test==Mode::LOST)require(!wizard.noise_cancel_ && uint32_t(micros()-lost_us)<=25000 &&
        Serial.log.find("sample stream stopped")!=std::string::npos,"lost sensor must fail promptly without cancellation");
      if(test==Mode::GAP)require(wizard.noise_capture_->status()==imu_cal::NoiseCaptureStatus::TIMING,
        "missing frame must remain a timing failure");
      if(test==Mode::FROZEN)require(wizard.noise_capture_->status()==imu_cal::NoiseCaptureStatus::FROZEN,
        "frozen registers must remain invalid");
      if(test==Mode::TEMPERATURE)require(wizard.noise_capture_->status()==imu_cal::NoiseCaptureStatus::TEMPERATURE,
        "thermal excursion must remain invalid");
      if(test==Mode::CANCEL)require(wizard.noise_cancel_,"operator cancellation was lost");
    }
  }
  std::cout<<"noise worker: 40 minutes, clock wrap, sample loss/gap/freeze, temperature and cancellation PASS\n";
}
'''.replace("WORKER", worker)

with tempfile.TemporaryDirectory() as directory:
    cpp_file = Path(directory) / "noise-task.cpp"
    cpp_file.write_text(cpp)
    eigen = root / "third_party/eigen"
    if not (eigen / "Eigen/Dense").exists():
        eigen = Path("/usr/include/eigen3")
    binary = Path(directory) / "noise-task"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-I" + str(root / "src"), "-I" + str(eigen),
        str(cpp_file), "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True, timeout=30)
