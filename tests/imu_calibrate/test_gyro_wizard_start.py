#!/usr/bin/env python3
"""Replay real rotation capture, with hardware/UI replaced by a virtual clock."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()
start = source.index("bool captureGyroRotations_(")
brace = source.index("{", start)
depth, end = 1, brace + 1
while depth:
    depth += (source[end] == "{") - (source[end] == "}")
    end += 1
method = source[start:end]

cpp = r'''
#define EIGEN_NON_ARDUINO
#define INPUT 0x01
#include "imu_calibrate/GyroCalFit.h"
#include "imu_calibrate/GyroCalCapture.h"
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <iostream>
using Eigen::Vector3f;
using RotationFit=imu_cal::GyroRotationFit<32768,24>;
struct ImuSample {Vector3f a,w,m;float tempC=25;uint32_t sample_us=0;};
struct Simulation {
  uint64_t us=0xffff0000ULL,last_sample=us,turn=0;
  uint32_t idle_ms=12000;
  bool latency=false,inject_gap=false,idle_gap=false,gapped=false,stale=false,small_turn=false;
  int captures=0,retries=0,starts=0;
  Eigen::Vector3d a{0,0,-9.8025605},m{24,5,40},axis{1,0,0};
  std::string reason;
  double angle() const {return turn && us>turn+uint64_t(idle_ms)*1000 ?
    std::min(small_turn?.05:.9,(us-turn-uint64_t(idle_ms)*1000)*.6e-6):0;}
  void next() {
    const Eigen::AngleAxisd r(-angle(),axis);a=r*a;m=r*m;turn=0;++starts;
    // Moving the device to press the button is unobserved. Its old endpoint
    // must not be reused as the starting gravity direction of this turn.
    const Eigen::AngleAxisd pickup(.3,Eigen::Vector3d(.2,-.4,1).normalized());a=pickup*a;m=pickup*m;
    axis=Eigen::Vector3d(std::sin(starts*.73),std::cos(starts*1.21),std::sin(starts*1.77+.5)).normalized();
    if(std::fabs(axis.dot(a.normalized()))>.85)axis=a.normalized().unitOrthogonal();
    us+=120000; // menu/display work must precede the sensor epoch
  }
} sim;
uint32_t micros(){return uint32_t(sim.us);}
uint32_t millis(){return uint32_t(sim.us/1000);}
void delay(unsigned ms){sim.us+=uint64_t(ms)*1000;}
struct Input {static void update(){}};
struct SerialMock {
  template<class... A> void printf(const char*,A...){}
} Serial;
template<class S> bool tryCalLogLine(S&,const char*){return true;}
void require(bool b,const char* text){if(!b){std::cerr<<text<<" reason="<<sim.reason<<'\n';std::exit(1);}}
class ImuCalWizard {
public:
  struct UI {
    void waitTap(const char*,const char*,const char*){sim.next();}
    void title(const char*){}
    void line(const char*){}
    void lineAt(int,const char* s){if(std::string(s)=="Turn now" && !sim.turn)sim.turn=sim.us;}
    int cursorY(){return 0;}
    void bar01(float){}
    void fail(const char*,const char*){}
    bool retryMenu(const char*,const char* why,const char*){++sim.retries;sim.reason=why;return false;}
  } ui_;
  enum class FitKind{GYRO};
  std::unique_ptr<RotationFit> gyro_rotations_{new RotationFit};
  imu_cal::GyroCalibrator<float,400,8> gyroCal_;
  imu_cal::GyroCalibration<float> gyr_out_;
  Vector3f endpoint_a_,endpoint_w_=Vector3f::Zero(),endpoint_m_;
  float endpoint_temp_=25;
  int alignment_pose_count_=0;
  struct Pose {Vector3f a,m_raw;} raw_alignment_[24];float alignment_temp_[24]{};
  bool prepareGyroRotationStorage_(){return true;}
  bool runFitTask_(FitKind,const char*,bool){return true;}
  bool captureGyro_(){
    ++sim.captures;endpoint_a_=sim.a.cast<float>();endpoint_m_=sim.m.cast<float>();
    gyroCal_.g=9.8025605f;gyr_out_.ok=gyr_out_.biasT.ok=true;return true;
  }
  bool readSample_(ImuSample& s){
    if(sim.stale)return false;
    if(sim.latency){sim.us+=55000;sim.latency=false;}
    if(sim.inject_gap && !sim.gapped && sim.angle()>.1){sim.us+=50000;sim.gapped=true;}
    if(sim.idle_gap && !sim.gapped && sim.turn && sim.us-sim.turn>500000){sim.us+=50000;sim.gapped=true;}
    if(sim.us-sim.last_sample<5000)return false;
    sim.last_sample=sim.us;s.sample_us=micros();
    const Eigen::AngleAxisd r(-sim.angle(),sim.axis);
    const double eps=std::sin(double(sim.us%1000000)*.000071);
    s.a=(r*sim.a).cast<float>()+Vector3f(float(.001*eps),0,0);
    s.w=Vector3f(float(.0001*eps),float(.00008*eps),float(-.00006*eps));
    if(sim.angle()>0 && sim.angle()<(sim.small_turn?.05:.9))s.w+=sim.axis.cast<float>()*.6f;
    // Independent 10 Hz magnetic observations, as on the real sensor.
    const double mag_angle=sim.turn && sim.us>sim.turn+uint64_t(sim.idle_ms)*1000 ?
      std::min(sim.small_turn?.05:.9,double((sim.us-sim.turn-uint64_t(sim.idle_ms)*1000)/100000)*.06):0;
    s.m=(Eigen::AngleAxisd(-mag_angle,sim.axis)*sim.m).cast<float>();
    s.m.x()+=float(.02*std::sin(double(sim.us/100000)));
    return true;
  }
METHOD
};
int main(){
  for(unsigned wait:{12000u,60000u}) {
    sim=Simulation{};sim.idle_ms=wait;sim.latency=true;
    ImuCalWizard w;
    require(w.captureGyroRotations_(),"waiting before motion must not produce RETRY");
    require(sim.retries==0 && sim.captures==1 && w.gyro_rotations_->segments()==20,
            "all twenty real captures must complete without restarting the bias stage");
    require(w.gyro_rotations_->samples()<32768,"stationary waiting must not exhaust trajectory RAM");
    imu_cal::AccelCalibration<float> acc;acc.ok=acc.biasT.ok=true;acc.g=9.8025605f;
    const bool fitted=w.gyro_rotations_->fit(acc,w.gyr_out_);
    std::cerr<<"wait="<<wait<<" fit="<<fitted<<" reason="<<imu_cal::gyroRotationReason(w.gyro_rotations_->reason)
      <<" error="<<(w.gyr_out_.S-Eigen::Matrix3f::Identity()).norm()<<" validation="<<w.gyr_out_.rotation_validation_rms
      <<" samples="<<w.gyro_rotations_->samples()<<'\n';
    require(fitted &&
      (w.gyr_out_.S-Eigen::Matrix3f::Identity()).norm()<5e-5f && w.gyr_out_.rotation_validation_rms<2e-5f,
      "waiting and unobserved menu motions must not bias the fitted matrix/held-out endpoints");
  }
  sim=Simulation{};sim.idle_gap=true;
  {ImuCalWizard w;require(w.captureGyroRotations_() && sim.retries==0 && sim.gapped,
      "a pre-motion gap must requalify the start without losing an active rotation");}
  sim=Simulation{};sim.inject_gap=true;
  {ImuCalWizard w;require(!w.captureGyroRotations_() && sim.retries==1 && sim.gapped &&
      sim.reason=="Rotation sample gap","a gap during actual motion must still reject the segment");}
  sim=Simulation{};sim.stale=true;
  {ImuCalWizard w;require(!w.captureGyroRotations_() && sim.retries==1 &&
      sim.reason=="No fresh samples","a stalled start must have its real failure reason");}
  sim=Simulation{};sim.small_turn=true;
  {ImuCalWizard w;require(!w.captureGyroRotations_() && sim.retries==1 &&
      sim.reason=="Rotate to another pose","an insufficient turn must not report OK");}
  std::cout<<"gyro wizard start: initial latency/wait/wrap/storage/motion gap/stale/angle PASS\n";
}
'''.replace("METHOD", method)
eigen = root / "third_party/eigen"
if not (eigen / "Eigen/Dense").exists():
    eigen = Path("/usr/include/eigen3")
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "gyro-wizard-start"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-fno-finite-math-only", "-Wall", "-Wextra", "-Werror",
        "-I", str(root / "src"), "-isystem", str(eigen),
        "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=60)
