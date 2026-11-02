#!/usr/bin/env python3
"""Run the production capture task concurrently with deliberately slow LCD work."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()

def block(marker):
    start = source.index(marker)
    end = source.index("{", start)
    depth = 1
    while depth:
        end += 1
        depth += (source[end] == "{") - (source[end] == "}")
    return source[start:end + 1]

cpp = r'''
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/GyroTurnCapture.h"
#include "imu_calibrate/AccelCalCapture.h"
#include <atomic>
#include <thread>
#include <chrono>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <cstring>
using V=Eigen::Vector3f;using Vector3f=V;
using Phase=imu_cal::GyroTurnPhase;
using TaskHandle_t=void*;constexpr int pdPASS=1;
struct ImuCalWizardCfg{static constexpr uint32_t FIT_STACK_BYTES=32768;};
struct ImuSample {V a,w,m;float tempC=25;uint32_t sample_us=0;};
std::thread worker;const auto ui_thread=std::this_thread::get_id();
std::atomic<uint64_t> clock_us{0xffff0000ULL};std::atomic<int> requested{-1};std::atomic<bool> twist_requested{false};
bool allocation_failure=false;int samples=0,redraws=0,max_gap=0,core=-1;unsigned bytes=0;
uint32_t micros(){return uint32_t(clock_us.load());}uint32_t millis(){return uint32_t(clock_us.load()/1000);}
void delay(unsigned){std::this_thread::sleep_for(std::chrono::microseconds(100));}
void vTaskDelay(unsigned){clock_us+=1000;std::this_thread::sleep_for(std::chrono::microseconds(30));}
void vTaskDelete(void*){}
struct Input{static void update(){}};
struct {template<class... A>void printf(const char*,A...){} } Serial;
template<class S>bool tryCalLogLine(S&,const char*){return true;}
void require(bool b,const char* s){if(!b){std::cerr<<s<<'\n';std::exit(1);}}
int xTaskCreatePinnedToCore(void(*fn)(void*),const char*,uint32_t stack,void* ctx,
                           unsigned,TaskHandle_t* task,int on_core) {
  bytes=stack;core=on_core;
  if(allocation_failure)return 0;
  *task=reinterpret_cast<void*>(1);worker=std::thread([=]{fn(ctx);});return pdPASS;
}
const V targets[4]={V(0,0,-1),V(-1,0,0),V(0,1,0),V(0,-1,0)};
class ImuCalWizard {
public:
  static constexpr uint8_t turn_pose_ids_[4]={0,2,4,5};
CONTEXT turn_;
  struct Accel {bool expectedUp(int pose,V& up){
    for(int i=0;i<4;++i)if(turn_pose_ids_[i]==pose){up=targets[i];return true;}return false;
  }} accel_;
  struct UI {
    void waitTap(const char*,const char*,const char*){}
    void title(const char* s){require(std::strlen(s)<=10,"headline does not fit");if(std::string(s)=="TWIST")twist_requested=true;}
    void line(const char* s){require(std::strlen(s)<=21,"instruction does not fit 128px display");
      for(int i=0;i<4;++i)if(std::string(s)==imu_cal::kAccelPoses[turn_pose_ids_[i]].name)requested=i;
    }
    void lineAt(int y,const char* s){require(y<=94 && std::strlen(s)<=21,"dynamic row overlaps or wraps");}
    void bar01(float){}
    void showOkAuto(const char*,const char*){}
  } ui_;
  std::unique_ptr<imu_cal::GyroRotationFit<131072,24>> gyro_rotations_{new imu_cal::GyroRotationFit<131072,24>};
  imu_cal::GyroCalibrator<float,400,8> gyroCal_;
  imu_cal::GyroCalibration<float> gyr_out_;
  V endpoint_a_,endpoint_w_,endpoint_m_;float endpoint_temp_=25;
  bool readSample_(ImuSample& s){
    require(std::this_thread::get_id()!=ui_thread,"sensor read ran on display task");
    static uint64_t last=0,start=0,twist=0;static float total=0,yaw=0;static V axis=V::UnitX();
    const uint64_t now=clock_us.load();
    if(samples)max_gap=std::max(max_gap,int(now-last));last=now;++samples;
    if(twist_requested && !twist)twist=now;
    const float previous_yaw=yaw;if(twist && !start)yaw=std::min(1.3f,float(now-twist)*.6e-6f);
    if(requested>=0 && !start){start=now;const V goal=targets[requested];
      total=std::acos(V(0,0,-1).dot(goal));axis=goal.cross(V(0,0,-1)).normalized();}
    const float angle=start?std::min(total,float(now-start)*.6e-6f):0;
    const float rate=start && angle<total?.6f:0;
    const float eps=std::sin(float(now%1000000)*.000071f);
    s.sample_us=micros();s.a=Eigen::AngleAxisf(-angle,axis)*V(0,0,-9.8025605f)+V(.001f*eps,0,0);
    s.w=axis*rate+V(.0001f*eps,.00008f*eps,-.00006f*eps);
    if(yaw>previous_yaw)s.w+=V(0,0,-.6f);
    const float ma=start?std::min(total,float((now-start)/100000)*.06f):0;
    s.m=Eigen::AngleAxisf(-ma,axis)*(Eigen::AngleAxisf(-yaw,V(0,0,-1))*V(24,5,40))+V(.02f*std::sin(float(now/100000)),0,0);
    return true;
  }
HELPER
RENDER
  void slowDraw_(int segment,Phase phase,int target,int progress,const char* reason,bool fresh) {
    require(std::this_thread::get_id()==ui_thread,"display ran on capture worker");
    // Real render instructions are exercised, then an 80ms SPI stall is injected.
    drawGyroTurn_(segment,phase,target,progress,reason,fresh);++redraws;
    std::this_thread::sleep_for(std::chrono::milliseconds(80));
  }
CALLBACK
RUN
};
class RetryWizard {
public:
  struct Fit {int n=0;void clear(){n=0;}int segments()const{return n;}};
  std::unique_ptr<Fit> gyro_rotations_{new Fit};
  enum class FitKind {GYRO};
  int captures=0,fits=0,calls=0,retries=0,alignment_pose_count_=0;
  struct UI {
    RetryWizard* owner;
    void fail(const char*,const char*){}
    bool retryMenu(const char* title,const char* why,const char* kept){
      require(std::string(title)=="Retry turn 7 / 20" && std::string(why)=="gyro scatter" &&
        std::string(kept)=="6 turns kept","timeout must explain the cause and retained progress");
      ++owner->retries;return true;
    }
  } ui_{this};
  imu_cal::GyroCalibration<float> gyr_out_;
  V endpoint_a_=V(0,0,-9.8f),endpoint_w_=V::Zero(),endpoint_m_=V(24,5,40);
  float endpoint_temp_=25;
  struct Pose {V a,m_raw;} raw_alignment_[24];float alignment_temp_[24]{};
  bool prepareGyroRotationStorage_(){return true;}
  bool captureGyro_(){++captures;return true;}
  bool runFitTask_(FitKind,const char*,bool){++fits;return true;}
  bool fakeTurn_(int segment,const char*& why){
    ++calls;if(segment==6 && retries==0){why="gyro scatter";return false;}
    ++gyro_rotations_->n;return true;
  }
WRAPPER
};
int main(){
  RetryWizard retry;
  require(retry.captureGyroRotations_() && retry.captures==1 && retry.fits==1 &&
    retry.calls==21 && retry.retries==1 && retry.gyro_rotations_->segments()==20 &&
    retry.alignment_pose_count_==21,"retry wrapper repeated completed turns or fitted bias");
  ImuCalWizard w;w.gyroCal_.g=9.8025605f;w.gyroCal_.max_gyro_norm=.12f;w.gyroCal_.max_accel_dev=.8f;
  w.gyr_out_.ok=w.gyr_out_.biasT.ok=true;const char* why=nullptr;
  require(w.runGyroTurn_(0,why),"slow display interrupted real turn capture");worker.join();
  require(samples>1000 && redraws>3 && max_gap<=20000 && w.gyro_rotations_->segments()==1,
    "80ms LCD redraw must preserve every sensor epoch and the accepted turn");
  require(core==0 && bytes>=32768,"capture worker core/stack settings incorrect");
  allocation_failure=true;require(!w.runGyroTurn_(1,why) && std::string(why)=="Capture task failed",
    "task allocation must fail without discarding accepted turns");
  require(w.gyro_rotations_->segments()==1,"failed worker start lost an accepted turn");
  std::cout<<"gyro worker: 80ms redraws, 200Hz capture, fitting storage, task failure and screen text PASS\n";
}
'''
cpp = (cpp.replace("CONTEXT", block("struct TurnCtx"))
       .replace("HELPER", block("static const char* gyroRestHelp_("))
       .replace("RENDER", block("void drawGyroTurn_("))
       .replace("CALLBACK", block("static void gyroTurnTask_("))
       .replace("WRAPPER", block("bool captureGyroRotations_(").replace("runGyroTurn_(", "fakeTurn_("))
       .replace("RUN", block("bool runGyroTurn_(").replace("drawGyroTurn_(", "slowDraw_(")))
eigen = root / "third_party/eigen"
if not (eigen / "Eigen/Dense").exists():
    eigen = Path("/usr/include/eigen3")
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "gyro-turn-task"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-pthread", "-fno-finite-math-only", "-Wall", "-Wextra", "-Werror",
        "-I", str(root / "src"), "-isystem", str(eigen), "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=30)
