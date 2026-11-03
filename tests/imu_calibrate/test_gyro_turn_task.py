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
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
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
struct M5Ui{enum class GyroFitAction {RETRY_FIT,ADD_TURNS,BIAS_ONLY,ABORT};};
struct M5UiCfg{static constexpr uint8_t ROT_READ=0;};
struct ImuCalWizardCfg{static constexpr uint32_t FIT_STACK_BYTES=32768;};
struct ImuSample {V a,w,m;float tempC=25;uint32_t sample_us=0;};
std::thread worker;const auto ui_thread=std::this_thread::get_id();
std::atomic<uint64_t> clock_us{0xffff0000ULL};std::atomic<int> requested{-1};std::atomic<bool> twist_requested{false};
std::atomic<bool> manual_retry{false};bool retry_on_hold=false;
bool allocation_failure=false;int samples=0,redraws=0,max_gap=0,core=-1;unsigned bytes=0;
uint32_t micros(){return uint32_t(clock_us.load());}uint32_t millis(){return uint32_t(clock_us.load()/1000);}
void delay(unsigned){std::this_thread::sleep_for(std::chrono::microseconds(100));}
void vTaskDelay(unsigned){clock_us+=1000;std::this_thread::sleep_for(std::chrono::microseconds(30));}
void vTaskDelete(void*){}
struct Input{static void update(){}static bool tapPressed(){return manual_retry.exchange(false);}};
struct {void println(const char*){}template<class... A>void printf(const char*,A...){} } Serial;
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
  int turn_total_=imu_cal::GyroCalibrationPlan::INITIAL;
  static constexpr uint8_t turn_pose_ids_[4]={0,2,4,5};
CONTEXT turn_;
  struct Accel {bool expectedUp(int pose,V& up){
    for(int i=0;i<4;++i)if(turn_pose_ids_[i]==pose){up=targets[i];return true;}return false;
  }} accel_;
  struct UI {
    uint8_t rot=0;bool finishing=false,saw_finish2=false;
    uint8_t rotation()const{return rot;}
    void setRotation(uint8_t r){rot=r;}
    void waitTapHere(const char*,const char*,const char*){}
    void title(const char* s){finishing=false;require(std::strlen(s)<=10,"headline does not fit");if(std::string(s)=="TWIST")twist_requested=true;}
    void line(const char* s){if(std::string(s).starts_with("Finish "))finishing=true;if(std::string(s)=="Finish 2/12")saw_finish2=true;require(std::strlen(s)<=21,"instruction does not fit 128px display");
      const char* captions[]={"Screen facing up","USB end up","Left edge down","Right edge down"};
      require(std::string(s).find("hands off")==std::string::npos,"hand support must be permitted");
      for(int i=0;i<4;++i)if(std::string(s)==captions[i])requested=i;
    }
    void lineAt(int y,const char* s){require(y<=100 && std::strlen(s)<=21,"dynamic row overlaps or wraps");}
    void bar01(float p){if(retry_on_hold && finishing && p>.25f){manual_retry=true;retry_on_hold=false;}}
    void showOkAuto(const char*,const char*,const char* = "OK"){}
  } ui_;
  std::unique_ptr<imu_cal::GyroRotationFit<131072,24>> gyro_rotations_{new imu_cal::GyroRotationFit<131072,24>};
  imu_cal::GyroCalibrator<float,400,8> gyroCal_;
  imu_cal::GyroCalibration<float> gyr_out_;
  V endpoint_a_,endpoint_w_,endpoint_m_;float endpoint_temp_=25;
  struct Physics {uint64_t last=0,start=0,twist=0;float total=0,yaw=0;V axis=V::UnitX();} physics_;
  bool readSample_(ImuSample& s){
    require(std::this_thread::get_id()!=ui_thread,"sensor read ran on display task");
    auto& last=physics_.last;auto& start=physics_.start;auto& twist=physics_.twist;
    auto& total=physics_.total;auto& yaw=physics_.yaw;auto& axis=physics_.axis;
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
ROTATION
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
  struct Fit {imu_cal::GyroRotationFail reason=imu_cal::GyroRotationFail::GEOMETRY;imu_cal::GyroRotationCfg cfg;int n=0;void clear(){n=0;}int segments()const{return n;}};
  std::unique_ptr<Fit> gyro_rotations_{new Fit};
  enum class FitKind {GYRO,GYRO_MATRIX,MAG_ALIGNMENT};
  struct {bool ok=false;imu_cal::FitFail reason=imu_cal::FitFail::OK;const char* detail="Mag poses too alike";} fit_;
  bool skip_before=false,skip_after=false,skip_turn=false,mag_bad=true,missing_memory=false;
  bool gyro_matrix_skipped_=false,mag_alignment_skipped_=false;
  imu_cal::GyroCalibration<float> gyro_bias_fit_;atoms3r_magdiag::Alignment alignment_;
  bool always_bad=false,initial_bias_fail=false;int matrix_fits=0,adds=0,intros=0,mag_fits=0,mag_captures=0;
  int captures=0,fits=0,calls=0,retries=0,alignment_pose_count_=0;
  int turn_total_=imu_cal::GyroCalibrationPlan::INITIAL;
  struct UI {
    RetryWizard* owner;
    void fail(const char*,const char*){}
    void waitTapHere(const char*,const char*,const char*){++owner->intros;}
    void showOkAuto(const char*,const char*,const char*){}
    M5Ui::GyroFitAction gyroTiltMenu(){return owner->skip_before?M5Ui::GyroFitAction::BIAS_ONLY:M5Ui::GyroFitAction::RETRY_FIT;}
    M5Ui::GyroFitAction gyroFitFailMenu(const char*,bool can_add=false){
      if(owner->skip_after)return M5Ui::GyroFitAction::BIAS_ONLY;
      if(can_add){++owner->adds;return M5Ui::GyroFitAction::ADD_TURNS;}
      return M5Ui::GyroFitAction::ABORT;
    }
    M5Ui::GyroFitAction gyroTurnFailMenu(const char* why,const char* kept){
      require(std::string(why)=="gyro scatter" && std::string(kept)=="6 turns kept",
        "timeout must explain the cause and retained progress");
      ++owner->retries;return owner->skip_turn?M5Ui::GyroFitAction::BIAS_ONLY:M5Ui::GyroFitAction::RETRY_FIT;
    }
    bool retryMenuHere(const char* title,const char* why,const char* kept=nullptr){
      if(std::string(title)=="Retry gyro bias")return true;
      if(std::string(title)=="Add 4 turns?"){++owner->adds;return true;}
      require(std::string(title)=="Retry turn 7 / 12" && std::string(why)=="gyro scatter" &&
        std::string(kept)=="6 turns kept","timeout must explain the cause and retained progress");
      ++owner->retries;return true;
    }
  } ui_{this};
  imu_cal::GyroCalibration<float> gyr_out_;
  V endpoint_a_=V(0,0,-9.8f),endpoint_w_=V::Zero(),endpoint_m_=V(24,5,40);
  float endpoint_temp_=25;
  struct Pose {V a,m_raw;} raw_alignment_[24];float alignment_temp_[24]{};
  bool prepareGyroRotationStorage_(){return !missing_memory;}
  bool captureGyro_(){++captures;return true;}
  bool runFitTask_(FitKind kind,const char*,bool){
    if(kind==FitKind::MAG_ALIGNMENT){++mag_fits;fit_.ok=!mag_bad;return fit_.ok;}
    if(kind==FitKind::GYRO){++fits;if(initial_bias_fail && fits==1){fit_.reason=imu_cal::FitFail::GYRO_NOT_STILL;return false;}return true;}
    ++matrix_fits;fit_.reason=imu_cal::FitFail::DEGENERATE_COVERAGE;
    fit_.ok=!always_bad && gyro_rotations_->segments()>=16;
    gyr_out_.S(0,0)=fit_.ok?1.03f:1.2f;gyr_out_.rotation_condition=999;return fit_.ok;
  }
  bool fakeTurn_(int segment,const char*& why){
    ++calls;if(segment==6 && retries==0){why="gyro scatter";return false;}
    ++gyro_rotations_->n;return true;
  }
  RetryWizard(){gyr_out_.ok=gyr_out_.biasT.ok=true;gyr_out_.biasT.b0=V(.01f,.02f,.03f);}
BASIC
COMPLETE_MAG
COMPLETE
WRAPPER
};
int main(){
  RetryWizard retry;
  require(retry.captureGyroRotations_() && retry.captures==1 && retry.fits==1 &&
    retry.calls==13 && retry.retries==1 && retry.gyro_rotations_->segments()==12 &&
    retry.alignment_pose_count_==13,"retry wrapper repeated completed turns or fitted bias");
  require(retry.captureGyroRotations_(16) && retry.captures==1 && retry.fits==1 &&
    retry.calls==17 && retry.gyro_rotations_->segments()==16 && retry.alignment_pose_count_==17,
    "extra turns must append without discarding captures, refitting bias or repeating instructions");
  require(retry.intros==3,"introductory instructions must not repeat when turns are appended");
  RetryWizard supported;
  const V strict_bias=supported.gyr_out_.biasT.b0;
  supported.endpoint_w_=strict_bias+V(0,0,.004f); // below the endpoint mean-rate limit
  require(supported.captureGyroRotations_(),"supported endpoint capture failed");
  const float yaw_error_45s=(supported.gyr_out_.biasT.b0-strict_bias).z()*45*57.2957795f;
  std::cout<<"supported-endpoint bias: induced_yaw_45s_deg="<<yaw_error_45s<<'\n';
  require(supported.gyr_out_.biasT.b0==strict_bias,
    "motion-tolerant endpoints must never replace the strict stationary bias");

  RetryWizard bias_retry;bias_retry.initial_bias_fail=true;
  require(bias_retry.captureGyroRotations_() && bias_retry.captures==2 && bias_retry.intros==3 &&
    bias_retry.gyro_rotations_->segments()==12,"bias-fit retry must preserve earlier stages and not repeat instructions");
  RetryWizard complete;
  require(complete.completeGyroCalibration_() && complete.captures==1 && complete.calls==17 &&
    complete.adds==1 && complete.matrix_fits==2 && complete.alignment_pose_count_==17,
    "failed final fit must add data to the existing set, not restart it");
  RetryWizard rejected;rejected.always_bad=true;
  require(!rejected.completeGyroCalibration_() && rejected.captures==1 && rejected.calls==21 &&
    rejected.gyro_rotations_->segments()==20 && rejected.alignment_pose_count_==21,
    "final rejection must keep all captures until explicit abort; no automatic reset");
  for(int phase:{0,1,2,3}) {
    RetryWizard basic;basic.skip_before=phase==0;basic.skip_after=phase==1 || phase==3;basic.skip_turn=phase==2;basic.missing_memory=phase==3;
    basic.always_bad=true;
    require(basic.completeGyroCalibration_() && basic.gyro_matrix_skipped_ && basic.gyr_out_.ok &&
      !basic.gyr_out_.matrix_ok && basic.gyr_out_.S==Eigen::Matrix3f::Identity() &&
      basic.gyr_out_.biasT.b0==V(.01f,.02f,.03f) && basic.gyr_out_.rotation_condition==0 &&
      basic.gyr_out_.rotation_train_segments==0,"skip must preserve this session's strict bias and discard unqualified matrix metadata");
    require(basic.gyro_rotations_->segments()==(phase==0 || phase==3?0:phase==1?12:6),
      "skipping tilt must not erase completed rotations or require extra turns");
    require(basic.matrix_fits==(phase==1?1:0),"skipping before capture/on timeout must avoid matrix fitting");
  }
  for(bool bad:{false,true}) {
    RetryWizard magnetic;magnetic.mag_bad=bad;magnetic.gyr_out_.scale_only=true;
    magnetic.gyro_rotations_->n=20;magnetic.alignment_pose_count_=21;
    magnetic.gyr_out_.S(0,0)=1.03f;magnetic.alignment_.ok=true;
    require(magnetic.completeMagAlignment_() && magnetic.mag_fits==1 &&
      magnetic.mag_captures==0 && magnetic.matrix_fits==0 && magnetic.captures==0 &&
      magnetic.gyro_rotations_->segments()==20 && magnetic.gyr_out_.S(0,0)==1.03f &&
      magnetic.gyr_out_.scale_only && magnetic.mag_alignment_skipped_==bad &&
      magnetic.alignment_.ok==!bad,
      "magnetic failure must skip automatically, retain gyro scales and continue without recapture");
  }
  ImuCalWizard w;w.gyroCal_.g=9.8025605f;w.gyroCal_.max_gyro_norm=.12f;w.gyroCal_.max_accel_dev=.8f;
  w.gyr_out_.ok=w.gyr_out_.biasT.ok=true;const char* why=nullptr;
  require(w.runGyroTurn_(1,why),"slow display interrupted real turn capture");worker.join();
  require(samples>1000 && redraws>3 && max_gap<=20000 && w.gyro_rotations_->segments()==1,
    "80ms LCD redraw must preserve every sensor epoch and the accepted turn");
  require(w.ui_.rotation()==2 && w.ui_.saw_finish2,"USB-up capture must use the accelerometer's upside-down correction");
  require(w.gyroReadRotation_(targets[0]*9.8f,targets,2)==0 &&
    w.gyroReadRotation_(targets[1]*9.8f,targets,0)==2 &&
    w.gyroReadRotation_(targets[2]*9.8f,targets,0)==1 &&
    w.gyroReadRotation_(targets[3]*9.8f,targets,0)==3,
    "all gyro pose rotations must match the existing accelerometer guidance");
  require(core==0 && bytes>=32768,"capture worker core/stack settings incorrect");
  w.physics_={};requested=-1;twist_requested=false;retry_on_hold=true;
  require(!w.runGyroTurn_(1,why) && std::string(why)=="Retry requested","tap during a hold must offer a retry");worker.join();
  require(w.gyro_rotations_->segments()==1,"manual retry must preserve the completed turn");
  allocation_failure=true;require(!w.runGyroTurn_(1,why) && std::string(why)=="Capture task failed",
    "task allocation must fail without discarding accepted turns");
  require(w.gyro_rotations_->segments()==1,"failed worker start lost an accepted turn");
  std::cout<<"gyro worker: 80ms redraws, 200Hz capture, fitting storage, task failure and screen text PASS\n";
}
'''
cpp = (cpp.replace("CONTEXT", block("struct TurnCtx"))
       .replace("BASIC", block("bool useGyroBiasOnly_("))
       .replace("COMPLETE_MAG", block("bool completeMagAlignment_("))
       .replace("ROTATION", block("static int gyroReadRotation_("))
       .replace("HELPER", block("static const char* gyroRestHelp_("))
       .replace("RENDER", block("void drawGyroTurn_("))
       .replace("CALLBACK", block("static void gyroTurnTask_("))
       .replace("COMPLETE", block("bool completeGyroCalibration_("))
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
    # Compile the full production UI header, not a copied menu implementation.
    # Its driver stubs only supply the LCD and clock unavailable on the host.
    (Path(directory) / "Arduino.h").write_text(r'''
#pragma once
#include <cstdint>
inline uint32_t now_ms=0;
inline uint32_t millis(){return now_ms;}
inline void delay(unsigned ms){now_ms+=ms;}
''')
    (Path(directory) / "M5Unified.h").write_text(r'''
#pragma once
#include "Arduino.h"
#include <vector>
#include <string>
#include <cstring>
#include <cassert>
constexpr int TFT_BLACK=0,TFT_WHITE=1,TFT_DARKGREY=2,TFT_GREEN=3;
struct DisplayMock {
  int rotation=0,y=0,size=1;std::vector<std::string> text;
  void setRotation(int r){rotation=r;}
  void setBrightness(int){} void setTextSize(int s){size=s;} void setTextColor(int,int){}
  void fillScreen(int){} void setCursor(int,int row){y=row;}
  void println(const char* s = ""){assert(std::strlen(s)*6*size<=128);assert(y+8*size<=128);text.emplace_back(s);y+=8*size;}
  void print(const char*){}
  void drawRect(int,int,int,int,int){} void fillRect(int,int,int,int,int){}
  void wakeup(){} int width(){return 128;}int height(){return 128;}int getCursorY(){return y;}
};
struct ButtonMock {bool pressed=false;uint32_t press_after=0;std::vector<uint32_t> scheduled;
  bool wasPressed(){bool p=pressed || (press_after && millis()>=press_after);
    if(!scheduled.empty() && millis()>=scheduled.front()){p=true;scheduled.erase(scheduled.begin());}
    if(p)press_after=0;pressed=false;return p;}};
inline struct {DisplayMock Display;ButtonMock BtnA;void update(){}} M5;
''')
    menus = r'''
#define EIGEN_NON_ARDUINO
#include "AtomS3R/AtomS3R_M5Ui.h"
#include "imu_calibrate/FrameAlignment.h"
#include <cassert>
#include <iostream>
int main(){
  atoms3r_ical::M5Ui ui;ui.begin();
  ui.setRotation(2);const uint32_t read_start=millis();M5.BtnA.press_after=read_start+300000;
  ui.waitTapHere("READ","Take your time");assert(millis()-read_start>=300000 && ui.rotation()==2);
  M5.BtnA.pressed=true;ui.waitTapHere("GYRO","Turn 2/12");
  assert(M5.Display.rotation==2 && ui.rotation()==2);
  M5.BtnA.pressed=true;assert(ui.retryMenuHere("Retry turn 2","Pose changed","1 turn kept"));
  assert(M5.Display.rotation==2 && ui.rotation()==2);
  ui.showOkAuto("Turn kept");assert(M5.Display.rotation==2);
  M5.BtnA.pressed=true;ui.waitTap("OTHER");assert(M5.Display.rotation==0);
  ui.setRotation(1);M5.BtnA.pressed=true;assert(ui.retryMenu("Other stage","Retry"));
  assert(M5.Display.rotation==0);
  uint32_t t=millis();M5.BtnA.scheduled={t+10,t+110};
  assert(ui.gyroFitFailMenu("Gyro check failed")==atoms3r_ical::M5Ui::GyroFitAction::BIAS_ONLY);
  t=millis();M5.BtnA.scheduled={t+10,t+110,t+210,t+1100};
  assert(ui.gyroFitFailMenu("Gyro check failed")==atoms3r_ical::M5Ui::GyroFitAction::RETRY_FIT);
  t=millis();M5.BtnA.scheduled={t+10,t+110,t+210,t+1100,t+1200,t+1300};
  assert(ui.gyroFitFailMenu("Gyro check failed")==atoms3r_ical::M5Ui::GyroFitAction::ABORT);
  t=millis();M5.BtnA.scheduled={t+10,t+110};
  assert(ui.gyroTiltMenu()==atoms3r_ical::M5Ui::GyroFitAction::BIAS_ONLY);
  M5.BtnA.pressed=true;
  assert(ui.gyroFitFailMenu("Need varied tilt axes",true)==atoms3r_ical::M5Ui::GyroFitAction::ADD_TURNS);
  t=millis();M5.BtnA.scheduled={t+10,t+110};
  assert(ui.gyroTurnFailMenu("gyro scatter","6 turns kept")==atoms3r_ical::M5Ui::GyroFitAction::BIAS_ONLY);
  auto gestures=[&](const std::vector<int>& groups){
    M5.BtnA.pressed=false;M5.BtnA.scheduled.clear();uint32_t when=millis()+10;
    for(int n:groups){for(int k=0;k<n;++k)M5.BtnA.scheduled.push_back(when+100*k);when+=1000;}
  };
  gestures({1});auto axes=ui.mountingMenu();
  assert(axes.up==imu_cal::FrameAxis::NEG_Z && axes.forward==imu_cal::FrameAxis::POS_X);
  for(int up=0;up<6;++up)for(int fwd=0;fwd<6;++fwd) {
    const auto u=imu_cal::FrameAxis(up),f=imu_cal::FrameAxis(fwd);
    if(!imu_cal::frameAxesPerpendicular(u,f))continue;
    std::vector<int> groups={2,1}; // change mount, acknowledge edge naming
    for(int n=0;n<up;++n)groups.push_back(1);
    groups.push_back(2);
    int current=imu_cal::frameAxesPerpendicular(u,imu_cal::FrameAxis::POS_X)?3:0;
    while(current!=fwd){groups.push_back(1);do{current=(current+1)%6;}while(current/2==up/2);}
    groups.push_back(2);groups.push_back(1); // choose forward, confirm mount
    gestures(groups);axes=ui.mountingMenu();assert(axes.up==u && axes.forward==f);
    Eigen::Matrix3f R;assert(imu_cal::mountingRotation(axes.up,axes.forward,R));
  }
  gestures({2,1,3,1});axes=ui.mountingMenu(); // back from up chooser, use default
  assert(axes.up==imu_cal::FrameAxis::NEG_Z);
  gestures({2,1,2,2,3,1});axes=ui.mountingMenu(); // cancel confirmation, use default
  assert(axes.forward==imu_cal::FrameAxis::POS_X);
  gestures({3,1,1});axes=ui.mountingMenu(); // untimed help returns without losing fits
  gestures({2,1});assert(ui.startMenu(true)==atoms3r_ical::M5Ui::StartAction::GYRO_BIAS);
  gestures({2,2});assert(ui.startMenu(true)==atoms3r_ical::M5Ui::StartAction::NOISE);
  std::cout<<"mounting: default, 24 orientations, cancel/back/help and saved-fit tools PASS\n";
  std::cout<<"gyro menus: real UI preserves pose rotation; other stages keep default behavior PASS\n";
}
'''
    menu_binary = Path(directory) / "gyro-menus"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-I", directory, "-I", str(root / "src"), "-isystem", str(eigen), "-x", "c++", "-",
        "-o", str(menu_binary),
    ], input=menus, text=True, check=True)
    subprocess.run([str(menu_binary)], check=True, timeout=10)
