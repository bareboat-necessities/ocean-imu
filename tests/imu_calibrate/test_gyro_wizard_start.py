#!/usr/bin/env python3
"""Replay the production turn controller and the real per-turn retry wrapper."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()

def block(marker):
    start = source.index(marker)
    end = brace = source.index("{", start)
    depth = 1
    while depth:
        end += 1
        depth += (source[end] == "{") - (source[end] == "}")
    return source[start:end + 1]

# The capture worker must stay free of display/input work. The independent
# scheduler replay below additionally runs the real worker and UI on threads.
worker = block("static void gyroTurnTask_(")
assert "ui_." not in worker and "Input::" not in worker

cpp = r'''
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/GyroTurnCapture.h"
#include <iostream>
#include <random>
#include <cstdlib>
using V=Eigen::Vector3f;
using Phase=imu_cal::GyroTurnPhase;
using Fit=imu_cal::GyroRotationFit<131072,24>;
void require(bool b,const char* text){if(!b){std::cerr<<text<<'\n';std::exit(1);}}
int main(){
  const V targets[4]={V(0,0,-1),V(-1,0,0),V(0,1,0),V(0,-1,0)};
  for(int mode=0;mode<8;++mode) {
    auto fit=std::make_unique<Fit>();
    imu_cal::GyroCalibrator<float,400,8> cal;cal.g=9.8025605f;cal.max_gyro_norm=.12f;cal.max_accel_dev=.8f;
    imu_cal::GyroCalibration<float> bias;bias.ok=bias.biasT.ok=true;
    std::mt19937 rng(811+mode);std::normal_distribution<float> noise(0,1);
    uint64_t clock=0xffff0000ULL;V base(0,0,-cal.g),field(24,5,40);
    bool gapped=false,saw_progress=false,saw_restart=false;int retries=0;
    for(int segment=0;segment<20;++segment) {
      bool accepted=false;
      while(!accepted) {
        // Button handling can change the device orientation before acquisition.
        Eigen::AngleAxisf pickup(.3f,V(.2f,-.4f,1).normalized());base=pickup*base;field=pickup*field;
        imu_cal::GyroTurnCapture<> capture(*fit,cal,bias,targets,segment);capture.begin(uint32_t(clock/1000));
        uint64_t turn=0,twist=0;float total=0,theta=0,yaw=0;const float speed=mode==5?.15f:.6f;V axis=V::UnitX(),a=base,m=field;
        V mag_noise=V::Zero();bool jolted=false,wrong_pose=false;
        for(int sample=0;sample<40000;++sample) {
          clock+=5000;
          if(capture.phase()==Phase::TWIST && !twist)twist=clock+(mode==0?60000000:1500000);
          if(capture.phase()==Phase::TURN && !turn) {
            turn=clock+500000;
            const V goal=targets[capture.target()];
            total=std::acos(std::clamp(base.normalized().dot(goal),-1.f,1.f));
            axis=goal.cross(base.normalized()).normalized();
          }
          const float ye=twist && clock>twist?float(clock-twist)*1e-6f:0;
          const float previous_yaw=yaw;
          if(!turn)yaw=(segment%2?-1.f:1.f)*std::min(1.3f,ye*speed);
          const float elapsed=turn && clock>turn?float(clock-turn)*1e-6f:0;
          float rate=elapsed>0?speed:0;
          if(mode==2 && elapsed>=.25f && elapsed<6.5f){theta=.15f;rate=0;}
          else theta=std::min(total,mode==2 && elapsed>=6.5f?.15f+(elapsed-6.5f)*speed:elapsed*speed);
          if(theta>=total)rate=0;
          if(mode==3 && segment==6 && theta>.2f && !gapped){clock+=50000;gapped=true;}
          const Eigen::AngleAxisf rotation(-theta,axis);
          a=rotation*base;V w=axis*rate;
          if(yaw!=previous_yaw)w+=base.normalized()*(segment%2?-speed:speed);
          const float eps=std::sin(float(clock%1000000)*.000071f);
          a.x()+=.001f*eps;w+=V(.0001f*eps,.00008f*eps,-.00006f*eps);
          if(mode==1){a+=.01f*V(noise(rng),noise(rng),noise(rng));w+=.0023f*V(noise(rng),noise(rng),noise(rng));}
          if(clock%100000<5000){
            const float me=turn && clock>turn?float((clock-turn)/100000)*.1f:0;
            const float mt=std::min(total,mode==2 && me>=.25f && me<6.5f?.15f:
              mode==2 && me>=6.5f?.15f+(me-6.5f)*speed:me*speed);
            m=Eigen::AngleAxisf(-mt,axis)*(Eigen::AngleAxisf(-yaw,base.normalized())*field);
            if(mode==1)mag_noise=.3f*V(noise(rng),noise(rng),noise(rng));
            m+=mag_noise+V(.02f*eps,0,0);
          }
          if(mode==4 && capture.phase()==Phase::HOLD && capture.progress()>.2f && !jolted) {
            a.x()+=1.5f;jolted=true;
          }
          if(mode==6 && capture.phase()==Phase::HOLD && capture.progress()>.8f && !jolted) {
            a=Eigen::AngleAxisf(.9f,a.normalized().unitOrthogonal())*a;jolted=true;
          }
          if(mode==7 && segment==1 && retries==0 && capture.phase()==Phase::HOLD && capture.progress()>.25f)wrong_pose=true;
          if(wrong_pose) {a=-targets[capture.target()]*cal.g+V(.001f*eps,0,0);w=V(.0001f*eps,0,0);}
          const auto phase=capture.update(uint32_t(clock/1000),uint32_t(clock),&a,&w,25,&m);
          if(!twist && capture.progress()>.01f)saw_progress=true;
          if((mode==4 || mode==6) && jolted && capture.resets() && capture.progress()<.1f){
            require(phase==Phase::HOLD,"a movement restarts the hold, not the turn instruction");saw_restart=true;
          }
          if(mode==6 && jolted)require(phase==Phase::HOLD || phase==Phase::DONE,
              "late USB-up jitter must not silently revert the finish screen to TIP");
          // Missing polls must not make the finish hold flicker back to TIP.
          if(phase==Phase::HOLD)require(capture.update(uint32_t(clock/1000),0,nullptr,nullptr,NAN,nullptr)==Phase::HOLD,
              "between sensor deliveries the hold instruction must remain stable");
          if(phase==Phase::FAILED){
            require(retries==0 && ((mode==3 && fit->segments()==6 && std::string(capture.reason())=="Rotation sample gap") ||
              (mode==7 && fit->segments()==1 && std::string(capture.reason())=="Target pose changed")),"only the interrupted turn may be discarded");
            ++retries;base=a;field=m;break;
          }
          if(phase==Phase::DONE){accepted=true;base=capture.accel();field=capture.mag();break;}
        }
        require(accepted || ((mode==3 || mode==7) && retries==1),"guided turn did not complete");
      }
    }
    require(fit->segments()==20 && fit->samples()<131072 && saw_progress,"twenty turns must retain their progress within storage");
    if(mode==4)require(saw_restart,"endpoint movement must exercise an explained hold restart");
    imu_cal::AccelCalibration<float> acc;acc.ok=acc.biasT.ok=true;acc.g=cal.g;
    const bool ok=fit->fit(acc,bias);
    std::cerr<<"mode="<<mode<<" fit="<<ok<<" reason="<<imu_cal::gyroRotationReason(fit->reason)
      <<" matrix_error="<<(bias.S-Eigen::Matrix3f::Identity()).norm()<<" samples="<<fit->samples()<<'\n';
    require(ok && (bias.S-Eigen::Matrix3f::Identity()).norm()<.005f,
      "guided faces must identify the gyro matrix and pass independent validation");
  }
  {
    auto fit=std::make_unique<Fit>();imu_cal::GyroCalibrator<float,400,8> cal;
    cal.max_gyro_norm=.12f;cal.max_accel_dev=.8f;
    imu_cal::GyroCalibration<float> bias;bias.ok=bias.biasT.ok=true;
    imu_cal::GyroTurnCapture<> capture(*fit,cal,bias,targets,0);capture.begin(0);
    uint32_t motion=0;Phase phase=Phase::REST;
    for(uint32_t ms=5;ms<45000 && phase!=Phase::FAILED;ms+=5) {
      if(phase==Phase::TWIST && !motion)motion=ms;
      const float eps=std::sin(ms*.071f),angle=motion?(ms-motion)*.00015f:0;
      V a(.001f*eps,0,-cal.g),w(.0001f*eps,0,motion?-.15f:0);
      V m=Eigen::AngleAxisf(angle,V::UnitZ())*V(24,5,40);
      m.x()+=.02f*std::sin(float(ms/100)*.7f);
      phase=capture.update(ms,ms*1000,&a,&w,25,&m);
    }
    std::cerr<<"budget phase="<<int(phase)<<" reason="<<capture.reason()<<" samples="<<fit->samples()<<'\n';
    require(phase==Phase::FAILED && std::string(capture.reason())=="Turn took too long" &&
      fit->samples()==0 && fit->segments()==0,
      "a long live motion must release only its own budget, without a false stale-sensor failure");
  }
  std::cout<<"gyro guide: named faces/wait/noise/pause/gap/restart/storage/matrix validation PASS\n";
}
'''
# Include memory explicitly for compilers where Eigen doesn't import it.
cpp = cpp.replace('#include <cstdlib>', '#include <cstdlib>\n#include <memory>')
eigen = root / "third_party/eigen"
if not (eigen / "Eigen/Dense").exists():
    eigen = Path("/usr/include/eigen3")
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "gyro-guide"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-fno-finite-math-only", "-Wall", "-Wextra", "-Werror",
        "-I", str(root / "src"), "-isystem", str(eigen), "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=60)
