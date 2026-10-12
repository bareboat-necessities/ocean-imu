#define EIGEN_NON_ARDUINO
#include "util/AttitudeProtection.h"
#include "util/MagneticInput.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cassert>
using namespace ocean_imu::attitude;
int main(){
 Protection ins{MotionUse::INS},compass;
 const Vec earth(20,0,47);bool compass_held=false;
 for(unsigned i=1;i<=24000;++i){
  const float t=float(i)*.01f;const Vec force(0,0,-9.8f-3*std::sin(2.513274f*t));
  const auto a=ins.step(.01f,Vec::Zero(),Vec::Zero(),force,true,9.8f,Quat::Identity(),earth,true,i%4==0,1000000+i*10000);
  const auto b=compass.step(.01f,Vec::Zero(),Vec::Zero(),force,true,9.8f,Quat::Identity(),earth,true,i%4==0,1000000+i*10000);
  assert(!a.dynamic&&!a.hold_attitude);if(i>20000)compass_held|=b.hold_attitude;
 }
 assert(compass_held);ins.reset();
 // Raw noncommuting gyro rotations survive the force window.
 GyroWindow window;const Vec x(1,0,0),y(0,2,0);window.push(x,.005f);window.push(y,.005f);
 Vec gyro;float dt;window.take(gyro,dt);
 assert(std::abs(dt-.01f)<1e-7f);
 assert(rotation(Vec(gyro*dt)).angularDistance(rotation(Vec(x*.005f))*rotation(Vec(y*.005f)))<1e-6f);
 window.take(gyro,dt);assert(!gyro.allFinite()&&!std::isfinite(dt));
 // A net zero window must not conceal a real high-rate hazard.
 ins.observeGyro(.005f,Vec(0,0,3));ins.observeGyro(.005f,Vec(0,0,-3));
 auto d=ins.step(.01f,Vec::Zero(),Vec::Zero(),Vec(0,0,-9.8f),true,9.8f,Quat::Identity(),earth,true,true,2000000);
 assert(d.dynamic);
 for(unsigned i=1;i<=350;++i){
  const Vec force(0,0,-9.8f-3*std::sin(float(i)*.02513274f));
  d=ins.step(.01f,Vec::Zero(),Vec::Zero(),force,true,9.8f,Quat::Identity(),earth,true,i%4==0,2000000+i*10000);
 }
 assert(!d.dynamic&&!d.hold_attitude);
 d=ins.step(.03f,Vec::Zero(),Vec::Zero(),Vec(0,0,-9.8f),true,9.8f,Quat::Identity(),earth,true,true,5530000);
 assert(d.dynamic);

 // Packets between force frames retain their own epoch, including zero lag.
 ocean_imu::magnetic::Input magnetic;magnetic.reset(0,Eigen::Matrix3f::Identity());
 atoms3r_ical::MagAcquisition timing;
 magnetic.advance(1000000,Vec::Zero(),0,timing,earth,Vec::Zero(),false);
 magnetic.capture(Quat::Identity(),Quat::Identity(),Vec(0,0,-9.8f),Vec::Zero());
 timing.acquired(1005,1005000,1005000,1005000);
 magnetic.advance(1005000,Vec::Zero(),0,timing,earth,Vec::Zero(),false);
 assert(magnetic.takeNewSource());assert(!magnetic.takeNewSource());
 magnetic.advance(1010000,Vec::Zero(),0,timing,earth,Vec::Zero(),false);
 const Quat end=rotation(Vec(0,0,.1f));magnetic.capture(end,end,Vec(0,0,-9.8f),Vec::Zero());
 assert(magnetic.observation().attitude_bw.angularDistance(rotation(Vec(0,0,.05f)))<1e-6f);
 assert(!magnetic.takeNewSource());
}
