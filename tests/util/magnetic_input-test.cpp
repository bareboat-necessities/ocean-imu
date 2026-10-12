#define EIGEN_NON_ARDUINO
#include "util/MagneticInput.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cassert>
using namespace ocean_imu::magnetic;
Q turn(float angle){return Q(Eigen::AngleAxisf(angle,V3::UnitZ()));}
struct Receiver {
 unsigned calls=0; Observation received;
 void updateMag(const V3&){++calls;}
 void updateMagTransported(const V3&,const Observation& o){++calls;received=o;}
 bool lastMagUpdateApplied()const{return false;} // reference learning can consume a packet
};
int main(){
 Input input;input.reset(.00135f,M3::Identity());
 atoms3r_ical::MagAcquisition source;
 const V3 rate(0,0,25),earth(20,0,47),acc(0,0,-9.8f);
 input.advance(1000000,rate,0,source,earth,V3::Zero(),false);
 input.capture(Q::Identity(),Q::Identity(),acc,rate);
 source.delay_us=10000;source.timing_allowance_s=.03f;
 source.acquired(1020,1020000,1020000,1020000);
 input.advance(1020000,rate,0,source,turn(.25f).conjugate()*earth,V3::Zero(),false);
 Receiver receiver;
 input.correct(receiver,1020000);assert(receiver.calls==0&&input.pending());
 // A real 0.5-radian rotation across a supported 20-ms history is not a
 // magnetic fault. Reference snapshots must not discard the valid packet.
 input.capture(turn(.5f),turn(.5f),acc,rate);
 input.correct(receiver,1020000);assert(receiver.calls==1&&!input.pending());
 assert(receiver.received.attitude_bw.angularDistance(turn(.25f))<1e-6f);
 input.correct(receiver,1020000);assert(receiver.calls==1);
 input.advance(1040000,rate,0,source,turn(.25f).conjugate()*earth,V3::Zero(),false);
 input.capture(turn(1),turn(1),acc,rate);
 // Even with zero calibrated lag, acquisition and estimator epochs differ
 // when the reader observes a packet between paired estimator frames.
 source.delay_us=0;source.acquired(1050,1050000,1050000,1050000);
 input.advance(1060000,rate,0,source,turn(1.25f).conjugate()*earth,V3::Zero(),false);
 input.capture(turn(1.5f),turn(1.5f),acc,rate);
 input.correct(receiver,1060000);assert(receiver.calls==2);
 assert(receiver.received.attitude_bw.angularDistance(turn(1.25f))<1e-6f);
}
