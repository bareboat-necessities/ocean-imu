// Copyright 2026, Mikhail Grushinskiy
#define EIGEN_NON_ARDUINO
// Arduino's GPIO macro must coexist with the host-tested calibration headers.
#define INPUT 0x01
#include "imu_calibrate/MagCalCapture.h"
#include "imu_calibrate/GyroCalCapture.h"
#include "imu_calibrate/AccelCalCapture.h"
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "AtomS3R/AtomS3R_CalLog.h"
#include <cstdio>
#include <cstring>
#include <random>
#include <string>

using V=Eigen::Vector3f;
using M=Eigen::Matrix3f;
using MC=imu_cal::MagCalibrator<float,400>;
using MS=imu_cal::MagCaptureStatus;
using GS=imu_cal::GyroCaptureStatus;
static int checks=0,failures=0;
static void check(bool ok,const char* name) {++checks;if(!ok){++failures;std::fprintf(stderr,"FAIL: %s\n",name);}}
static V direction(float i) {
  const float z=.95f*std::sin(.071f*i),p=.19f*i;
  return V(std::sqrt(1-z*z)*std::cos(p),std::sqrt(1-z*z)*std::sin(p),z);
}
static M distortion() {M d;d<<1.2,.18,.09,.18,.9,.12,.09,.12,1.05;return d;}
static const V offset(8,-5,3);

static void testMagneticQuality() {
  for(int mode=0;mode<5;++mode) {
    MC cal;imu_cal::MagCapture<float,400> cap(cal);cap.begin(0);
    std::mt19937 rng(9421);std::normal_distribution<float> noise(0,1);
    MS status=MS::CAPTURING;
    for(int i=0;i<=563;++i) {
      V m=50*distortion()*direction(i)+offset;
      if(mode==1)m+=.25f*V(noise(rng),noise(rng),noise(rng));
      if(mode==2 && i>=280)m+=V(30,0,0);
      if(mode==3)m=V(40,0,20)+20*V(noise(rng),noise(rng),noise(rng));
      if(mode==4 && i%20==0)m+=V(12,-9,7);
      status=cap.update(uint32_t(80*i),&m);
    }
    check(status==MS::READY,"full rotations complete capture");
    imu_cal::MagCalibration<float> out;
    const bool ok=cal.fit(out);
    if(mode==2 || mode==3) {
      check(!ok && !out.ok,"changing field/noise cloud cannot become a calibration");
      std::printf("mag rejected mode=%d gate=%s\n",mode,imu_cal::magFitGateText(cal.quality.gate));
      continue;
    }
    check(ok,"clean/noisy/outlier magnetic captures qualify");
    check(out.rms<=out.quality.rms+1e-7 && std::fabs(out.rms-out.quality.trimmed_rms)<1e-6,
          "legacy trimmed report is distinct from full-inlier acceptance RMS");
    check(out.quality.iterations>0 && out.quality.refined_cost<=out.quality.initial_cost,"geometric refinement decreases robust objective");
    if(!ok) {std::printf("mag unexpected gate=%s inliers=%d rms=%g cells=%d\n",imu_cal::magFitGateText(cal.quality.gate),cal.quality.inliers,cal.quality.rms,cal.quality.cells);continue;}
    float heading=0;
    float angle_error=0;
    for(int i=0;i<360;++i) {
      const float p=float(M_PI)*i/180;V u(std::cos(p),std::sin(p),0);
      const V v=out.apply(50*distortion()*u+offset);
      heading=std::max(heading,std::fabs(std::remainder(std::atan2(v.y(),v.x())-p,2*float(M_PI)))*180/float(M_PI));
      const V tilted=direction(i+137),corrected=out.apply(50*distortion()*tilted+offset).normalized();
      angle_error=std::max(angle_error,std::atan2(corrected.cross(tilted).norm(),corrected.dot(tilted))*180/float(M_PI));
    }
    check(heading<(mode==1?1.0f:.25f),"calibration retains heading accuracy on unseen clean rotations");
    check(angle_error<(mode==1?1.0f:.25f),"calibration retains vector direction at appreciable pitch and roll");
    // A later sweep checks frozen coefficients; it never enters a refit.
    const M saved_A=out.A;const V saved_b=out.b;
    imu_cal::MagCaptureCfg cfg;cfg.min_time_ms=12000;cfg.timeout_ms=60000;cfg.required_samples=140;
    imu_cal::MagCapture<float,400> verify(cal,cfg,&out);verify.begin(100000);
    for(int i=0;i<=170;++i) {V m=50*distortion()*direction(i+765)+offset;status=verify.update(100000+80*i,&m);}
    check(status==MS::READY,"independent sweep obtains coverage");
    imu_cal::MagFitQuality q;
    check(cal.geometric.check(cal.buf.v,cal.buf.n,out.field_uT,out.A,out.b,q,cal.sample_ms),"independent float validation passes clean data");
    for(int i=0;i<cal.buf.n;++i)cal.buf.v[i]+=V(10,0,0);
    check(!cal.geometric.check(cal.buf.v,cal.buf.n,out.field_uT,out.A,out.b,q,cal.sample_ms),"field change between fit and check rejects candidate");
    check(out.A==saved_A && out.b==saved_b,"verification never modifies candidate coefficients");
    std::printf("mag mode=%d heading_max=%g rms=%g refine=%d\n",mode,heading,out.rms,out.quality.iterations);
  }
  // Uneven time spent on an easy direction must not swamp the remaining 3D data.
  MC cal;
  for(int i=0;i<400;++i) {V m=50*distortion()*(i<180?V(1,0,0):direction(i))+offset;cal.addSample(m);}
  imu_cal::MagCalibration<float> out;
  check(cal.fit(out),"uneven directions remain identifiable with balanced fitting");
  if(out.ok)check((out.b-offset).norm()<.01f,"balanced fit retains offset accuracy");
  // Stored coefficients are checked at runtime too, without new blob layout.
  atoms3r_ical::ImuCalBlobV4 b;b.mag_ok=1;b.mag_field_uT=50;b.mag_A[0]=b.mag_A[4]=b.mag_A[8]=1;
  check(atoms3r_ical::magSetValid(b),"ordinary stored magnetic coefficients validate");
  b.mag_rms=3.22f;check(!atoms3r_ical::magSetValid(b),"poor stored fit RMS cannot bypass acceptance");
  b.mag_rms=0;b.mag_A[8]=.001f;check(!atoms3r_ical::magSetValid(b),"implausible stored matrix cannot bypass acceptance");

  // A short field disturbance concentrated in time cannot hide in the trim.
  cal.clear();
  for(int i=0;i<400;++i)cal.addSample((i<18?54.f:50.f)*direction(i),80*i);
  imu_cal::MagFitQuality q;
  check(!cal.geometric.check(cal.buf.v,400,50.f,M::Identity().eval(),V::Zero().eval(),q,cal.sample_ms) &&
        q.gate==imu_cal::MagFitGate::FIELD_CHANGED,"temporally concentrated interference fails despite acceptable global fraction");

  // Smooth, slow sweeps with differing starting attitudes, noise and isolated
  // outliers exercise the workflow without relying on exact orientations.
  for(int seed=0;seed<20;++seed) {
    cal.clear();imu_cal::MagCapture<float,400> capture(cal);capture.begin(0);
    std::mt19937 rng(471+seed);std::normal_distribution<float> noise(0,1);
    MS status=MS::CAPTURING;
    for(int i=0;i<1000 && status!=MS::READY;++i) {
      V m=50*distortion()*direction(.2f*i+seed*14.3f)+offset+.2f*V(noise(rng),noise(rng),noise(rng));
      if(i%67==0)m+=V(8,-6,4);
      status=capture.update(80*i,&m);
    }
    check(status==MS::READY,"slow arbitrary sweep reaches complete coverage");
    check(cal.fit(out) && (out.b-offset).norm()<.3f,"multi-seed slow/noisy sweep yields accurate offsets");
    if(!out.ok)std::printf("slow sweep seed=%d reason=%s gate=%s rms=%g drift=%g\n",seed,
        imu_cal::fitFailStr(cal.lastFail()),imu_cal::magFitGateText(cal.quality.gate),cal.quality.rms,cal.quality.time_drift);
  }
}

static void testGuidance() {
  const char* wanted[3]={"Tilt forward/back","Roll left/right","Flip over"};
  for(int missing=0;missing<3;++missing) {
    MC cal;imu_cal::MagCapture<float,400> capture(cal);capture.begin(0);
    MS status=MS::CAPTURING;
    for(int i=0;i<=563;++i) {
      V v=V::Zero();v[(missing+1)%3]=50*std::cos(.19f*i);v[(missing+2)%3]=50*std::sin(.19f*i);
      status=capture.update(80*i,&v);
    }
    check(status==MS::COVERAGE_LOW,"planar capture requires more directions");
    check(std::strcmp(capture.hint(45040),wanted[missing])==0,"prompt follows missing coverage");
    for(int i=564;i<850;++i) {V v=50*direction(i);status=capture.update(80*i,&v);}
    check(status==MS::READY,"late 3D motion repairs coverage without restarting capture");
    check(std::strcmp(capture.hint(68000),"Turn slowly")==0,"full coverage restores broad turn guidance");
    imu_cal::MagCalibration<float> out;check(cal.fit(out),"late retained directions also pass fit information gates");
  }
}

static void testGyroCapture() {
  for(int mode=0;mode<10;++mode) {
    imu_cal::GyroCalibrator<float,400> cal;cal.max_gyro_norm=.12f;cal.max_accel_dev=.8f;
    imu_cal::GyroCapture<float,400> capture(cal);
    const uint32_t start=mode==5?0xfffff000u:0u;capture.begin(start);
    std::mt19937 rng(900+mode);std::normal_distribution<float> noise(0,1);
    const V bias(.002f,-.003f,.001f);
    GS status=GS::SETTLING;float angle=0;int ready_ms=-1;
    for(int ms=0;ms<25000;ms+=5) {
      float rate=0;
      if(mode==1 && ms<12000)rate=.04f; // steady yaw invisible to gravity
      if(mode==2 && ms<550)rate=.09f;  // original half-moving capture
      if(mode==3 && ms>=3000 && ms<6000)rate=.05f; // interruption after progress
      if(mode==4 && ms<12000)rate=.04f; // roll detected without magnetometer
      angle+=rate*.005f;
      V a=V(0,0,-cal.g)+.01f*V(noise(rng),noise(rng),noise(rng));
      V w=bias+(mode==6?.0023f:.0015f)*V(noise(rng),noise(rng),noise(rng));
      if(mode==4) {a=Eigen::AngleAxisf(-angle,V::UnitX())*a;w.x()+=rate;} else w.z()+=rate;
      V m=Eigen::AngleAxisf(-angle,V::UnitZ())*V(20,0,42)+.08f*V(noise(rng),noise(rng),noise(rng));
      if(mode==7 && ms<12000)m=V(20,0,42); // frozen magnetic register
      const bool missing=mode==8 && ms>=3000 && ms<3300;
      const bool missing_mag=mode==4 || (mode==9 && ms<3000);
      status=capture.update(start+uint32_t(ms),missing?nullptr:&a,missing?nullptr:&w,31+.002f*noise(rng),missing_mag?nullptr:&m);
      if(status==GS::READY){ready_ms=ms;break;}
    }
    check(status==GS::READY,"gyro completes after a sustained quiet interval");
    check(cal.buf.n>=220 && cal.buf.n<=245,"retained gyro samples respect spacing across block boundaries");
    if(mode==1 || mode==4)check(ready_ms>=18000,"steady yaw/roll is not saved as bias");
    if(mode==2)check(ready_ms>=6500,"half-moving capture is discarded before still hold");
    if(mode==3)check(ready_ms>=12000 && capture.resets()>0,"movement resets accumulated quiet capture");
    if(mode==7)check(ready_ms>=18000,"frozen magnetometer cannot establish stationary yaw");
    if(mode==8)check(ready_ms>=9300,"sample gap resets continuous quiet interval");
    if(mode==9)check(ready_ms>=9000,"late magnetic reference starts a new complete hold");
    imu_cal::GyroCalibration<float> out;
    check(cal.fit(out) && (out.biasT.b0-bias).norm()<.0007f,"recovered gyro calibration estimates bias only");
    std::printf("gyro mode=%d ready=%dms n=%d restarts=%u bias_error=%g\n",mode,ready_ms,cal.buf.n,capture.resets(),(out.biasT.b0-bias).norm());
  }
  imu_cal::GyroCalibrator<float,400> cal;cal.max_gyro_norm=.12f;
  for(int i=0;i<220;++i)cal.addSample(V(0,0,i<110?.09f:0),V(0,0,-cal.g),25);
  imu_cal::GyroCalibration<float> out;imu_cal::FitFail why;
  check(!cal.fit(out,&why) && why==imu_cal::FitFail::GYRO_NOT_STILL,"direct narrow-temperature fit rejects motion averaged into bias");
  imu_cal::GyroCapture<float,400> cap(cal);cap.begin(0);
  V a(0,0,-cal.g),w(.01f,0,0);
  for(uint32_t t=0;t<12000;t+=5)cap.update(t,&a,&w,25);
  check(cap.update(12001,&a,&w,25)==GS::STALE && cal.buf.n==0,"frozen register values never qualify a still hold");
  cap.begin(0);check(cap.update(12001,nullptr,nullptr,NAN)==GS::STALE,"missing gyro data reaches stale failure");
}

static void testGyroMagneticNoise() {
  // The IMU is polled at 200 Hz, but magnetic registers change more slowly.
  // Repeated magnetic words are one observation, not independent noise draws.
  // A stationary device must complete without waiting for a lucky quiet run.
  int completed=0,first_hold=0;
  for(const int mag_ms : {10,50,100}) for(int seed=0;seed<30;++seed) {
    imu_cal::GyroCalibrator<float,400> cal;cal.max_gyro_norm=.12f;cal.max_accel_dev=.8f;
    imu_cal::GyroCapture<float,400> capture(cal);capture.begin(0);
    std::mt19937 rng(27100+seed);std::normal_distribution<float> noise(0,1);
    const V bias(.002f,-.003f,.001f);
    V m=V::Zero();GS status=GS::SETTLING;
    for(int ms=0;ms<=8500 && status!=GS::READY;ms+=5) {
      V a=V(0,0,-cal.g)+.01f*V(noise(rng),noise(rng),noise(rng));
      V w=bias+.0023f*V(noise(rng),noise(rng),noise(rng));
      if(ms%mag_ms==0)m=V(20,0,42)+.35f*V(noise(rng),noise(rng),noise(rng));
      status=capture.update(uint32_t(ms),&a,&w,31,&m);
    }
    completed+=status==GS::READY;
    first_hold+=status==GS::READY && capture.resets()==0;
    check(status==GS::READY && capture.resets()==0,"stationary multi-rate magnetic noise does not restart gyro hold");
    imu_cal::GyroCalibration<float> out;
    check(cal.fit(out) && (out.biasT.b0-bias).norm()<.0007f,"noisy magnetic qualification retains gyro bias accuracy");
  }
  std::printf("gyro magnetic noise: %d/90 ready, %d/90 uninterrupted\n",completed,first_hold);
}

static void testGyroMagneticInterruptions() {
  // Longer magnetic averaging must still detect a fixed-reference drift and
  // must not treat a lost/frozen yaw reference as an absent sensor.
  for(int mode=0;mode<4;++mode) for(int seed=0;seed<10;++seed) {
    imu_cal::GyroCalibrator<float,400> cal;cal.max_gyro_norm=.12f;cal.max_accel_dev=.8f;
    imu_cal::GyroCapture<float,400> capture(cal);capture.begin(0);
    std::mt19937 rng(49100+seed);std::normal_distribution<float> noise(0,1);
    const V bias(.002f,-.003f,.001f);
    V m=V::Zero();GS status=GS::SETTLING;
    float angle=0;int ready_ms=-1;bool magnetic_restart=false;
    for(int ms=0;ms<25000;ms+=5) {
      const float rate=mode==0 && ms<12000 ? .015f : 0;
      angle+=rate*.005f;
      V a=V(0,0,-cal.g)+.01f*V(noise(rng),noise(rng),noise(rng));
      V w=bias+V(0,0,rate)+.0015f*V(noise(rng),noise(rng),noise(rng));
      const bool interrupted=ms>=3000 && ms<12000;
      if(ms%100==0 && !(mode==2 && interrupted)) {
        m=Eigen::AngleAxisf(-angle,V::UnitZ())*V(20,0,42)+.35f*V(noise(rng),noise(rng),noise(rng));
        if(mode==1 && interrupted)m+=V(6,0,0);
      }
      status=capture.update(uint32_t(ms),&a,&w,31,(mode==3 && interrupted)?nullptr:&m);
      magnetic_restart=magnetic_restart || std::strcmp(capture.resetReason(),"magnetic change")==0 ||
          std::strcmp(capture.resetReason(),"magnetic stream lost")==0;
      if(status==GS::READY) {ready_ms=ms;break;}
    }
    // After the field step it is legitimate to capture in the new stable
    // field, but samples from before the step must not contribute.
    const int earliest_ready=mode==1?9000:18000;
    if(status!=GS::READY || ready_ms<earliest_ready || !magnetic_restart)
      std::printf("gyro magnetic interruption mode=%d seed=%d ready=%d restarts=%u reason=%s\n",
                  mode,seed,ready_ms,capture.resets(),capture.resetReason());
    check(status==GS::READY && ready_ms>=earliest_ready && magnetic_restart,
          "slow yaw, field change and lost/frozen magnetic stream require a new quiet hold");
    imu_cal::GyroCalibration<float> out;
    check(cal.fit(out) && (out.biasT.b0-bias).norm()<.0007f,"interrupted hold never contaminates saved gyro bias");
  }
}

static void testAccelPoseProgress() {
  using Phase=imu_cal::AccelHoldPhase;
  imu_cal::AccelCaptureCfg cfg;
  imu_cal::AccelRegionSpec region;region.expected=V(0,0,-1);
  // A slow placement used to fill the bar at 6.5 s with no accepted data,
  // leaving a full bar until the 30 s timeout or the user found the pose.
  for(const bool recheck : {false,true}) {
    imu_cal::AccelObs obs[40];int n=0;
    imu_cal::AccelHoldEngine<40> cap;
    cap.begin(cfg,region,obs,&n,0,cfg.hold_useful_ms,recheck?cfg.recheck_verify_ms:0,V::Zero(),true);
    uint32_t t=0;
    auto feed=[&](bool right_pose, bool moving=false) {
      imu_cal::AccelRawSample s;s.t_us=t*1000;s.tempC=25;
      s.a=V(.002f*std::sin(.07f*t),0,right_pose?-cfg.g:cfg.g);
      s.w=V(moving?.4f:.0002f*std::cos(.03f*t),0,0);
      cap.feed(s);t+=5;
    };
    while(t<8000)feed(false);
    check(n==0 && cap.phase()==Phase::PLACING && cap.view().progress01==0,
          "eight-second placement never looks like a completed pose");
    check(cap.view().hint==imu_cal::AccelHoldHint::CHECK_POSE,"wrong pose has actionable guidance");
    while(t<15000 && n<19)feed(true);
    const float before=cap.view().progress01;
    check(n==19 && before<1,"nearly complete pose still needs accepted data");
    const uint32_t disturbance_end=t+1500;
    while(t<disturbance_end)feed(true,true);
    check(n==19 && cap.phase()==Phase::PLACING && cap.view().progress01==before,
          "late motion preserves accepted progress without claiming completion");
    bool honest=true;
    while(t<22000 && cap.phase()!=Phase::DONE && cap.phase()!=Phase::FAILED) {
      feed(true);
      honest=honest && (cap.view().progress01<1 || cap.phase()==Phase::DONE);
    }
    const int expected=recheck?28:20;
    check(cap.phase()==Phase::DONE && n==expected && honest && cap.view().progress01==1,
          "100 percent coincides with completion including verification blocks");
    int verify=0;
    for(int i=0;i<n;++i)verify+=obs[i].role==(uint8_t)imu_cal::AccelObsRole::VERIFY;
    check(verify==(recheck?8:0),"five useful seconds and independent verification are retained");
  }
}

struct CalSerial {
  bool connected=true;
  int room=256,writes=0,waits=0;
  std::string output;
  explicit operator bool() const {return connected;}
  int availableForWrite() const {return room;}
  size_t write(const uint8_t* data,size_t n) {
    ++writes;
    if(room<0 || n>(size_t)room) {++waits;return 0;}
    room-=(int)n;output.append((const char*)data,n);return n;
  }
};

static void testAccelCompletionLogging() {
  using atoms3r_ical::tryCalLogLine;
  CalSerial serial;
  check(tryCalLogLine(serial,"[ACC] OK") && serial.output=="[ACC] OK\r\n",
        "available serial output keeps complete diagnostic lines");
  serial.room=0;
  for(int i=0;i<28;++i)tryCalLogLine(serial,"[ACCBLK] retained observation");
  check(serial.writes==1 && serial.waits==0,"pose-end burst never writes into a full USB buffer");
  serial.room=100;serial.connected=false;
  check(!tryCalLogLine(serial,"[ACC] OK") && serial.writes==1,"disconnected USB does not delay the wizard");
  serial.connected=true;serial.room=9;
  check(!tryCalLogLine(serial,"[ACC] OK") && serial.writes==1,"capacity check includes the line terminator");
  serial.room=10;
  check(tryCalLogLine(serial,"[ACC] OK"),"an exactly fitting line is delivered");
  serial.room=-1;
  check(!tryCalLogLine(serial,"[ACC] OK"),"serial error cannot become an unsigned capacity");

  // Run a real first pose with a serial host that never drains its buffer.
  // The next preparation screen must follow the final accepted block, with
  // no fit or diagnostic transport wait in between.
  struct Io : imu_cal::AccelCalIo {
    uint32_t t=0,done_at=0,ok_at=0,next_at=0;int preps=0;
    CalSerial serial;
    Io(){serial.room=0;}
    bool prep(const imu_cal::AccelStepView&) override {if(++preps==1)return true;next_at=t;return false;}
    bool sample(imu_cal::AccelRawSample& s) override {
      s.t_us=t*1000;s.a=V(.002f*std::sin(.07f*t),0,-9.80665f);
      s.w=V(.0002f*std::cos(.03f*t),0,0);s.tempC=25;return true;
    }
    uint32_t nowMs() override {return t;}
    void capture(const imu_cal::AccelStepView&,const imu_cal::AccelHoldView& h) override {
      if(h.phase==imu_cal::AccelHoldPhase::DONE)done_at=t;
    }
    void holdOk(const imu_cal::AccelStepView&) override {ok_at=t;t+=980;}
    bool holdRetry(const imu_cal::AccelStepView&,const char*) override {return false;}
    bool runFit(imu_cal::AccelFitJob&,const char*) override {check(false,"no fit between individual poses");return false;}
    void log(const char* line) override {tryCalLogLine(serial,line);}
    void idle() override {t+=5;}
  } io;
  imu_cal::AccelCalProcedure<> proc;
  proc.begin(imu_cal::AccelCaptureCfg{},imu_cal::AccelFitCfg{},imu_cal::AccelThermalPrior{},V::Zero(),true);
  check(!proc.runMainStage(io) && io.preps==2,"timing probe reaches the next pose");
  check(proc.nObs()==20 && io.done_at<7000 && io.ok_at-io.done_at<=5,
        "quiet pose finishes promptly with every required block");
  check(io.next_at-io.done_at<=985 && io.serial.writes==0 && io.serial.waits==0,
        "full USB buffer cannot add a pause before the next pose");
}

int main() {
  testMagneticQuality();testGuidance();testGyroCapture();testGyroMagneticNoise();testGyroMagneticInterruptions();
  testAccelPoseProgress();testAccelCompletionLogging();
  std::printf("calibration_workflow-test: %d/%d checks passed\n",checks-failures,checks);
  return failures?1:0;
}
