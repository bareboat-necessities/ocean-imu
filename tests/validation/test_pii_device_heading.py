"""Exercise the PII sketch update and the real shared compass renderer."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from test_device_compass_startup import function

ROOT = Path(__file__).resolve().parents[2]
SKETCH = ROOT / 'sensors/full_marine_ins/atomS3R_ins_pii_observer/atomS3R_ins_pii_observer.ino'


class PiiDeviceHeadingTest(unittest.TestCase):
    def test_invalid_compass_does_not_report_north(self):
        stub = r'''
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
constexpr int TFT_BLACK=0,TFT_DARKGREY=1,TFT_LIGHTGREY=2,TFT_WHITE=3,
              TFT_RED=4,TFT_GREEN=5,TFT_ORANGE=6;
inline std::vector<std::string> printed;
inline int needles=0;
struct Screen {int width(){return 128;}int height(){return 128;}};
inline struct {Screen Display;} M5;
struct M5Canvas {
  std::vector<uint8_t> buffer;
  void setColorDepth(int){}
  void* createSprite(int w,int h){buffer.resize(unsigned(w*h*2));return buffer.data();}
  void deleteSprite(){buffer.clear();}
  void* getBuffer(){return buffer.data();}
  void pushSprite(Screen*,int,int){}
  void fillSprite(int){}
  void drawCircle(int,int,int,int){}
  void drawLine(int,int,int,int,int){}
  void fillTriangle(int,int,int,int,int,int,int){++needles;}
  void fillCircle(int,int,int,int){}
  void fillRect(int,int,int,int,int){}
  void setTextColor(int,int){}
  void setTextSize(int){}
  int textWidth(const char* s){return int(std::strlen(s))*6;}
  int fontHeight(){return 8;}
  void setCursor(int,int){}
  void print(const char* s){printed.emplace_back(s);}
};
'''
        cpp = r'''
#include "AtomS3R/AtomS3R_CompassUI.h"
#include <algorithm>
#include <cassert>
#include <limits>
bool text(const char* s){return std::find(printed.begin(),printed.end(),s)!=printed.end();}
int main(){
  CompassUI ui;ui.begin();assert(ui.ok());
  printed.clear();needles=0;ui.draw(137,true,51,false);
  assert(needles==1&&text("137")&&!text("WAIT"));
  for(float heading:{0.f,137.f,std::numeric_limits<float>::quiet_NaN(),
                     std::numeric_limits<float>::infinity()}){
    printed.clear();needles=0;ui.draw(heading,false,51,false);
    assert(needles==0&&text("---")&&text("WAIT")&&!text("0")&&!text("N"));
  }
  printed.clear();needles=0;ui.draw(std::numeric_limits<float>::quiet_NaN(),true);
  assert(needles==0&&text("WAIT"));
  printed.clear();needles=0;ui.draw(0,true);
  assert(needles==1&&text("0")&&text("N"));
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / 'M5Unified.h').write_text(stub)
            source = folder / 'compass.cpp'
            binary = folder / 'compass'
            source.write_text(cpp)
            subprocess.run(['g++', '-std=c++20', '-Wall', '-Wextra', '-Werror',
                            '-I'+str(folder), '-I'+str(ROOT/'src'), str(source),
                            '-o', str(binary)], check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True, timeout=5)

    def test_pii_sketch_heading_through_rotation_and_field_outage(self):
        eigen = Path(os.environ.get('EIGEN_INCLUDE_DIR', '/usr/include/eigen3'))
        if not (eigen / 'Eigen/Dense').is_file():
            eigen = ROOT / 'third_party/eigen'
        self.assertTrue((eigen / 'Eigen/Dense').is_file(), 'Install Eigen or set EIGEN_INCLUDE_DIR')
        source = SKETCH.read_text()
        cpp = r'''
#define EIGEN_NON_ARDUINO
extern const float g_std=9.80665f;
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "pii_observer/AdaptiveVerticalPIIMahony.h"
#include "util/MagneticInput.h"
#include "util/MagneticHeading.h"
#include "util/ImuLoopHelpers.h"
#include "detrend/AdaptiveWaveDetrender.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <cassert>
using Vector3f=Eigen::Vector3f;using Matrix3f=Eigen::Matrix3f;using Quaternionf=Eigen::Quaternionf;
namespace ins=ocean_imu::ins;
constexpr float APP_G_LOCAL=9.8025605f,APP_FREQ_GUESS=.3f,MAG_PRESENT_MIN_UT=5,
 MAG_PRESENT_MAX_UT=200,MAG_FIELD_MIN_UT=20,MAG_FIELD_MAX_UT=80,MAHONY_SEED_G_TOL_FRAC=.2f;
constexpr uint32_t HEADING_MAG_TIMEOUT_MS=2000;
#define SEA_STATE_USE_STRICT_MAG_FIELD_GATE 0
#define SEA_STATE_OUTPUT_TRUE_HEADING 0
#define SEA_STATE_MAG_HEADING_USER_OFFSET_DEG 0.f
uint32_t clock_us=1000000;uint32_t micros(){return clock_us;}uint32_t millis(){return clock_us/1000;}
struct ImuSample {Vector3f a=Vector3f::Zero(),w=Vector3f::Zero(),m=Vector3f::Zero();float tempC=25;uint32_t sample_us=0;bool accel_fresh=true;
 atoms3r_ical::MagAcquisition mag_timing;};
struct App {
 using Fusion=marine_obs::AdaptiveVerticalPIIMahony<float,true,TrackerType::PLL>;
 Fusion fusion_;
 bool have_blob_=true,mahony_seeded_=false,mag_present_=false,mag_field_sane_=false,
 mag_fresh_=false,mag_used_=false,heading_valid_=false;
 struct Runtime {Vector3f applyAccel(const Vector3f& a,float){return a;}
 Vector3f applyGyro(const Vector3f& w,float){return w;}Vector3f applyMag(const Vector3f& m){return m;}} runtime_;
 ocean_imu::magnetic::Input mag_input_;
 ins::MagFreshGate mag_gate_;
 ins::SampleDtTracker sample_dt_{.005f};
 ins::SampleRateMeter rates_;
 ins::StillGyroBiasEma gyro_bias_{5,.12f,.15f};ins::RateOfTurnFilter rot_;
 AdaptiveWaveDetrender z_detrender_;
 uint32_t last_mag_correction_ms_=0;
 Vector3f a_cal_,w_cal_,m_cal_;
 float dt_=0,imu_temp_c_=0,a_raw_norm_=0,mag_norm_uT_=0,roll_deg_=0,pitch_deg_=0,
 yaw_mahony_deg_=0,heading_mag_deg_=0,heading_deg_=0,heave_m_=0,heave_speed_mps_=0,
 heave_raw_m_=0,wave_hz_=APP_FREQ_GUESS,wave_envelope_m_=0,heave_baseline_m_=0,
 heave_wave_raw_m_=0,heave_wave_clean_m_=0;
 App(){Fusion::Config cfg;cfg.gravity_mps2=APP_G_LOCAL;fusion_.configure(cfg);fusion_.reset();
 mag_input_.reset(0,Eigen::Matrix3f::Identity());}
''' + function(source, 'ned_to_mahony_body_') + '\n' + function(source, 'seedMahonyFromAccelMag_') + '\n' + function(source, 'updateFilter_') + r'''
};
int main(){
 App app;ImuSample s;
 Quaternionf truth(Eigen::AngleAxisf(137.f/57.295779513f,Vector3f::UnitZ()));
 bool timed_out=false,resumed=false;
 for(int i=1;i<=5000;++i){
   clock_us+=5000;s.sample_us=clock_us;s.accel_fresh=i%2==0;
   const float rate=i>1000&&i<=1400?3.14159265f:0;
   const Vector3f w(0,0,rate);
   truth=(truth*Quaternionf(Eigen::AngleAxisf(rate*.005f,Vector3f::UnitZ()))).normalized();
   s.a=truth.conjugate()*Vector3f(0,0,-APP_G_LOCAL);s.w=w;
   const bool field_available=i<=4000||i>4800;
   if(i%8==0&&field_available){s.m=truth.conjugate()*Vector3f(20.3f,0,46.8f);
     s.mag_timing.acquired(clock_us/1000,clock_us,clock_us,clock_us);}
   if(s.accel_fresh) app.updateFilter_(s);
   if(i==999)assert(app.heading_valid_&&std::abs(app.heading_deg_-137)<1);
   if(i>1000&&i<=2000){
     assert(app.mag_used_&&app.heading_valid_);
     const float actual=std::atan2((truth*Vector3f::UnitX()).y(),(truth*Vector3f::UnitX()).x())*57.295779513f;
     assert(std::abs(std::remainder(app.heading_deg_-actual,360.f))<2);
   }
   if(i>4600&&i<=4800){assert(!app.mag_used_&&!app.heading_valid_);timed_out=true;}
   if(i>4900&&app.heading_valid_)resumed=true;
 }
 assert(timed_out&&resumed&&app.mahony_seeded_);
 assert(std::abs(std::remainder(app.heading_deg_-137,360.f))<2);
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'pii.cpp'
            binary = Path(tmp) / 'pii'
            path.write_text(cpp)
            subprocess.run(['g++', '-std=c++20', '-O1', '-I'+str(ROOT/'src'),
                            '-isystem', str(eigen), str(path),
                            '-o', str(binary)],
                           check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True, timeout=10)
