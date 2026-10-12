"""All estimator sketches must carry the shared motion/freshness policy.

The acquisition test compiles definitions extracted from the shipping header;
its fake driver controls arrival masks and AUX bursts, not estimator behavior.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from test_device_compass_startup import function

ROOT = Path(__file__).resolve().parents[2]


class SensorDeliveryTest(unittest.TestCase):
    def test_diagnostic_rotation_and_freshness(self):
        """Exercise the diagnostic's actual integration and source-epoch code."""
        eigen = Path(os.environ.get('EIGEN_INCLUDE_DIR', '/usr/include/eigen3'))
        if not (eigen / 'Eigen/Dense').is_file():
            eigen = ROOT / 'third_party/eigen'
        self.assertTrue((eigen / 'Eigen/Dense').is_file(), 'Install Eigen or set EIGEN_INCLUDE_DIR')
        src = (ROOT / 'sensors/imu_basic/atomS3R_mag_diagnostics/'
               'atomS3R_mag_diagnostics.ino').read_text()
        sample = function(src, 'onSample')
        integration = sample[:sample.index('  // Use the same acquisition epoch')] + '}\n'
        fresh_start = sample.index('  const bool fresh_m =')
        fresh_end = sample.index('  if (fresh_m)', fresh_start)
        baseline_start = sample.index('    last_m_sequence =')
        baseline_end = sample.index('\n  }', baseline_start)
        cpp = r'''
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cstdlib>
#include <limits>
#include "util/ImuLoopHelpers.h"
using V3=Eigen::Vector3f;
struct ImuSample {
  V3 a=V3::Zero(), w=V3::Zero(), m=V3::Zero();
  float tempC=25; uint32_t sample_us=1000;
  struct {uint32_t sequence=0;} mag_timing;
};
struct Calibration {
  V3 applyAccel(const V3& v,float){return v;}
  V3 applyGyro(const V3& v,float){return v;}
} cals;
float last_temp_c=0;
V3 a_cal=V3::Zero(), w_cal=V3::Zero();
Eigen::Quaternionf q_gyro=Eigen::Quaternionf::Identity();
ocean_imu::ins::SampleDtTracker sample_dt{0.005f};
int ring_n=5, ring_head=4;
bool have_m=false;
uint32_t last_m_sequence=0, last_m_change_us=0;
''' + integration + '\nbool observe(const ImuSample& s) {\n' + sample[fresh_start:fresh_end] + r'''
  if (fresh_m) {
''' + sample[baseline_start:baseline_end] + r'''
  }
  return fresh_m;
}
void require(bool ok){if(!ok)std::abort();}
int main(){
  ImuSample s;
  s.w=V3(0,0,10);
  onSample(s);
  require(ring_n==5);
  s.sample_us+=75000;
  onSample(s);
  const Eigen::Quaternionf expected(Eigen::AngleAxisf(0.8f,V3::UnitZ()));
  require(q_gyro.angularDistance(expected)<1e-6f);
  require(ring_n==0 && ring_head==0);
  const auto saved=q_gyro;
  s.w.x()=std::numeric_limits<float>::quiet_NaN();
  s.sample_us+=5000; ring_n=4;
  onSample(s);
  require(q_gyro.coeffs().allFinite() && q_gyro.angularDistance(saved)<1e-6f);
  require(ring_n==0);
  // An invalid first sample cannot poison future source comparisons.
  s.m.setConstant(std::numeric_limits<float>::quiet_NaN());
  s.mag_timing.sequence=1;
  require(!observe(s));
  s.m=V3(20,0,40);
  require(observe(s) && !observe(s));
  // Equal quantized field values may still belong to a fresh AUX frame.
  s.mag_timing.sequence=2;
  require(observe(s));
  s.mag_timing.sequence=3; s.m.x()=std::numeric_limits<float>::quiet_NaN();
  require(!observe(s));
  s.m.x()=20;
  require(observe(s));
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'diagnostic.cpp'
            path.write_text(cpp)
            binary = Path(tmp) / 'diagnostic'
            subprocess.run(['g++', '-std=c++20', '-O1', '-I'+str(ROOT/'src'), '-I'+str(eigen),
                            str(path), '-o', str(binary)], check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True)

    def test_shipping_acquisition_mask_and_single_burst(self):
        """Compile shipping acquisition helpers and check delivery masks and AUX freshness."""
        eigen = Path(os.environ.get('EIGEN_INCLUDE_DIR', '/usr/include/eigen3'))
        if not (eigen / 'Eigen/Dense').is_file():
            eigen = ROOT / 'third_party/eigen'
        self.assertTrue((eigen / 'Eigen/Dense').is_file(), 'Install Eigen or set EIGEN_INCLUDE_DIR')
        src = (ROOT / 'src/AtomS3R/AtomS3R_ImuCal.h').read_text()
        def struct(name):
            """Extract the named source struct through its closing declaration."""
            start = src.index('struct ' + name + ' {')
            return src[start:src.index('\n};', start) + 3]
        cpp = r'''
#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "util/ImuFreshness.h"
#include "AtomS3R/AtomS3R_MagAcquisition.h"
using Vector3f=Eigen::Vector3f;
uint32_t clock_us=1000;
uint32_t micros(){return clock_us;}
uint32_t millis(){return clock_us/1000;}
namespace m5 {
struct imu_3d_t {float x=0,y=0,z=0;};
struct IMU_Base {
  uint8_t bytes[8]={9,2,3,4,5,6,7,8};
  bool good=true; unsigned reads=0;
  bool readRegister(int,uint8_t* out,size_t n) {++reads;if(good)std::memcpy(out,bytes,n);return good;}
};
}
struct Data {m5::imu_3d_t accel,gyro,mag;};
struct Driver {Data d;Data getImuData(){return d;} bool getTemp(float* t){*t=35;return true;}};
struct {Driver Imu;} M5;
namespace atoms3r_ical {
struct Timing {uint32_t delay_us=10000,configuration=7;float fit_width_s=.001f,residual_time_s=.002f;};
Timing timing;
const Timing& activeMagTiming(){return timing;}
bool magTimingValid(const Timing&){return true;}
uint32_t activeMagTimingConfiguration(){return 7;}
struct Bmm150Trim {};
struct Bmm150Raw {uint8_t x; static Bmm150Raw fromData(const uint8_t* d){return {d[0]};}};
bool bmm150Compensate(Bmm150Raw r,const Bmm150Trim&,float* v){v[0]=r.x;v[1]=20;v[2]=40;return true;}
void bmm150SensorToAtomS3RBody(const float* a,float* b){std::memcpy(b,a,3*sizeof(float));}
Vector3f map_acc_to_body_ned_(m5::imu_3d_t v){return {v.x,v.y,v.z};}
Vector3f map_gyr_to_body_ned_(m5::imu_3d_t v){return {v.x,v.y,v.z};}
Vector3f map_mag_to_body_uT_(m5::imu_3d_t v){return {v.x,v.y,v.z};}
#define ATOMS3R_IMU_MASK_ACCEL 1u
#define ATOMS3R_IMU_MASK_GYRO 2u
constexpr uint32_t kImuMaskAccelGyro=3u;
''' + struct('ImuSample') + '\n' + struct('AtomS3RMagSource') + r'''
AtomS3RMagSource source;
AtomS3RMagSource& atoms3rMagSource(){return source;}
''' + function(src, 'readMagBody_') + '\n' + function(src, 'readImuGyroMapped') + '\n' + function(src, 'readImuMapped') + r'''
}
void require(bool ok){if(!ok)std::abort();}
int main(){
  using namespace atoms3r_ical;
  m5::IMU_Base aux;source.compensated=true;source.bmi270=&aux;
  ImuSample s;
  require(!readImuMapped(M5.Imu,2,clock_us,s));
  require(!readImuMapped(M5.Imu,1,clock_us,s));
  require(!readImuMapped(M5.Imu,2,clock_us,s)&&s.accel_fresh);
  require(s.m.x()==8 && source.acquisition.sequence==1);
  require(s.mag_timing.delay_us==10000 && std::abs(s.mag_timing.timing_allowance_s-.003f)<1e-6f);
  clock_us+=5000;
  require(!readImuMapped(M5.Imu,2,clock_us,s)&&!s.accel_fresh);
  require(source.acquisition.sequence==1);
  // M5's magnetic cache did not change; the actual AUX burst did.
  aux.bytes[0]=17;clock_us+=5000;
  require(readImuMapped(M5.Imu,3,clock_us,s)&&s.accel_fresh);
  require(s.m.x()==16 && source.acquisition.sequence==2);
  require(source.acquisition.frame_us==clock_us);
  // AUX status bits are not a new magnetic payload.
  aux.bytes[0]=19;clock_us+=5000;
  require(!readImuMapped(M5.Imu,2,clock_us,s));
  require(source.acquisition.sequence==2 && s.m.x()==16);
  timing.configuration=8;clock_us+=5000;
  require(!readImuMapped(M5.Imu,2,clock_us,s) && s.mag_timing.delay_us==0);
  aux.good=false;clock_us+=5000;
  require(!readImuMapped(M5.Imu,2,clock_us,s));
  require(source.acquisition.sequence==2);
  clock_us+=300000;
  require(readImuMapped(M5.Imu,3,clock_us,s));
  require(!s.m.allFinite());
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'acquisition.cpp'
            path.write_text(cpp)
            binary = Path(tmp) / 'acquisition'
            subprocess.run(['g++', '-std=c++20', '-O1', '-I'+str(ROOT/'src'), '-I'+str(eigen),
                            str(path), '-o', str(binary)], check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    unittest.main()
