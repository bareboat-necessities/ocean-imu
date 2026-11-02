#!/usr/bin/env python3
"""Exercise the actual compass backend's reference/gating order on a host."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "sensors/compass_ahrs/atomS3R_compass_qmekf/atomS3R_compass_qmekf.ino").read_text()
start = source.index("class QmekfBackend")
brace = source.index("{", start)
depth, end = 1, brace + 1
while depth:
    depth += (source[end] == "{") - (source[end] == "}")
    end += 1
backend = source[start:end] + ";"
eigen = next(p for p in (root / "third_party/eigen", Path(os.environ.get("EIGEN_INCLUDE_DIR", "/usr/include/eigen3")))
             if (p / "Eigen/Dense").exists())
cpp = r'''
#define EIGEN_NON_ARDUINO
#include <memory>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include "AtomS3R/AtomS3R_ImuCalBlob.h"
#include "ahrs/KalmanQMEKF.h"
namespace atoms3r_compass {
using namespace atoms3r_ical;
struct CalibratedSample {
  Vector3f a_cal=Vector3f(0,0,-ImuCalCfg::g_cal_local),w_cal=Vector3f::Zero(),m_cal=Vector3f(20,0,43);
  float dt=.005f;bool mag_ok=true,mag_fresh=true;
};
struct AttitudeSolution {float x=0,y=0,z=0,w=1;};
AttitudeSolution makeAttitudeFromQuat(float x,float y,float z,float w){return {x,y,z,w};}
struct IAttitudeBackend {
  virtual ~IAttitudeBackend()=default;
  virtual void reset()=0;
  virtual void configureNoise(const RuntimeCals&) {}
  virtual void step(const CalibratedSample&,AttitudeSolution&)=0;
  virtual bool isValid() const=0;
};
}
using namespace atoms3r_compass;
''' + backend + r'''
void require(bool ok,const char* message){if(!ok){std::cerr<<message<<'\n';std::exit(1);}}
Eigen::Quaternionf quat(const AttitudeSolution& a){return {a.w,a.x,a.y,a.z};}
int main(){
  QmekfBackend backend;backend.reset();CalibratedSample sample;AttitudeSolution out;
  for(int i=0;i<2000;++i)backend.step(sample,out);
  const auto before=quat(out);
  const unsigned rejected=backend.magneticHealth().rejected;
  sample.m_cal=Vector3f(20,0,73); // a dip/norm disturbance that a per-reading reference could absorb
  for(int i=0;i<100;++i)backend.step(sample,out);
  require(backend.magneticHealth().rejected==rejected+100,"reference adaptation concealed a field disturbance");
  require(quat(out).angularDistance(before)<1e-5f,"disturbance changed the compass attitude");
  sample.m_cal=Vector3f(0,-20,43);
  for(int i=0;i<100;++i)backend.step(sample,out);
  require(quat(out).angularDistance(before)<1e-5f,"90-degree disturbance flipped the consumer heading");
  sample.m_cal=Vector3f(20,0,43);
  for(int i=0;i<80;++i)backend.step(sample,out);
  require(backend.magneticHealth().accepted,"consumer did not recover from a clean field");
  backend.reset();sample.mag_ok=false;
  for(int i=0;i<100;++i)backend.step(sample,out);
  sample.mag_ok=true;sample.m_cal=Vector3f(0,-20,43);backend.step(sample,out);
  require(backend.isValid() && std::isfinite(out.w),"late first magnetic acquisition is invalid");
  std::cout<<"qMEKF consumer: reference protected, dip/norm and 90-degree disturbances rejected, late acquisition valid\n";
}
'''
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "qmekf-reference"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-DEIGEN_UNROLLING_LIMIT=0", "-I" + str(eigen),
        "-I" + str(root / "src"), "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True)
