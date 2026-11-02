#!/usr/bin/env python3
"""Exercise the production alignment worker and its stage-specific diagnostics."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / 'src/AtomS3R/AtomS3R_ImuCalWizard.h').read_text()
def block(marker):
    start = source.index(marker)
    end = source.index('{', start)
    depth = 1
    while depth:
        end += 1
        depth += (source[end] == '{') - (source[end] == '}')
    return source[start:end + 1]
worker = block('static void fitTaskMagAlignment_(')
cpp = r'''
#define EIGEN_NON_ARDUINO
#include "imu_calibrate/GyroCalFit.h"
#include <memory>
#include "AtomS3R/AtomS3R_MagDiagnostics.h"
#include <random>
#include <cstdlib>
#include <iostream>
using V=Eigen::Vector3f;
uint32_t micros(){return 0;}void vTaskDelete(void*){}void vTaskDelay(unsigned){}
struct {template<class... T>void printf(const char*,T...){} } Serial;
struct AccelProc {struct Fitter {
  static void toFloat(const imu_cal::AccelCalibration<float>& in,float,imu_cal::AccelCalibration<float>& out){out=in;}
};};
void require(bool b,const char* s){if(!b){std::cerr<<s<<'\n';std::exit(1);}}
struct ImuCalWizard {
  struct FitCtx {ImuCalWizard* wiz;bool ok=false,done=false;const char* detail="";imu_cal::FitFail reason;};
  struct {imu_cal::AccelCalibration<float> fit;const auto& result(){return fit;}} accel_;
  struct {float g=9.80256f;} accel_fcfg_;
  atoms3r_magdiag::Pose raw_alignment_[24];float alignment_temp_[24]{};
  int alignment_pose_count_=0,alignment_train_count_=0,alignment_validation_count_=0;
  atoms3r_magdiag::Alignment alignment_;float alignment_validation_=INFINITY;
  imu_cal::MagCalibration<float> mag_out_;
  imu_cal::GyroCalibration<float> gyr_out_;
  struct RotationMock {
    bool full_ok=true,scale_ok=true;int full_calls=0,scale_calls=0;
    imu_cal::GyroRotationFail reason=imu_cal::GyroRotationFail::GEOMETRY;
    bool fit(const imu_cal::AccelCalibration<float>&,imu_cal::GyroCalibration<float>& g,void(*)(void*)) {
      ++full_calls;g.matrix_ok=full_ok;g.scale_only=false;return full_ok;
    }
    bool fitScaleOnly(const imu_cal::AccelCalibration<float>&,imu_cal::GyroCalibration<float>& g,void(*)(void*)) {
      ++scale_calls;g.matrix_ok=scale_ok;g.scale_only=scale_ok;
      if(scale_ok){g.S.setIdentity();g.S(0,0)=1.04f;}return scale_ok;
    }
  };
  std::unique_ptr<RotationMock> gyro_rotations_{new RotationMock};
  static uint32_t hwmBytes_(){return 0;}
GYRO_WORKER
WORKER
};
int main(){
  for(int scenario:{0,1,2}) {
    ImuCalWizard gyro;gyro.gyr_out_.ok=gyro.gyr_out_.biasT.ok=true;auto& r=*gyro.gyro_rotations_;r.full_ok=scenario==0;r.scale_ok=scenario==1;
    ImuCalWizard::FitCtx result{&gyro};gyro.fitTaskGyroMatrix_(&result);
    require(result.done && result.ok==(scenario<2) && r.full_calls==1 && r.scale_calls==(scenario==0?0:1),
      "production gyro worker must try three scales automatically only after full fitting fails");
    require(gyro.gyr_out_.scale_only==(scenario==1),"worker must publish the correct qualified model");
  }
  ImuCalWizard w;auto& acc=w.accel_.fit;acc.ok=acc.biasT.ok=true;
  acc.biasT.b0=V(.1f,-.2f,.15f);acc.S.diagonal()<<1.02f,.98f,1.01f;
  w.mag_out_.ok=true;w.mag_out_.field_uT=V(20,0,43).norm();
  w.gyr_out_.S(0,0)=1.03f;w.gyr_out_.rotation_train_segments=14;w.gyr_out_.rotation_validation_segments=6;
  const auto rotation=Eigen::AngleAxisf(1.5707963f,V(.3f,1,.2f).normalized()).toRotationMatrix();
  auto populate=[&](int count){
    std::mt19937 rng(19);std::normal_distribution<float> normal;
    w.alignment_pose_count_=count;
    for(int i=0;i<count;++i){
      Eigen::Quaternionf q(normal(rng),normal(rng),normal(rng),normal(rng));q.normalize();
      w.raw_alignment_[i].a=acc.S.inverse()*(q.conjugate()*V(0,0,-9.80256f))+acc.biasT.b0;
      w.raw_alignment_[i].m_raw=rotation.transpose()*(q.conjugate()*V(20,0,43));
      w.alignment_temp_[i]=25;
    }
  };
  for(int count:{21,13}) {
    populate(count);ImuCalWizard::FitCtx ctx{&w};w.fitTaskMagAlignment_(&ctx);
    require(ctx.done && ctx.ok && w.alignment_.R.isApprox(rotation,.001f),"worker must qualify actual calibrated endpoints");
    require(w.alignment_train_count_==(count==13?9:15) && w.alignment_validation_count_==(count==13?4:6),
      "magnetic fit must persist its actual counts independently of the gyro set");
    require(w.gyr_out_.S(0,0)==1.03f && w.gyr_out_.rotation_validation_segments==6,
      "magnetic worker changed successful gyro data");
  }
  populate(13);w.raw_alignment_[0].m_raw*=1.3f;
  ImuCalWizard::FitCtx field{&w};w.fitTaskMagAlignment_(&field);
  require(!field.ok && std::string(field.detail)=="Mag field changed","field disturbance must be named separately");
  populate(13);for(int i=1;i<13;++i)w.raw_alignment_[i]=w.raw_alignment_[0];
  ImuCalWizard::FitCtx geometry{&w};w.fitTaskMagAlignment_(&geometry);
  require(!geometry.ok && std::string(geometry.detail)=="Mag poses too alike","missing geometry must be named separately");
  populate(13);w.raw_alignment_[3].m_raw=-w.raw_alignment_[3].m_raw;
  ImuCalWizard::FitCtx validation{&w};w.fitTaskMagAlignment_(&validation);
  require(!validation.ok && std::string(validation.detail)=="Mag check failed" && w.alignment_.ok,
    "independent validation must reject a direction error after training succeeds");
  std::cout<<"mag worker: full registration, calibrated endpoints, independent counts and distinct failure causes PASS\n";
}
'''.replace('GYRO_WORKER', block('static void fitTaskGyroMatrix_(')).replace('WORKER', worker)
eigen = root / 'third_party/eigen'
if not (eigen / 'Eigen/Dense').exists():
    eigen = Path('/usr/include/eigen3')
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / 'mag-alignment-task'
    subprocess.run(shlex.split(os.environ.get('CXX', 'g++')) + [
        '-std=c++20', '-O1', '-fno-finite-math-only', '-I', str(root / 'src'),
        '-isystem', str(eigen), '-x', 'c++', '-', '-o', str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=30)
