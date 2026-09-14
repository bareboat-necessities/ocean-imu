#!/usr/bin/env python3
"""Native shipping correspondence for the Live preserve-yaw hard reset.

This is an implementation regression only. It compiles and executes the actual
Kalman3D_Wave_OU_III header on two deterministic cases chosen to have simple
ideal-real witnesses:

* a retained-yaw plus 24-7-25 accelerometer tilt, checking the final quaternion,
  the accel-only intermediate covariance axis, and zero attitude cross blocks;
* a nonzero but sub-1e-8 rotation axis, checking shipping's near-parallel branch.

It does not qualify libm for all admitted inputs, does not bound roundoff, and
cannot enable any ALT theorem or storage gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
HEADER = ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'


def eigen_include() -> Path:
    candidates = [os.environ.get('EIGEN_INCLUDE_DIR',''), '/usr/include/eigen3',
                  '/usr/local/include/eigen3',
                  '/opt/pyvenv/lib/python3.13/site-packages/casadi/include/eigen3']
    for raw in candidates:
        p=Path(raw)
        if raw and (p/'Eigen/Core').is_file(): return p
    raise RuntimeError('Eigen3 not found; set EIGEN_INCLUDE_DIR')


SOURCE = r'''
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
const float g_std = 9.80665f;

int main() {
    using Filter = Kalman3D_Wave_OU_III<float,true,true>;
    using V3 = Eigen::Vector3f;
    Filter m(V3::Constant(.0148f), V3::Constant(.00157f), V3::Constant(.3f));

    // BODY->WORLD old attitude: pure yaw with half-angle (4/5,3/5).
    const Eigen::Quaternionf q_old(4.f/5.f, 0.f, 0.f, 3.f/5.f);
    m.initialize_from_attitude(q_old, .035f, 1.5708f);

    // At this exact-real geometry quaternion_from_acc installs
    // qref_tilt=(4/5,-3/5,0,0), hence u_down=(0,24/25,7/25).
    const V3 acc(0.f, -24.f/25.f*g_std, -7.f/25.f*g_std);
    m.initialize_from_acc_preserve_yaw(acc);

    Eigen::Quaternionf got = m.quaternion_boat();
    got.normalize();
    Eigen::Quaternionf expected_q(16.f/25.f, 12.f/25.f, 9.f/25.f, 12.f/25.f);
    expected_q.normalize();
    if (got.coeffs().dot(expected_q.coeffs()) < 0.f) got.coeffs() *= -1.f;
    const float qerr = (got.coeffs()-expected_q.coeffs()).cwiseAbs().maxCoeff();

    const V3 u(0.f,24.f/25.f,7.f/25.f);
    const Eigen::Matrix3f I=Eigen::Matrix3f::Identity();
    const Eigen::Matrix3f expected_P =
        .035f*.035f*(I-u*u.transpose()) + 1.5708f*1.5708f*(u*u.transpose());
    const auto P=m.covariance_full();
    const float perr=(P.block<3,3>(0,0)-expected_P).cwiseAbs().maxCoeff();
    const float cross=P.block<3,18>(0,3).cwiseAbs().maxCoeff();

    // Nonzero transverse component, but below shipping norm_axis < 1e-8.
    const float t=1e-9f;
    const float den=1.f+t*t;
    const V3 target(0.f,2.f*t/den,(1.f-t*t)/den);
    const Eigen::Quaternionf near=Filter::quaternion_from_acc(-target);
    const float near_err=(near.coeffs()-Eigen::Quaternionf::Identity().coeffs()).cwiseAbs().maxCoeff();

    std::cout << std::setprecision(9)
              << qerr << ' ' << perr << ' ' << cross << ' ' << near_err << '\n';
    return (qerr <= 2e-6f && perr <= 2e-6f && cross == 0.f && near_err == 0.f) ? 0 : 2;
}
'''


def run(output: Path) -> dict:
    if not HEADER.is_file(): raise RuntimeError('shipping OU-III header missing')
    with tempfile.TemporaryDirectory(prefix='ou3-alt-tilt-native-') as td:
        td=Path(td); src=td/'tilt.cpp'; exe=td/'tilt'
        src.write_text(SOURCE)
        command=[os.environ.get('CXX','g++'),'-std=c++20','-O1','-fno-fast-math',
                 '-ffp-contract=off','-DEIGEN_NON_ARDUINO','-I'+str(eigen_include()),
                 '-I'+str(ROOT/'src'),str(src),'-o',str(exe)]
        subprocess.run(command,check=True,cwd=ROOT)
        proc=subprocess.run([str(exe)],check=True,capture_output=True,text=True,cwd=ROOT)
    vals=proc.stdout.strip().split()
    if len(vals)!=4: raise RuntimeError('unexpected native tilt-reset output')
    qerr,perr,cross,near_err=map(float,vals)
    result={
      'qualification':'OU3_ALT_NATIVE_TILT_RESET_IMPLEMENTATION_REGRESSION_ONLY',
      'shipping_header_sha256':hashlib.sha256(HEADER.read_bytes()).hexdigest(),
      'preserve_yaw_quaternion_max_abs_error':qerr,
      'intermediate_axis_covariance_max_abs_error':perr,
      'attitude_cross_covariance_max_abs':cross,
      'near_parallel_identity_max_abs_error':near_err,
      'same_shipping_header_executed':True,
      'all_input_libm_binary32_correspondence_qualified':False,
      'source_uniform_word_qualified':False,
      'storage_feasibility_attempted':False,
      'deployment_roundoff_enclosed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
    if qerr>2e-6 or perr>2e-6 or cross!=0 or near_err!=0:
        raise AssertionError(result)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    print(json.dumps(run(args.output),indent=2,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
