"""First-sample operand ancestry and native profile checks, not source admission."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import os, shutil, subprocess, tempfile, unittest

from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as X
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V

ROOT=Path(__file__).resolve().parents[2]
CFG=V.Config(X.rn(F(1,5)),X.rn(F(1,50)),X.rn(F(196133,20000)),20)
DT=X.rn(F(1,200))


class Tests(unittest.TestCase):
    def test_first_valid_seed_and_update_share_accelerometer_and_no_seed_port(self):
        acc=tuple(map(X.rn,(F(2,5),F(1,5),-CFG.gravity)))
        out=X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=acc)
        out.startup.validate()
        self.assertEqual(out.startup.acc,acc)
        self.assertEqual(out.startup.branch,'ordinary-FromTwoVectors')
        self.assertTrue(out.vertical.state.initialized)
        self.assertEqual(out.vertical.state.elapsed,DT)
        self.assertEqual(sum(x.kind=='sqrt' for x in out.startup.operations),4)
        with self.assertRaises(TypeError):
            X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=acc,seed=(1,0,0,0))
        with self.assertRaisesRegex(ValueError,'seed detached'):
            replace(out.startup,quaternion=(1,0,0,0)).validate()

    def test_small_first_sample_keeps_reset_and_next_valid_sample_seeds(self):
        s=V.State()
        out=X.step(s,CFG,dt=DT,gyro=(0,0,0),acc=(0,0,0))
        self.assertIs(out.vertical.state,s)
        self.assertEqual(out.startup.branch,'acc-norm-too-small')
        out.startup.validate()
        second=X.step(out.vertical.state,CFG,dt=DT,gyro=(0,0,0),acc=(0,0,-CFG.gravity))
        self.assertTrue(second.vertical.state.initialized)
        self.assertEqual(second.vertical.state.elapsed,DT)

    def test_initialized_successor_does_not_reseed(self):
        a=(0,0,-CFG.gravity)
        first=X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=a)
        second=X.step(first.vertical.state,CFG,dt=DT,gyro=(0,0,0),acc=a)
        self.assertIsNone(second.startup)
        self.assertEqual(second,M.step_initialized(first.vertical.state,CFG,dt=DT,gyro=(0,0,0),acc=a))

    def test_near_antiparallel_branch_computes_and_binds_solver_axis(self):
        acc=(X.rn(0),X.rn(0),CFG.gravity)
        computed=X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=acc)
        self.assertIsNotNone(computed.startup.solver)
        self.assertEqual(computed.startup.solver.sweep_count,1)
        computed.startup.validate()
        svd=X.svd_witness((F(0),F(0),F(-1)),(F(0),F(0),F(1)),
                          (X.rn(1),X.rn(0),X.rn(0)))
        out=X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=acc,svd=svd)
        self.assertEqual(out.vertical,computed.vertical)
        self.assertEqual(out.startup.branch,'near-antiparallel-JacobiSVD')
        self.assertEqual(out.startup.quaternion,(F(0),F(1),F(0),F(0)))
        out.startup.validate()
        bad=replace(svd,v0_dot=F(1))
        with self.assertRaisesRegex(ValueError,'residuals detached'):
            X.step(V.State(),CFG,dt=DT,gyro=(0,0,0),acc=acc,svd=bad)

    def test_uninitialized_nonreset_memory_and_unstored_inputs_rejected(self):
        with self.assertRaisesRegex(ValueError,'reset state'):
            X.step(V.State(up=1),CFG,dt=DT,gyro=(0,0,0),acc=(0,0,-CFG.gravity))
        with self.assertRaisesRegex(ValueError,'binary32'):
            X.step(V.State(),CFG,dt=F(1,200),gyro=(0,0,0),acc=(0,0,-CFG.gravity))

    def test_sqrt_cells_and_residuals_cover_zero_subnormal_and_exponent_edges(self):
        for e in (0,1,2,60,126,127,128,200,253,254):
            for m in (0,1,0x7fffff):
                x=M.value((e<<23)|m); y=X.sqrt32(x)
                X.Operation('sqrt',(x,),y,y*y-x).validate()
        with self.assertRaisesRegex(ValueError,'residual'):
            X.Operation('sqrt',(F(2),),X.sqrt32(2),F(0)).validate()
        with self.assertRaises(ValueError): X.sqrt32(-1)

    def test_readiness_keeps_all_source_and_target_obligations_open(self):
        r=X.readiness()
        self.assertTrue(r['both_shipping_accelerometer_normalizations_retained'])
        self.assertTrue(r['near_antiparallel_JacobiSVD_branch_topology_materialized_with_solver_witness'])
        for key in ('near_antiparallel_JacobiSVD_solver_correspondence_qualified',
                    'target_sqrt_Eigen_and_compiler_correspondence_closed',
                    'every_admitted_startup_history_covered','storage_search_allowed',
                    'ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[key])


HARNESS=r'''
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#define private public
#include "tuner/VerticalAccelComplementary.h"
#undef private
uint32_t bits(float f) { uint32_t b; std::memcpy(&b,&f,4); return b; }
void emit(float f) { std::cout << ' ' << bits(f); }
int main() {
    const float g=9.80665f;
    const Eigen::Vector3f inputs[]={
      {0,0,-g},{.4f,.2f,-g},{.4f,.2f,g},{.5f,0,0},
      {0,0,0},{0,0,-1e-5f},{0,0,-1e-3f},
      {0,0,-std::nextafter(1e-3f,1.0f)}};
    for (const auto& acc : inputs) {
        VerticalAccelComplementary v(.2f,.02f,20.0f);
        const Eigen::Vector3f gyro(.003f,-.004f,.002f);
        for (int k=0;k<2;++k) {
            std::cout << k;
            for (int i=0;i<3;++i) emit(gyro(i));
            for (int i=0;i<3;++i) emit(acc(i));
            v.update(.005f,gyro,acc,g);
            emit(v.ahrs_.q0);emit(v.ahrs_.q1);emit(v.ahrs_.q2);emit(v.ahrs_.q3);
            emit(v.ahrs_.integralFBx);emit(v.ahrs_.integralFBy);emit(v.ahrs_.integralFBz);
            emit(v.elapsed_sec_);emit(v.up_ms2_);
            std::cout << ' ' << v.initialized_ << '\n';
        }
    }
}
'''


class NativeTests(unittest.TestCase):
    def test_public_first_sample_update_matches_scalar_Eigen_profile(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        candidates=[Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3')),ROOT/'third_party/eigen']
        candidates+=list(Path('/opt/pyvenv').glob('lib/python*/site-packages/casadi/include/eigen3'))
        eigen=next((p for p in candidates if (p/'Eigen/Dense').is_file()),None)
        if not cxx or eigen is None:
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required native first-sample prerequisites missing')
            self.skipTest('native first-sample correspondence needs g++ and Eigen')
        with tempfile.TemporaryDirectory() as td:
            cpp=Path(td)/'seed.cpp'; exe=cpp.with_suffix('');cpp.write_text(HARNESS)
            subprocess.run([cxx,'-std=c++20','-O2','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_NON_ARDUINO','-DEIGEN_DONT_VECTORIZE',f'-I{ROOT / "src"}',f'-I{eigen}',
                str(cpp),'-o',str(exe)],check=True,capture_output=True,text=True,timeout=90)
            rows=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=15).stdout.splitlines()
        s=V.State()
        for row in rows:
            words=list(map(int,row.split())); k=words[0]
            if k==0: s=V.State()
            vals=list(map(M.value,words[1:-1])); gyro=tuple(vals[:3]);acc=tuple(vals[3:6])
            out=X.step(s,CFG,dt=DT,gyro=gyro,acc=acc)
            expected=V.State(tuple(vals[6:10]),tuple(vals[10:13]),bool(words[-1]),vals[13],vals[14])
            self.assertEqual(out.vertical.state,expected,msg=f'first/update row: {row}')
            if out.startup is not None: out.startup.validate()
            s=out.vertical.state


if __name__=='__main__': unittest.main()
