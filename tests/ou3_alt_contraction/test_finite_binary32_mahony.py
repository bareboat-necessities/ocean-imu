"""Exact-rounding regressions and passive shipping correspondence, not admission."""
from dataclasses import replace
from fractions import Fraction as F
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as B
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


HARNESS = r'''
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
// Observe private state only; never assign synthetic runtime roots.
#define private public
#include "tuner/VerticalAccelComplementary.h"
#undef private
uint32_t bits(float f) { uint32_t b; std::memcpy(&b,&f,4); return b; }
float val(uint32_t b) { float f; std::memcpy(&f,&b,4); return f; }
void emit(float x) { std::cout << ' ' << bits(x); }
void state(const VerticalAccelComplementary& v) {
    emit(v.ahrs_.q0); emit(v.ahrs_.q1); emit(v.ahrs_.q2); emit(v.ahrs_.q3);
    emit(v.ahrs_.integralFBx); emit(v.ahrs_.integralFBy); emit(v.ahrs_.integralFBz);
    emit(v.elapsed_sec_); emit(v.up_ms2_);
}
int main() {
    for (uint32_t e=1;e<255;++e) {
        for (uint32_t m : {0U,1U,0x7fffffU}) {
            uint32_t b=(e<<23)|m;
            std::cout << "N " << b << ' ' << bits(Mahony_AHRS<float>::invSqrt(val(b))) << '\n';
        }
    }
    for (uint32_t b : {0U,1U,2U,0x3fffffU,0x7fffffU})
        std::cout << "N " << b << ' ' << bits(Mahony_AHRS<float>::invSqrt(val(b))) << '\n';
    for (int mode=0;mode<3;++mode) {
        float kp=0.2f, ki=mode==0 ? 0.0f : 0.02f, gravity=9.80665f, dt=0.005f;
        VerticalAccelComplementary v(kp,ki,20.0f);
        for (int k=0;k<18;++k) {
            if (mode==2 && k==8) { ki=0.0f; v.setGains(kp,ki); }
            Eigen::Vector3f gyro(0.003f*std::cos(float(k)), -0.004f, 0.002f);
            Eigen::Vector3f acc(0.4f*std::sin(0.2f*k),0.2f*std::cos(0.3f*k),-gravity+0.1f*std::sin(0.4f*k));
            // Low-level zero-feedback branch regression, not an admitted sea.
            if (k==6 || k==8) acc.setZero();
            if (k) {
                std::cout << "S " << mode << ' ' << k;
                emit(kp); emit(ki); emit(gravity); emit(dt);
                for (int i=0;i<3;++i) emit(gyro(i));
                for (int i=0;i<3;++i) emit(acc(i));
                state(v);
            }
            v.update(dt,gyro,acc,gravity);
            if (k) { state(v); std::cout << '\n'; }
        }
    }
}
'''


def quantized(x):
    return B.value(B.round_bits(x))


class ExactRoundingTests(unittest.TestCase):
    def test_decode_roundtrip_all_exponent_edges(self):
        for e in range(255):
            for m in (0,1,2,(1<<22),(1<<23)-1):
                for sign in (0,B.SIGN):
                    b=sign | (e<<23) | m
                    expected=b if b != B.SIGN else 0
                    self.assertEqual(B.round_bits(B.value(b)),expected)

    def test_ties_even_at_every_normal_binade_boundary(self):
        for e in range(1,254):
            for m in (0,1,(1<<23)-1):
                low=(e<<23)|m
                middle=(B.value(low)+B.value(low+1))/2
                expected=low+(low & 1)
                self.assertEqual(B.round_bits(middle),expected)
                self.assertEqual(B.round_bits(-middle),expected|B.SIGN)

    def test_independent_rounding_cells_reject_wrong_tie_and_outside_neighbor(self):
        one=B.exact_bits(1)
        middle=(B.value(one)+B.value(one+1))/2
        self.assertTrue(B.rounding_cell_contains(middle,one))
        self.assertFalse(B.rounding_cell_contains(middle,one+1))
        self.assertFalse(B.rounding_cell_contains(F(1),one+1))
        self.assertFalse(B.rounding_cell_contains(F(-1),one))
        self.assertFalse(B.rounding_cell_contains(B.pow2(128)-B.pow2(103),B.MAX_FINITE))

    def test_gradual_underflow_ties_and_signed_zero(self):
        tiny=B.pow2(-149)
        self.assertEqual(B.round_bits(tiny/2),0)
        self.assertEqual(B.round_bits(3*tiny/2),2)
        self.assertEqual(B.round_bits(-tiny/2),B.SIGN)
        a=B.Arithmetic()
        self.assertEqual(a.add(B.SIGN,B.SIGN),B.SIGN)
        self.assertEqual(a.add(B.SIGN,0),0)
        self.assertEqual(a.mul(B.SIGN,B.exact_bits(2)),B.SIGN)
        a.verify()

    def test_overflow_and_nonfinite_fail_closed(self):
        threshold=B.pow2(128)-B.pow2(103)
        self.assertEqual(B.round_bits(threshold-B.pow2(80)),B.MAX_FINITE)
        with self.assertRaises(OverflowError): B.round_bits(threshold)
        with self.assertRaises(ValueError): B.value(0x7f800000)
        with self.assertRaises(ValueError): B.value(0x7fc00000)
        with self.assertRaises(ValueError): B.Arithmetic().invsqrt(B.SIGN)

    def test_runtime_boundary_does_not_silently_quantize(self):
        with self.assertRaises(ValueError): B.exact_bits(F(1,10))
        with self.assertRaises(TypeError): B.exact_bits(0.1)
        with self.assertRaises(TypeError): B.exact_bits(True)

    def test_normalizer_is_not_ideal_sqrt_and_keeps_zero_subnormal(self):
        a=B.Arithmetic(); out=a.invsqrt(B.exact_bits(1))
        self.assertNotEqual(B.value(out),1)
        self.assertEqual(a.normalizations[0].seed_bits,B.MAGIC-(B.exact_bits(1)>>1))
        self.assertGreater(B.value(a.invsqrt(0)),0)
        self.assertGreater(B.value(a.invsqrt(1)),0)
        a.verify()

    def test_defect_and_normalizer_mutations_rejected(self):
        a=B.Arithmetic(); a.invsqrt(B.exact_bits(3))
        self.assertTrue(any(op.defect for op in a.operations))
        bad=replace(a.operations[0],defect=a.operations[0].defect+1)
        with self.assertRaises(ValueError): bad.validate()
        a.normalizations[0]=replace(a.normalizations[0],seed_bits=1)
        with self.assertRaises(ValueError): a.verify()

    def test_profile_seed_and_free_reciprocal_rejected(self):
        cfg=V.Config(quantized(F(1,5)),quantized(F(1,50)),quantized(F(196133,20000)),20)
        kw=dict(dt=quantized(F(1,200)),gyro=(0,0,0),acc=(0,0,-cfg.gravity))
        with self.assertRaises(ValueError): B.step_initialized(V.State(),cfg,**kw)
        with self.assertRaises(ValueError): B.step_initialized(V.State(initialized=True),cfg,profile='fma',**kw)
        with self.assertRaises(TypeError): B.step_initialized(V.State(initialized=True),cfg,quat_invnorm=1,**kw)
        result=B.step_initialized(V.State(initialized=True),cfg,**kw)
        self.assertIsInstance(result.vertical,V.Result)
        B.Arithmetic(list(result.operations),list(result.normalizations)).verify()
        # Actual invSqrt leaves a norm defect; exact normalization would be wrong.
        self.assertNotEqual(result.vertical.vertical_accel,0)

    def test_existing_vertical_API_dispatches_without_free_reciprocals(self):
        cfg=V.Config(quantized(F(1,5)),quantized(F(1,50)),quantized(F(196133,20000)),20)
        state=V.State(initialized=True)
        kw=dict(dt=quantized(F(1,200)),gyro=(0,0,0),acc=(0,0,-cfg.gravity))
        expected=B.step_initialized(state,cfg,**kw).vertical
        self.assertEqual(V.step(state,cfg,arithmetic_profile=B.PROFILE,**kw),expected)
        with self.assertRaises(ValueError):
            V.step(state,cfg,arithmetic_profile=B.PROFILE,quat_invnorm=V.InvSqrtWitness(1,1),**kw)
        with self.assertRaises(ValueError):
            V.step(V.State(),cfg,arithmetic_profile=B.PROFILE,**kw)

    def test_readiness_does_not_promote_profile_to_deployment(self):
        r=B.readiness()
        self.assertTrue(r['inverse_sqrt_seed_and_Newton_operands_bound'])
        for k in ('actual_target_compiler_profile_qualified','startup_seed_qualified',
                  'nonfinite_intermediate_branches_qualified','source_uniform_word_qualified',
                  'complete_word_finite_identity','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


class ShippingCorrespondenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        candidates=[Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3')),
                    ROOT/'third_party/eigen']
        # Local packaging sometimes provides Eigen here; CI uses libeigen3-dev.
        candidates += list(Path('/opt/pyvenv').glob('lib/python*/site-packages/casadi/include/eigen3'))
        eigen=next((p for p in candidates if (p/'Eigen/Dense').is_file()),None)
        if not cxx or eigen is None:
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE') == '1':
                raise RuntimeError('required native Mahony correspondence prerequisites missing')
            raise unittest.SkipTest('shipping correspondence needs g++ and Eigen; exact algebra tests still run')
        cls.tmp=tempfile.TemporaryDirectory()
        cpp=Path(cls.tmp.name)/'observer.cpp'; binary=cpp.with_suffix('')
        cpp.write_text(HARNESS)
        command=[cxx,'-std=c++20','-O2','-ffp-contract=off','-fno-fast-math',
                 '-DEIGEN_NON_ARDUINO','-DEIGEN_DONT_VECTORIZE',f'-I{ROOT / "src"}',
                 f'-I{eigen}',str(cpp),'-o',str(binary)]
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=60)
        cls.rows=subprocess.run([str(binary)],check=True,capture_output=True,text=True,timeout=15).stdout.splitlines()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls,'tmp'): cls.tmp.cleanup()

    def test_scalar_normalizer_matches_shipping_across_every_exponent(self):
        rows=[r.split() for r in self.rows if r.startswith('N ')]
        self.assertEqual(len(rows),767)
        for _,n,expected in rows:
            self.assertEqual(B.Arithmetic().invsqrt(int(n)),int(expected),n)

    def test_initialized_observer_matches_passive_startup_execution(self):
        rows=[r.split() for r in self.rows if r.startswith('S ')]
        self.assertEqual(len(rows),51)
        for row in rows:
            mode,k=map(int,row[1:3]); vals=list(map(lambda s:B.value(int(s)),row[3:]))
            kp,ki,gravity,dt=vals[:4]; gyro=vals[4:7]; acc=vals[7:10]
            prev=vals[10:19]; expected=vals[19:28]
            state=V.State(tuple(prev[:4]),tuple(prev[4:7]),True,prev[7],prev[8])
            out=B.step_initialized(state,V.Config(kp,ki,gravity,20),dt=dt,gyro=gyro,acc=acc)
            s=out.vertical.state
            got=list(s.q)+list(s.integral)+[s.elapsed,s.up]
            self.assertEqual(got,expected,(mode,k,[(i,g,e) for i,(g,e) in enumerate(zip(got,expected)) if g!=e]))
            B.Arithmetic(list(out.operations),list(out.normalizations)).verify()


if __name__ == '__main__': unittest.main()
