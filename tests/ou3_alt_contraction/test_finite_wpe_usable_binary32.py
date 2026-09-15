from fractions import Fraction as F
import unittest
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
from tools.stability.ou3_alt_contraction import finite_wpe_usable_binary32 as U
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M


class Tests(unittest.TestCase):
    def test_both_inclusive_boundaries_are_required(self):
        w=U.PeriodWitness(0,1)
        kw=dict(produced_period=True,log_period=0,lambda_=1,witness=w)
        self.assertFalse(U.update(False,elapsed=U.B.rn32(F(399,100)),**kw).after)
        r=U.update(False,elapsed=4,**kw)
        self.assertTrue(r.after)
        self.assertEqual((r.moment_start,r.usable_floor,r.moment_history),(3,4,1))
        self.assertFalse(U.update(False,elapsed=4,produced_period=True,log_period=1,
            lambda_=1,witness=U.PeriodWitness(1,2)).after)
        self.assertTrue(U.update(False,elapsed=5,produced_period=True,log_period=1,
            lambda_=1,witness=U.PeriodWitness(1,2)).after)

    def test_already_latched_and_nonproducing_paths_consume_no_getter(self):
        for before in (False,True):
            for produced in (False,True):
                if produced and not before: continue
                kw=dict(produced_period=produced,log_period=None,elapsed=0,lambda_=1)
                self.assertEqual(U.update(before,**kw).after,before)
                with self.assertRaises(ValueError):
                    U.update(before,witness=U.PeriodWitness(0,1),**kw)

    def test_post_update_log_ancestry_and_invalid_period_branches(self):
        with self.assertRaisesRegex(ValueError,'post-update log'):
            U.update(False,produced_period=True,log_period=1,elapsed=5,lambda_=1,
                witness=U.PeriodWitness(0,1))
        for period in (None,0,-1):
            r=U.update(False,produced_period=True,log_period=0,elapsed=5,lambda_=1,
                witness=U.PeriodWitness(0,period))
            self.assertFalse(r.after)
            self.assertEqual(r.branch,'invalid-period')

    def test_second_moment_retains_literal_multiply_order(self):
        a=F(8589935,34359738368); p=F(8589935,8589934592); v=F(9369095,8388608)
        literal=M.B.add(M.B.mul(M.B.sub(1,a),p),M.B.mul(M.B.mul(a,v),v))
        self.assertEqual(literal,F(2816655,2147483648))
        self.assertIn(literal,M._second_moment(p,v,a))
        self.assertNotIn(literal,M._ema(p,M.B.mul(v,v),a))

    def test_latch_relation_does_not_promote_libm_or_source_deadline(self):
        r=U.readiness()
        self.assertTrue(r['post_update_binary32_usability_predicates_materialized'])
        self.assertFalse(r['target_period_exp_correspondence_closed'])
        self.assertFalse(r['source_uniform_takeover_deadline_closed'])


class NativeTests(unittest.TestCase):
    def test_actual_header_latch_and_second_moment_order(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        if not cxx:
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++ unavailable')
            self.skipTest('g++ unavailable')
        root=Path(__file__).resolve().parents[2]
        source=r"""
#include <algorithm>
#include <cmath>
#include <limits>
#include <bit>
#include <cstdint>
#include <iostream>
#define private public
#include "tuner/WavePeriodEstimator.h"
#undef private
unsigned bits(float x) { return std::bit_cast<uint32_t>(x); }
int main() {
    for (float t : {std::nextafter(4.0f,0.0f),4.0f,5.0f}) {
        for (float lp : {0.0f, std::log(2.0f), NAN}) {
            for (bool before : {false,true}) {
                WavePeriodEstimator w;
                w.lambda_=1.0f; w.elapsed_sec_=t;
                w.log_period_sec_=lp; w.usable_period_=before;
                float p=w.getPeriodSec();
                w.update_usable_period_();
                std::cout << "L " << bits(t) << " " << bits(lp) << " "
                          << bits(p) << " " << before << " " << w.usable_period_ << "\n";
            }
        }
    }
    for (int k=1;k<=200;++k) {
        WavePeriodEstimator w;
        w.lambda_=1.0f; w.elapsed_sec_=4.0f; w.weight_=1.0f;
        w.velocity_=float(k)/77.0f; w.elevation_=0.25f;
        w.velocity_sq_=0.001f;
        float previous=w.velocity_sq_;
        w.update(0.005f,0.0f);
        float alpha=1.0f-std::exp(-0.005f/w.last_moment_horizon_sec_);
        std::cout << "M " << bits(previous) << " " << bits(w.velocity_) << " "
                  << bits(alpha) << " " << bits(w.velocity_sq_) << "\n";
    }
}
"""
        with tempfile.TemporaryDirectory() as td:
            cpp=Path(td)/'wpe.cpp'; exe=Path(td)/'wpe'
            cpp.write_text(source)
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                f'-I{root / "src"}',str(cpp),'-o',str(exe)],check=True,
                capture_output=True,text=True,timeout=120)
            result=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=30)
        def decode(bits):
            value=struct.unpack('!f',struct.pack('!I',int(bits)))[0]
            return F(value) if value==value and abs(value)!=float('inf') else None
        mismatches=0; latch_cases=0; moment_cases=0
        for line in result.stdout.splitlines():
            kind,*fields=line.split()
            if kind=='L':
                t,lp,period=map(decode,fields[:3]); before,after=map(int,fields[3:])
                decision=U.update(bool(before),produced_period=True,log_period=lp,
                    elapsed=t,lambda_=1,witness=None if before else U.PeriodWitness(lp,period))
                self.assertEqual(decision.after,bool(after)); latch_cases+=1
            else:
                previous,value,alpha,actual=map(decode,fields)
                self.assertIn(actual,M._second_moment(previous,value,alpha))
                mismatches+=actual not in M._ema(previous,M.B.mul(value,value),alpha)
                moment_cases+=1
        self.assertEqual((latch_cases,moment_cases),(18,200))
        self.assertGreater(mismatches,0,'native cases must discriminate the wrong multiplication order')


if __name__=='__main__': unittest.main()
