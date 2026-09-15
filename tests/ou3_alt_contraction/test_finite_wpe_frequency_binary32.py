"""WPE binary32 getter/tuner-store deployment regressions."""
from fractions import Fraction as F
import unittest
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
from dataclasses import replace

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as STORE
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as X
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE


class Tests(unittest.TestCase):
    def test_eager_machine_getter_and_exact_shadow_select_independent_branches(self):
        for exact_usable in (False,True):
            shadow=WPE.WPEState(log_period=F(0) if exact_usable else None,usable_period=exact_usable)
            for machine_usable in (False,True):
                for returned in (None,0,-1,B.rn32(F(1,2))):
                    out=X.machine_frequency(shadow,log_period=B.rn32(0),usable=machine_usable,
                        getter=X.FrequencyExp(0,returned),min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
                        shadow_frequency=F(1) if exact_usable else None)
                    expected=returned if machine_usable and returned is not None and returned>0 else X.PRIOR
                    self.assertEqual(out.stored.input_hz,expected)
                    self.assertEqual(out.exact_shadow_frequency,F(1) if exact_usable else F(1,5))
                    self.assertEqual(out.machine_minus_shadow,out.stored.stored_hz-out.exact_clamped_frequency)
                    self.assertEqual(out.machine_read.exp.result,returned)
        # A prior-selected result does not erase the eager getter's ancestry.
        with self.assertRaisesRegex(ValueError,'eager frequency exp detached'):
            X.machine_frequency(WPE.WPEState(),log_period=0,usable=False,getter=None,min_hz=1,max_hz=2)
        with self.assertRaisesRegex(ValueError,'eager frequency exp detached'):
            X.MachineRead(0,False,X.FrequencyExp(1,1))
        self.assertEqual(X.MachineRead(None,True,None).frequency,X.PRIOR)

    def test_outer_supply_includes_nonfinite_and_retained_state(self):
        lo,hi=F(3,100),F(6,5)
        bound=X.outer_supply_bound(exact_min_hz=lo,exact_max_hz=hi)
        ml,mh=map(B.rn32,(lo,hi))
        for exact in (None,-1,0,F(1,5),100):
            for machine in (None,-1,0,B.rn32(F(1,5)),100):
                e=X.final_tuning_clamp(exact,min_hz=lo,max_hz=hi)
                m=X.final_tuning_clamp(machine,min_hz=ml,max_hz=mh)
                self.assertEqual(bound.check(e,m),m-e)
        self.assertEqual(X.final_tuning_clamp(None,min_hz=ml,max_hz=mh),ml)


    def test_uniform_supply_includes_rounded_endpoints_and_opposite_branches(self):
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
        cfg=R.StatsConfig(4,F(3,10),60,F(1,20),5)
        bound=X.statistics_supply_bound(cfg,exact_min_hz=F(3,100),exact_max_hz=F(6,5))
        self.assertEqual(bound.exact_interval,(F(1,20),F(6,5)))
        self.assertEqual(bound.machine_interval,(B.rn32(F(1,20)),B.rn32(F(6,5))))
        lo,hi=bound.residual_interval
        self.assertLess(lo,0); self.assertGreater(hi,F(23,20))
        for e in bound.exact_interval:
            for m in bound.machine_interval:
                self.assertEqual(bound.check(e,m),m-e)
        with self.assertRaises(ValueError): bound.check(6,B.rn32(F(1,5)))
        with self.assertRaises(ValueError): bound.check(F(1,5),6)

    def test_uniform_supply_handles_disjoint_clamp_ranges(self):
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
        cfg=R.StatsConfig(4,F(3,10),60,2,5)
        b=X.statistics_supply_bound(cfg,exact_min_hz=F(1,10),exact_max_hz=1)
        self.assertEqual(b.exact_interval,(1,1))
        self.assertEqual(b.machine_interval,(1,1))
        self.assertEqual(b.residual_interval,(0,0))
    def test_statistics_clamp_precedes_distinct_outer_tuning_clamp(self):
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
        from dataclasses import replace
        cfg=R.StatsConfig(4,F(3,10),60,F(1,20),5)
        shadow=self.shadow(); log=X.bind_log_state(shadow,B.rn32(shadow.log_period))
        # The getter value remains a conditional libm witness. This tests the
        # literal downstream clamp graph, not exp(log) accuracy.
        getter=X.getters(log,period_exp=25,frequency_exp=B.rn32(F(1,25)))
        external=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
            getter=getter,shadow_frequency=F(1,25))
        out=X.through_statistics(external,cfg,exact_min_hz=F(3,100),exact_max_hz=F(6,5))
        self.assertEqual(out.stats_stored.stored_hz,B.rn32(cfg.f_min))
        self.assertEqual(out.stored.stored_hz,B.rn32(cfg.f_min))
        self.assertNotEqual(out.stored.stored_hz,external.stored.stored_hz)
        self.assertEqual(out.exact_clamped_frequency,cfg.f_min)
        self.assertEqual(out.machine_minus_shadow,B.rn32(cfg.f_min)-cfg.f_min)
        self.assertEqual(out.supply_bound.check(out.exact_clamped_frequency,out.stored.stored_hz),out.machine_minus_shadow)
        with self.assertRaisesRegex(ValueError,'ordered statistics'):
            replace(out,stored=external.stored)

    def test_stats_upper_and_tuning_upper_bounds_are_not_identified(self):
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
        shadow=self.shadow(); log=X.bind_log_state(shadow,B.rn32(shadow.log_period))
        external=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
            getter=X.getters(log,period_exp=1,frequency_exp=8),shadow_frequency=8)
        out=X.through_statistics(external,R.StatsConfig(4,F(3,10),60,F(1,20),5),exact_min_hz=F(3,100),exact_max_hz=F(6,5))
        self.assertEqual(out.stats_stored.input_hz,8)
        self.assertEqual(out.stats_stored.stored_hz,5)
        self.assertEqual(out.stored.input_hz,5)
        self.assertEqual(out.stored.stored_hz,B.rn32(F(6,5)))
        self.assertEqual(out.exact_clamped_frequency,F(6,5))
        self.assertNotEqual(out.machine_minus_shadow,0)
        from dataclasses import replace
        with self.assertRaisesRegex(ValueError,'outer tuning bounds detached'):
            replace(out,exact_tune_bounds=(F(3,100),F(7,5)))

    def shadow(self): return WPE.WPEState(log_period=F(7,10),usable_period=True)
    def exact_frequency(self): return F(1,2)
    def getter(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        return X.getters(log,period_exp=B.rn32(F(201,100)),frequency_exp=B.rn32(F(497,1000)))

    def test_two_shipping_exp_calls_share_one_stored_log_but_not_bit_reciprocity(self):
        out=self.getter()
        self.assertEqual(out.period_argument,out.log.stored_log_period)
        self.assertEqual(out.frequency_argument,-out.log.stored_log_period)
        self.assertNotEqual(out.period_result*out.frequency_result,1)

    def test_usable_preupdate_WPE_retains_machine_minus_exact_frequency_supply(self):
        shadow=self.shadow(); out=self.getter()
        q=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
                            getter=out,shadow_frequency=self.exact_frequency())
        self.assertEqual(q.branch,'wpe'); self.assertIsInstance(q.stored,STORE.StoredFrequency)
        self.assertEqual(q.stored.input_hz,out.frequency_result)
        self.assertEqual(q.exact_clamped_frequency,self.exact_frequency())
        self.assertEqual(q.machine_minus_shadow,q.stored.stored_hz-q.exact_clamped_frequency)
        self.assertNotEqual(q.machine_minus_shadow,0)

    def test_preusable_WPE_uses_literal_prior_and_retains_constant_rounding_supply(self):
        shadow=WPE.WPEState()
        q=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)))
        self.assertEqual(q.branch,'prior'); self.assertIsNone(q.getter)
        self.assertEqual(q.exact_shadow_frequency,X.PRIOR_EXACT)
        self.assertEqual(q.stored.input_hz,X.PRIOR)
        self.assertEqual(q.machine_minus_shadow,q.stored.stored_hz-X.PRIOR_EXACT)
        with self.assertRaisesRegex(ValueError,'consumes no getter'):
            X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
                              getter=self.getter(),shadow_frequency=self.exact_frequency())

    def test_log_state_and_getter_arguments_cannot_be_spliced(self):
        shadow=self.shadow()
        with self.assertRaisesRegex(ValueError,'deployment residual detached'):
            X.StoredLogPeriod(shadow.log_period,B.rn32(F(7,10)),0)
        log=X.bind_log_state(shadow,B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'arguments detached'):
            X.GetterResult(log,log.stored_log_period,log.stored_log_period,B.rn32(2),B.rn32(F(1,2)))
        other=WPE.WPEState(log_period=F(4,5),usable_period=True)
        with self.assertRaisesRegex(ValueError,'preupdate canonical'):
            X.tuner_frequency(other,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),
                              getter=self.getter(),shadow_frequency=self.exact_frequency())

    def test_getter_requires_binary32_results_without_claiming_libm_correctness(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'binary32 witnesses'):
            X.getters(log,period_exp=F(2),frequency_exp=F(1,3))
        r=X.readiness()
        self.assertTrue(r['shipping_WPE_dual_exp_and_tuner_source_shape_matches'])
        self.assertTrue(r['machine_minus_exact_tuner_frequency_supply_exposed'])
        self.assertTrue(r['WPE_frequency_getter_to_tuner_binary32_store_topology_closed'])
        self.assertFalse(r['binary32_period_frequency_bit_reciprocity_assumed'])
        self.assertFalse(r['WPE_binary32_log_period_production_closed'])
        self.assertFalse(r['WPE_frequency_exp_target_libm_correspondence_closed'])
        self.assertTrue(r['source_uniform_WPE_frequency_supply_bound_closed'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_LIVE_PASS'])


class NativeTests(unittest.TestCase):
    def test_actual_wrapper_selection_and_statistics_clamp(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++/Eigen unavailable')
            self.skipTest('g++/Eigen unavailable')
        root=Path(__file__).resolve().parents[2]
        source=r"""
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>
#include <bit>
#include <cstdint>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std=9.80665f;
unsigned bits(float x) { return std::bit_cast<uint32_t>(x); }
int main() {
    SeaStateFusionFilter_OU_III<TrackerType::KALMANF> f;
    for (float lp : {NAN,INFINITY,-INFINITY,0.0f,std::log(2.0f),100.0f,-100.0f,1000.0f}) {
        for (bool usable : {false,true}) {
            f.wave_period_.log_period_sec_=lp;
            f.wave_period_.usable_period_=usable;
            const float returned=f.wave_period_.getFrequencyHz();
            const float chosen=f.tuner_frequency_hz_();
            f.tuner_.update(0.005f,0.0f,chosen);
            std::cout << bits(lp) << " " << usable << " " << bits(returned) << " "
                      << bits(chosen) << " " << bits(f.tuner_.getFrequencyHz()) << "\n";
        }
    }
}
"""
        with tempfile.TemporaryDirectory() as td:
            cpp=Path(td)/'frequency.cpp';exe=Path(td)/'frequency';cpp.write_text(source)
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_DONT_VECTORIZE',f'-I{eigen}',f'-I{root / "src"}',str(cpp),'-o',str(exe)],
                check=True,capture_output=True,text=True,timeout=180)
            run=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=30)
        def decode(bits):
            value=struct.unpack('!f',struct.pack('!I',int(bits)))[0]
            return F(value) if value==value and abs(value)!=float('inf') else None
        subnormal=False
        rows=run.stdout.splitlines();self.assertEqual(len(rows),16)
        for row in rows:
            lb,ub,rb,cb,sb=row.split();lp,returned,chosen,stored=map(decode,(lb,rb,cb,sb))
            exp=None if lp is None else X.FrequencyExp(-lp,returned)
            read=X.MachineRead(lp,bool(int(ub)),exp)
            self.assertEqual(read.frequency,chosen)
            clamp=STORE.store(chosen,B.rn32(F(1,20)),B.rn32(5))
            self.assertEqual(clamp.stored_hz,stored)
            subnormal |= chosen>0 and chosen<F(1,1<<126)
        self.assertTrue(subnormal,'native input class must include a selected subnormal getter')


if __name__=='__main__': unittest.main()
