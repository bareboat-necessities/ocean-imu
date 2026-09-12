from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as X
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as WC
import test_finite_live_magnetic_word as OLD


class Tests(unittest.TestCase):
    def _core_at(self,core,time):
        return replace(core,reference=replace(core.reference,time=F(time)))

    def test_startup_call_uses_wrapper_clock_for_tuner_and_continuous_timestamp(self):
        bridge=OLD.FIRST.startup()
        cfg=OLD.X.Config(gate=OLD.X.GATE.Config(mag_delay=0),
                         gravity=OLD.X.GRAVITY.Config(mag_delay=0),
                         tuner=OLD.X.TUNER.Config(min_samples=1,min_window=0),
                         refinement_start=90,refinement_window=0,
                         continuous_enabled=True)
        model=OLD.X.SOURCE.Model((30,0,0),(0,0,0),'dual-clock-model')
        proxy=bridge.frontend_before.tuner.vertical
        word=OLD.X.START.State(OLD.X.GRAVITY.State(gravity_good=2,aligned_branch=True),
                               OLD.PREFIX.State(proxy))
        s=OLD.X.begin_startup(word,cfg,model,bridge.state.mekf.reference.history_id,'dual-clock-ref')
        ts=WC.at_physical_time(bridge.state.mekf.reference.time)
        out=X.startup_call(s,bridge.state.mekf.reference,residual_body=(0,0,0),packet_id='startup',
            proxy_q_norm=OLD.X.TILT.SqrtWitness(1,1),proxy_yaw_half=OLD.zero_yaw(),
            hi_decay=OLD.HI.Decay(cfg.sample_dt,cfg.continuous.memory,1),
            mag_norm=OLD.X.TUNER.SqrtWitness(900,30),mean_norm=OLD.X.TUNER.SqrtWitness(900,30),
            horizontal_sqrt=OLD.GAUGE.HorizontalSqrt(900,30),ready_yaw_half=OLD.zero_yaw(30))
        self.assertEqual(out.state.memory.last_hi_time,ts.wrapper_time)
        self.assertEqual(out.startup.admission.wrapper_time,ts.wrapper_time)

    def test_90_second_refinement_uses_outer_binary32_not_physical_clock(self):
        # At sample zero physical and wrapper clocks agree, so this conditional
        # startup fixture is also a valid seed for the dual-clock Live branch.
        bridge,start=OLD.startup(refinement_start=90,continuous=False)
        live=OLD.X.enter_live(start)

        core90=self._core_at(bridge.state.mekf,90)
        out90=X.live_call(live,core90,bridge.state.tuner.vertical,
                          **OLD.live_kwargs(live))
        self.assertIsNone(out90.refinement)
        self.assertFalse(out90.state.refinement_started)
        self.assertLess(WC.at_physical_time(90).wrapper_time,F(90))

        core90010=self._core_at(out90.filter,F(9001,100))
        kw=OLD.live_kwargs(out90.state); kw.update(OLD.refine_kwargs())
        out90010=X.live_call(out90.state,core90010,bridge.state.tuner.vertical,**kw)
        self.assertIsNotNone(out90010.refinement)
        self.assertTrue(out90010.state.refinement_started)
        self.assertGreaterEqual(WC.at_physical_time(F(9001,100)).wrapper_time,F(90))

    def test_inner_measurement_time_remains_physical_coordinate(self):
        bridge,start=OLD.startup(refinement_start=200,continuous=False)
        live=OLD.X.enter_live(start)
        core=self._core_at(bridge.state.mekf,90)
        out=X.live_call(live,core,bridge.state.tuner.vertical,**OLD.live_kwargs(live))
        self.assertEqual(out.measurement.state.filter.reference.time,F(90))
        self.assertEqual(out.state.memory.applied.last_time,None)  # blocked before refinement completes

    def test_readiness_closes_dual_clock_topology_only(self):
        r=X.readiness()
        for key in ('startup_outer_delay_uses_exact_binary32_wrapper_clock',
                    'startup_continuous_statistics_use_binary32_wrapper_elapsed_time',
                    'startup_tuner_packet_clock_uses_binary32_wrapper_time',
                    'live_refinement_start_and_elapsed_use_binary32_wrapper_clock',
                    'live_continuous_sample_and_apply_clocks_use_binary32_wrapper_clock',
                    'inner_MEKF_magnetic_call_retains_physical_inner_time',
                    'dual_clock_magnetic_word_composed',
                    'canonical_prefix_wrapper_clock_arithmetic_closed'):
            self.assertTrue(r[key])
        self.assertFalse(r['binary32_exp_solver_roundoff_closed'])
        self.assertFalse(r['source_uniform_complete_magnetic_word_qualified'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
