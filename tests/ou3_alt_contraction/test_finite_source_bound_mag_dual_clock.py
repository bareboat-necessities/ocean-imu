import unittest
from dataclasses import replace

from tools.stability.ou3_alt_contraction import finite_source_bound_mag_dual_clock as X
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as WC
import test_finite_source_bound_live_word as BASE


class Tests(unittest.TestCase):
    def test_saturated_machine_count_keeps_exact_total_event_ledger(self):
        s=BASE.root_state(); inter=s.live
        maximum=X.LIVE.GATE.SIGNED_COUNTER_MAX
        t=inter.live.live.mekf.reference.time
        control=replace(inter.magnetic.control,updates=maximum,first_time=t)
        inter=replace(inter,magnetic=replace(inter.magnetic,control=control),
            clock=replace(inter.clock,calls=maximum,last_time=t))
        s=replace(s,live=inter)
        out=X.mag_step(s,**BASE.mag_kwargs(s))
        self.assertEqual(out.state.live.clock.calls,maximum+1)
        self.assertEqual(out.state.live.magnetic.control.updates,maximum)
        self.assertIsNotNone(out.forcing)
        after=BASE.imu(out.state)
        self.assertEqual(after.state.live.clock.calls,maximum+1)
        self.assertEqual(after.state.live.magnetic.control.updates,maximum)

    def test_sample_zero_dual_clock_mag_preserves_source_and_forcing(self):
        s=BASE.root_state(); before=s.source
        out=X.mag_step(s,**BASE.mag_kwargs(s))
        self.assertEqual(out.state.source,before)
        self.assertEqual(out.state.live.live.live.mekf.reference,
                         s.live.live.live.mekf.reference)
        self.assertIsNotNone(out.forcing)

    def test_physical_schedule_clock_and_wrapper_binary32_clock_are_not_identified(self):
        s=BASE.root_state(); ref=s.live.live.live.mekf.reference
        ts=WC.at_physical_time(ref.time)
        self.assertEqual(s.live.clock.live_time,ref.live_origin)
        self.assertEqual(ts.physical_time,ref.time)
        if ref.time:
            self.assertNotEqual(ts.wrapper_time,ref.time)

    def test_readiness_closes_live_edge_not_startup_history(self):
        r=X.readiness()
        self.assertTrue(r['source_checked_async_endpoint_required_before_dual_clock_magnetic_edge'])
        self.assertTrue(r['outer_binary32_and_inner_physical_magnetic_clocks_composed_on_Live_edge'])
        self.assertTrue(r['physical_MAG_CALL_SCHEDULE_clock_kept_distinct_from_wrapper_binary32_clock'])
        self.assertTrue(r['same_event_correlated_magnetic_ISS_forcing_retained'])
        self.assertFalse(r['startup_dual_clock_history_feeds_this_master_edge'])
        self.assertFalse(r['binary32_exp_solver_roundoff_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
