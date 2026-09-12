from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as X

class Tests(unittest.TestCase):
    def test_default_25hz_schedule_forces_250_call_guard_within_10s(self):
        r=X.release_reachability()
        self.assertEqual(r.unlock_count,250)
        self.assertEqual(r.first_call_after_live_max,F(1,25))
        self.assertEqual(r.elapsed_first_to_unlock_max,F(249,25))
        self.assertEqual(r.live_to_unlock_max,F(10))
        self.assertTrue(r.strict_one_second_guard_satisfied)
        self.assertTrue(X.require_release_reachable(r))

    def test_schedule_is_explicit_assumption_and_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'wrong magnetometer'): X.Schedule(F(1,25),F(1,25),'OTHER')
        fast=X.Schedule(F(1,2),F(1,1000))
        r=X.release_reachability(fast,2)
        self.assertEqual(r.elapsed_first_to_unlock_max,F(1001,1000))
        self.assertEqual(r.live_to_unlock_max,F(1501,1000))
        self.assertTrue(X.require_release_reachable(r))
        for n in (0,True,F(3,2)):
            with self.assertRaisesRegex(ValueError,'positive integer'):
                X.release_reachability(unlock_count=n)

    def test_clustered_250_calls_do_not_fake_the_strict_time_guard(self):
        from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as G
        import test_finite_core as CORE
        core=CORE.root('H')
        cfg=G.Config(mag_delay=0)
        control=G.State()
        first=F(1,100)
        # All 250 calls occur within 0.249 s, although every gap is <=40 ms.
        for j in range(250):
            out=G.update_mag_call(control,core,cfg,time=first+F(j,1000),
                                  live=True,measurement_state=core)
            control,core=out.state,out.filter_state
        self.assertTrue(control.locked)
        # Continue the schedule to the strict boundary without missing a gap.
        for j in range(250,1001):
            out=G.update_mag_call(control,core,cfg,time=first+F(j,1000),
                                  live=True,measurement_state=core)
            control,core=out.state,out.filter_state
        self.assertTrue(control.locked)  # exactly first+1 is NOT enough
        out=G.update_mag_call(control,core,cfg,time=first+F(1001,1000),
                              live=True,measurement_state=core)
        self.assertFalse(out.state.locked)
        self.assertEqual(out.filter_state.mode,'A')
        self.assertGreater(out.state.updates,250)
        self.assertLessEqual(first+F(1001,1000),X.release_reachability().live_to_unlock_max)

    def test_prefix_deadlines_and_equal_timestamp_calls(self):
        s=X.default_schedule();c=X.Clock(F(7))
        self.assertTrue(X.check_prefix(c,s,time=F(176,25)))
        c=X.record_call(c,s,time=F(176,25))
        c=X.record_call(c,s,time=c.last_time)
        self.assertEqual(c.calls,2)
        with self.assertRaisesRegex(ValueError,'deadline missed'):
            X.check_prefix(c,s,time=c.last_time+s.gap_max+F(1,1000))
        with self.assertRaisesRegex(ValueError,'backwards'):
            X.record_call(c,s,time=c.last_time-F(1,1000))
        with self.assertRaisesRegex(ValueError,'must agree'):
            X.Clock(7,calls=1)

    def test_current_schedule_does_not_prove_signed_counter_lifetime(self):
        c=X.counter_lifetime()
        self.assertEqual(c.horizon,F(3))
        self.assertEqual(c.signed_max,(1<<31)-1)
        self.assertFalse(c.schedule_supplies_positive_min_gap)
        self.assertIsNone(c.uniform_call_count_upper)
        self.assertFalse(c.no_signed_overflow_proved)
        with self.assertRaisesRegex(ValueError,'does not bound signed counter lifetime'):
            X.require_counter_lifetime_closed(c)
        # Equal-timestamp calls are legal in V1, so no rate ceiling can be
        # inferred from the 40 ms maximum-gap requirement.
        s=X.default_schedule(); clock=X.Clock(0)
        for _ in range(1000):
            clock=X.record_call(clock,s,time=0)
        self.assertEqual(clock.calls,1000)

    def test_readiness_does_not_claim_complete_edge(self):
        r=X.readiness()
        self.assertTrue(r['named_deterministic_mag_call_schedule_declared'])
        self.assertTrue(r['strict_one_second_guard_forced'])
        self.assertFalse(r['schedule_supplies_positive_minimum_call_gap'])
        self.assertIsNone(r['uniform_call_count_upper_on_canonical_3s_window'])
        self.assertFalse(r['shipping_signed_mag_counter_lifetime_closed'])
        self.assertFalse(r['external_acc_bias_hold_excluded_here'])
        self.assertFalse(r['H18_A21_complete_word_edge_attached'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
