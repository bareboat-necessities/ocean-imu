from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_handoff_control as X
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as BRIDGE
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AW
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
import test_finite_startup_live_runtime_bridge as FIXTURE


class Tests(unittest.TestCase):
    def decision(self, **changes):
        d = X.Decision(X.Config(), CLOCK.TIMEOUT_CROSSING, True, True,
                       'TunerWarm', True, F(0), False)
        return replace(d, **changes)

    def test_timeout_retains_gravity_and_proxy_but_does_not_require_tuner_or_north(self):
        for stage in ('Cold', 'TunerWarm', 'TunerReady'):
            d = self.decision(tuner_stage=stage)
            self.assertTrue(d.by_timeout)
            self.assertFalse(d.by_quality)
        for key in ('begun', 'proxy_initialized', 'aligned_branch'):
            self.assertFalse(self.decision(**{key: False}).handoff_due)
        self.assertFalse(self.decision(wrapper_time=CLOCK.TIMEOUT_PREDECESSOR).handoff_due)

    def test_quality_and_extended_magnetic_deadline_follow_literal_predicates(self):
        d = self.decision(wrapper_time=F(8), tuner_stage='TunerReady',
                          gravity_good=F(2), north_reference=True)
        self.assertTrue(d.by_quality)
        self.assertFalse(d.by_timeout)
        self.assertFalse(replace(d, tuner_stage='TunerWarm').handoff_due)
        self.assertFalse(replace(d, north_reference=False).handoff_due)
        self.assertTrue(replace(d, config=X.Config(with_mag=False), north_reference=False).by_quality)
        self.assertEqual(X.Config().acquisition_deadline, 60)
        extended = X.Config(settle=F(140))
        self.assertEqual(extended.deadline, 200)
        self.assertFalse(self.decision(config=extended).handoff_due)

    def test_timeout_bridge_preserves_memory_and_rejects_detached_control(self):
        front, entry, fresh, active, scheduler = FIXTURE.objects(False)
        t = CLOCK.STARTUP_TIMEOUT_STEPS * CLOCK.DT_REAL
        front = replace(front, tuner=replace(front.tuner, stage='TunerWarm', time=t,
                        vertical=replace(front.tuner.vertical, initialized=True)))
        ref = replace(fresh.reference, time=t, live_origin=t)
        fresh = replace(fresh, reference=ref)
        kwargs = dict(scope=SCOPE.certified_scope(), commit_cfg=FIXTURE.cfg(),
                      bench_noise_sigma=0, noise_sqrt=BAND.NoiseSqrtWitness(0),
                      scheduler=scheduler, racc=RACC.State(), aw_sync=AW.State())
        with self.assertRaisesRegex(ValueError, 'attached timeout'):
            BRIDGE.bridge(entry, fresh, front, **kwargs)
        d = self.decision()
        out = BRIDGE.bridge(entry, fresh, front, handoff_decision=d, **kwargs)
        self.assertIs(out.frontend_live.tuner.vertical, front.tuner.vertical)
        self.assertIs(out.frontend_live.tuner.wpe, front.tuner.wpe)
        self.assertEqual(out.frontend_live.tuner.stage, 'Live')
        self.assertEqual(out.state.aw_sync.last_sync_time, t)
        self.assertFalse(out.frontend_live.tuner.pending)
        for bad in (replace(d, tuner_stage='Cold'), replace(d, wrapper_time=F(151))):
            with self.assertRaisesRegex(ValueError, 'detached'):
                BRIDGE.bridge(entry, fresh, front, handoff_decision=bad, **kwargs)

    def test_clock_crossing_and_all_machine_budgets_include_last_two_samples(self):
        from tools.stability.ou3_alt_contraction import finite_tuner_tau_600_domain as TAU
        from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TL
        from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SL
        from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RL
        from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as BANDL
        from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as STATS
        from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
        from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
        from tools.stability.ou3_alt_contraction import finite_aw_sync_clock_binary64 as AWC
        self.assertEqual(CLOCK.STARTUP_TIMEOUT_STEPS, 30002)
        self.assertLess(CLOCK.TIMEOUT_PREDECESSOR, 150)
        self.assertGreaterEqual(CLOCK.TIMEOUT_CROSSING, 150)
        for count in (TAU.TOTAL_UPDATES, TL.MAX_UPDATES, SL.MAX_UPDATES, RL.MAX_UPDATES,
                      BANDL.MAX_SAMPLES, STATS.MAX_SAMPLES, STILL.MAX_SAMPLES,
                      GUARD.MAX_SAMPLES, AWC.MAX_STEPS):
            self.assertEqual(count, 30602)
        self.assertTrue(TAU.assert_induction_closes())
        self.assertFalse(CLOCK.readiness()['universal_startup_deadline_closed'])


if __name__ == '__main__':
    unittest.main()
