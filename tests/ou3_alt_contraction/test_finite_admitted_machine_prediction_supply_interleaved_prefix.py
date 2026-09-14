"""Same-event full machine prediction-supply carrier regressions."""
import unittest
from dataclasses import replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as X
import test_finite_admitted_machine_scheduler_interleaved_prefix as BASE
import test_finite_admitted_machine_prediction_interleaved_prefix as PBASE
import test_finite_source_bound_live_word as LBASE


def executed_supply(*,pending):
    import test_finite_admitted_machine_tunestate_interleaved_prefix as T
    from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
    from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as MF
    from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as BOUND
    from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
    s=T.state(usable=True,pending=pending)
    if pending:
        # Keep old active Sigma; change only the pending candidate. This must
        # change the exact prediction roots on this same IMU, not the next one.
        prefix=s.base.base.prefix; word=prefix.prefix.live.live_word; imu=word.live.live.live
        imu=replace(imu,tuner=replace(imu.tuner,tune=replace(imu.tuner.tune,sigma_applied=F(5,4))))
        word=replace(word,live=replace(word.live,live=replace(word.live.live,live=imu)))
        admitted=replace(prefix.prefix.live,live_word=word)
        s=replace(s,base=replace(s.base,base=replace(s.base.base,prefix=replace(prefix,prefix=replace(prefix.prefix,live=admitted)))))
    kwargs,machine=T.event_operands(s)
    runtime=s.base.base.prefix.prefix.live.live_word.runtime
    floors=MF.boundary_floors(s.frontends,bench_noise_sigma=runtime.boundary_bench_noise_sigma,required=pending)
    transaction=BOUND.imu_boundary(s.machine,runtime.commit_cfg,live=True,
        separate_band_noise_floor_sigma=None if not pending else floors[0].noise_sigma,
        fma_band_noise_floor_sigma=None if not pending else floors[1].noise_sigma)
    roots={}
    for mode in ('separate','fma'):
        active=getattr(s,mode+'_active') if not pending else ACTIVE.ActiveParameters.from_commit(getattr(transaction,mode+'_commit'))
        # Conditional exact-root outer witnesses, not target-libm evidence.
        x=F(1,200)/active.tau
        roots[mode+'_machine_root_kwargs']=dict(ou_alpha=1-x,ou_em1=-x,
            qaxis_marginal_psd=kwargs['qaxis_marginal_psd'],qaxis_final_psd=kwargs['qaxis_final_psd'])
    pred=PBASE.X.begin(s)
    sched=BASE.X.begin(pred,POST.Scheduler(s.separate_active.pseudo_period,0),POST.Scheduler(s.fma_active.pseudo_period,0))
    before=X.begin(sched)
    return before,X.imu_step(before,**kwargs,**machine,**roots)


class Tests(unittest.TestCase):
    def test_executed_nonpending_event_carries_frontend_through_full_prediction_supply(self):
        before,out=executed_supply(pending=False)
        self.assertEqual(out.state.supply_steps,1)
        self.assertEqual(out.state.base.scheduler_steps,1)
        self.assertEqual(out.state.base.base.prediction_steps,1)
        self.assertEqual(out.state.base.base.base.frontends.samples,1)
        self.assertFalse(out.lower.separate.retargeted)
        self.assertEqual(out.separate.supply.z,(0,)*24)
        self.assertIs(out.lower.lower.lower.separate_frontend.before,before.base.base.base.frontends.separate)

    def test_executed_changed_pending_sigma_rebuilds_postcommit_exact_prediction(self):
        before,out=executed_supply(pending=True)
        previous=before.base.base.base
        self.assertNotEqual(out.separate.exact_roots.qaxis.sigma_aw,previous.separate_active.Sigma_aw)
        self.assertEqual(out.separate.exact_roots.qaxis.sigma_aw[2][2],F(25,16))
        self.assertEqual(out.separate.machine_roots.active,out.state.base.base.base.separate_active)
        self.assertNotEqual(out.separate.supply.covariance,((0,)*21,)*21)
        self.assertEqual(len(out.separate.supply.z),24)
        self.assertEqual(len(out.separate.supply.covariance),21)
        self.assertTrue(out.lower.separate.retargeted)
        self.assertTrue(out.lower.fma.retargeted)
        self.assertEqual(out.state.supply_steps,1)
        self.assertEqual(out.state.base.base.base.frontends.samples,1)

    def test_begin_roots_supply_counter_at_actual_admitted_ordinal(self):
        s=X.begin(BASE.state())
        self.assertEqual(s.supply_steps,0)
        self.assertEqual(X._imu_steps(s.base),s.entry_imu_steps)

    def test_counter_cannot_detach_from_lower_word(self):
        base=BASE.state(); s=X.begin(base)
        with self.assertRaisesRegex(ValueError,'count detached'):
            X.State(base,s.entry_imu_steps,1)

    def test_exact_root_witness_extractor_fails_closed(self):
        with self.assertRaisesRegex(TypeError,'witnesses missing'):
            X._exact_root_kwargs({'ou_alpha':1})

    def test_MAG_and_HOLD_preserve_prediction_supply_count(self):
        s=X.begin(BASE.state()); n=s.supply_steps
        word=PBASE.X._word(s.base.base.base)
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertEqual(s.supply_steps,n)
        s,_=X.set_hold(s,hold=False)
        self.assertEqual(s.supply_steps,n)

    def test_complete_cannot_skip_full_prediction_supply(self):
        with self.assertRaises((ValueError,TypeError)):
            X.complete(X.begin(BASE.state()))

    def test_readiness_stops_before_postprediction_measurement_and_storage(self):
        r=X.readiness()
        self.assertTrue(r['machine_root_effect_injected_into_joint24_prediction_relation'])
        self.assertTrue(r['rebuilt_exact_prediction_required_equal_executed_shipping_prediction'])
        self.assertTrue(r['complete_word_requires_full_prediction_supply_on_all_600_IMU_edges'])
        for k in ('post_prediction_floor_and_S_service_propagate_machine_prediction_supply',
                  'machine_RS_measurement_effect_attached','accelerometer_measurement_propagates_machine_prediction_supply',
                  'source_uniform_machine_prediction_supply_bound_closed','all_target_libm_and_Eigen_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
