"""Machine post-prediction/S/accelerometer same-event supply regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as PRED
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
import test_finite_admitted_machine_scheduler_interleaved_prefix as SBASE
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE
import test_finite_source_bound_prediction_word as PBASE
import test_finite_live_magnetic_word as MAG


def state(*,machine_due=False):
    # Executed-event fixtures require an actually usable WPE/log state.  The
    # scheduler-only helper intentionally starts pre-usable because its own unit
    # tests never execute the tuner event, so build the same stack explicitly.
    p=SBASE.PRED.begin(TBASE.state(usable=True))
    sa=p.base.separate_active; fa=p.base.fma_active
    ss=POST.Scheduler(sa.pseudo_period,F(0),F(0))
    fs=POST.Scheduler(fa.pseudo_period,F(0),F(0))
    s=SBASE.X.begin(p,ss,fs)
    if machine_due:
        h=F(1,200)
        ss=POST.Scheduler(s.separate_scheduler.period,s.separate_scheduler.period-h,s.separate_scheduler.tolerance)
        fs=POST.Scheduler(s.fma_scheduler.period,s.fma_scheduler.period-h,s.fma_scheduler.tolerance)
        s=replace(s,separate_scheduler=ss,fma_scheduler=fs)
    return X.begin(PRED.begin(s))


def event_operands(s):
    sched=s.base.base
    mtune=sched.base.base
    kwargs,machine=TBASE.event_operands(mtune)
    roots=PBASE.root_args(); roots.pop('temperature_c')
    return dict(machine,**kwargs,
        separate_machine_root_kwargs=dict(roots),fma_machine_root_kwargs=dict(roots))


class Tests(unittest.TestCase):
    def test_identical_machine_coefficients_propagate_zero_supply_through_rejected_measurements(self):
        s=state(); kw=event_operands(s)
        out=X.imu_step(s,separate_accel_ldlt=MAG.REJECT,fma_accel_ldlt=MAG.REJECT,**kw)
        self.assertEqual(out.state.measurement_steps,1)
        for mode in (out.separate,out.fma):
            self.assertFalse(mode.scheduler_event.due)
            self.assertIsNone(mode.S_service.measurement)
            self.assertFalse(mode.accelerometer.accepted)
            self.assertTrue(all(v==0 for v in mode.accel_supply.z))
            self.assertTrue(all(v==0 for row in mode.accel_supply.covariance for v in row))

    def test_machine_due_branch_is_reexecuted_even_when_exact_shadow_is_not_due(self):
        s=state(machine_due=True); kw=event_operands(s)
        out=X.imu_step(s,separate_S_ldlt=MAG.REJECT,fma_S_ldlt=MAG.REJECT,
            separate_accel_ldlt=MAG.REJECT,fma_accel_ldlt=MAG.REJECT,**kw)
        exact_live=X._live_result(out.lower)
        self.assertFalse(exact_live.post_prediction.S_service_due)
        self.assertTrue(out.separate.scheduler_event.due)
        self.assertTrue(out.fma.scheduler_event.due)
        self.assertIsNotNone(out.separate.S_service.measurement)
        self.assertFalse(out.separate.S_service.measurement.accepted)
        self.assertEqual(out.separate.post_prediction.scheduler,out.separate.scheduler_event.after_step)

    def test_due_machine_branch_cannot_borrow_missing_or_exact_ldlt_implicitly(self):
        s=state(machine_due=True); kw=event_operands(s)
        with self.assertRaisesRegex(TypeError,'due S branch requires shipping safe-LDLT witness'):
            X.imu_step(s,separate_S_ldlt=None,fma_S_ldlt=MAG.REJECT,
                separate_accel_ldlt=MAG.REJECT,fma_accel_ldlt=MAG.REJECT,**kw)

    def test_accepted_machine_only_S_service_retains_nonzero_full_covariance_supply(self):
        s=state(machine_due=True); kw=event_operands(s)
        out=X.imu_step(s,separate_S_ldlt=MAG.PASS,fma_S_ldlt=MAG.PASS,
            separate_accel_ldlt=MAG.REJECT,fma_accel_ldlt=MAG.REJECT,**kw)
        exact_live=X._live_result(out.lower)
        self.assertFalse(exact_live.post_prediction.S_service_due)
        self.assertTrue(out.separate.S_service.measurement.accepted)
        self.assertTrue(any(v!=0 for row in out.separate.S_supply.covariance for v in row))
        self.assertTrue(any(v!=0 for row in out.separate.accel_supply.covariance for v in row))

    def test_MAG_and_HOLD_preserve_measurement_supply_counter(self):
        s=state(); n=s.measurement_steps
        word=PRED._preword(s.base.base)
        import test_finite_source_bound_live_word as LBASE
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertEqual(s.measurement_steps,n)
        s,_=X.set_hold(s,hold=False)
        self.assertEqual(s.measurement_steps,n)

    def test_complete_cannot_skip_measurement_supply_recurrence(self):
        with self.assertRaises((ValueError,TypeError)):
            X.complete(state())

    def test_readiness_closes_local_RS_and_accel_propagation_but_not_storage(self):
        r=X.readiness()
        self.assertTrue(r['machine_RS_measurement_effect_attached'])
        self.assertTrue(r['accelerometer_measurement_propagates_machine_prediction_supply'])
        self.assertTrue(r['full_joint24_and_21x21_supply_retained_after_post_S_and_accelerometer'])
        for k in ('machine_Racc_coefficient_displacement_attached',
                  'machine_floor_eigensolver_and_measurement_LDLT_finite_precision_closed',
                  'source_uniform_machine_event_supply_bound_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
