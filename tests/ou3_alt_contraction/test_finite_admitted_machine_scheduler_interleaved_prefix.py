"""Persistent per-compiler machine S-scheduler regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_scheduler_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_interleaved_prefix as PRED
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_scheduler_nextafter_binary32 as NEXT
import test_finite_admitted_machine_tunestate_interleaved_prefix as BASE
import test_finite_source_bound_live_word as LBASE


def state():
    p=PRED.begin(BASE.state())
    sa=p.base.separate_active; fa=p.base.fma_active
    ss=POST.Scheduler(sa.pseudo_period,F(0),F(0))
    fs=POST.Scheduler(fa.pseudo_period,F(0),F(0))
    return X.begin(p,ss,fs)


class Tests(unittest.TestCase):
    def test_begin_binds_scheduler_periods_to_applied_machine_states(self):
        s=state()
        self.assertEqual(s.separate_scheduler.period,s.base.base.separate_active.pseudo_period)
        self.assertEqual(s.fma_scheduler.period,s.base.base.fma_active.pseudo_period)
        self.assertEqual(s.scheduler_steps,0)

    def test_state_rejects_scheduler_period_detached_from_machine_active(self):
        s=state(); bad=POST.Scheduler(s.separate_scheduler.period+F(1,100),0)
        with self.assertRaisesRegex(ValueError,'scheduler period detached'):
            X.State(s.base,bad,s.fma_scheduler,s.entry_imu_steps,s.scheduler_steps)

    def test_no_boundary_advance_uses_same_scheduler_without_retarget(self):
        s=state(); active=s.base.base.separate_active
        out=X._advance_one(s.separate_scheduler,active,boundary_consumed=False,h=F(1,200),park=None)
        self.assertFalse(out.retargeted)
        self.assertEqual(out.after_retarget,s.separate_scheduler)
        self.assertEqual(out.after_step.elapsed,F(1,200))

    def test_boundary_retarget_preserves_credit_when_below_new_period(self):
        s=state(); before=POST.Scheduler(s.separate_scheduler.period,F(1,100))
        active=ACTIVE.ActiveParameters(s.base.base.separate_active.tau,
             s.base.base.separate_active.Sigma_aw,s.separate_scheduler.period+F(1,10),
             s.base.base.separate_active.R_S)
        out=X._advance_one(before,active,boundary_consumed=True,h=F(1,200),park=None)
        self.assertTrue(out.retargeted)
        self.assertEqual(out.after_retarget.elapsed,before.elapsed)
        self.assertEqual(out.after_retarget.period,active.pseudo_period)
        self.assertEqual(out.after_step.elapsed,F(3,200))

    def test_overdue_retarget_requires_exact_nextafter_witness(self):
        s=state(); old=s.separate_scheduler.period
        before=POST.Scheduler(old,old-F(1,1000))
        newp=B.rn32(old/F(2))
        active=ACTIVE.ActiveParameters(s.base.base.separate_active.tau,
             s.base.base.separate_active.Sigma_aw,newp,s.base.base.separate_active.R_S)
        with self.assertRaisesRegex(ValueError,'nextafter witness'):
            X._advance_one(before,active,boundary_consumed=True,h=F(1,200),park=None)
        park=ACTIVE.NextafterParkWitness(NEXT.predecessor_positive(newp))
        out=X._advance_one(before,active,boundary_consumed=True,h=F(1,200),park=park)
        self.assertEqual(out.after_retarget.elapsed,park.parked_elapsed)
        with self.assertRaisesRegex(ValueError,'exact binary32 predecessor'):
            X._advance_one(before,active,boundary_consumed=True,h=F(1,200),
                           park=ACTIVE.NextafterParkWitness(newp-F(1,10000)))

    def test_MAG_and_HOLD_preserve_both_machine_schedulers(self):
        s=state(); ss=s.separate_scheduler; fs=s.fma_scheduler
        word=PRED._word(s.base.base)
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertIs(s.separate_scheduler,ss); self.assertIs(s.fma_scheduler,fs)
        s,_=X.set_hold(s,hold=False)
        self.assertIs(s.separate_scheduler,ss); self.assertIs(s.fma_scheduler,fs)

    def test_complete_cannot_skip_scheduler_recurrence(self):
        with self.assertRaises((ValueError,TypeError)):
            X.complete(state())

    def test_readiness_closes_scheduler_effect_only(self):
        r=X.readiness()
        self.assertTrue(r['machine_pseudo_period_scheduler_effect_attached'])
        self.assertTrue(r['machine_due_not_due_branch_retained_per_compiler_history'])
        self.assertTrue(r['scheduler_nextafter_binary32_correspondence_closed'])
        self.assertFalse(r['machine_RS_measurement_effect_attached'])
        self.assertFalse(r['machine_due_branch_S_measurement_composed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
