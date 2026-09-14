"""Same-event guard/Mahony/frontend/sigma -> Racc measurement join regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_scheduler_interleaved_prefix as SCHED
from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as PRED
from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as M
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_machine_racc_supply_interleaved_prefix as RBASE
import test_finite_admitted_machine_vertical_stillness_interleaved_prefix as VBASE
import test_finite_source_bound_prediction_word as PBASE
import test_finite_live_magnetic_word as MAG
import test_finite_source_bound_live_word as LBASE


def _measurement_state_from_mtune(mt):
    """Build the existing prediction/scheduler/measurement wrappers around mt."""
    p=SCHED.LOWER.begin(mt)
    sa=p.base.separate_active; fa=p.base.fma_active
    ss=POST.Scheduler(sa.pseudo_period,F(0),F(0))
    fs=POST.Scheduler(fa.pseudo_period,F(0),F(0))
    sched=SCHED.begin(p,ss,fs)
    return M.begin(PRED.begin(sched))


def fixture():
    # First construct the stronger source fixture.  Its machine frontend uses
    # the guard/private-Mahony vertical output rather than the old x=0 component
    # placeholder.  Then wrap THAT SAME MTUNE predecessor in the existing
    # prediction/S/Racc layers; no physical event has yet been executed.
    vertical,tkw,machine,join,guarded,sep,fma=VBASE.executed_fixture()
    mt=vertical.base
    r=RBASE.X.begin(_measurement_state_from_mtune(mt))

    # RBASE supplies prediction/S/accelerometer witnesses. Replace only its
    # old conditional frontend/sigma witnesses with the source-owned witnesses
    # from VBASE so the nested TuneState event and the joined source relation
    # consume exactly the same machine band input/stillness history.
    kw=RBASE.event_operands(r)
    kw.update(machine)
    s=X.begin(r,guard=vertical.guard,guard_cfg=vertical.guard_cfg,
              separate_source=vertical.separate_source,fma_source=vertical.fma_source)
    return s,kw,join,guarded,sep,fma


class Tests(unittest.TestCase):
    def test_same_Racc_event_consumes_frontend_and_sigma_from_joined_machine_source(self):
        s,kw,join,guarded,sep,fma=fixture()
        out=X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        mt=X.LOWER._mtune_result(out.lower.lower)
        self.assertEqual(out.guard,guarded)
        self.assertEqual(out.separate_source,sep); self.assertEqual(out.fma_source,fma)
        self.assertEqual(mt.separate_frontend.band.envelope.x,out.separate_source.band_input)
        self.assertEqual(mt.fma_frontend.band.envelope.x,out.fma_source.band_input)
        self.assertEqual(mt.separate_sigma_join.machine.still_time,out.separate_source.stillness.state.still_time)
        self.assertEqual(out.state.source_steps,1)
        self.assertEqual(out.state.base.racc_steps,s.base.racc_steps+1)

    def test_frontend_input_splice_is_rejected_against_nested_Racc_TuneState_event(self):
        s,kw,join,_,_,_=fixture()
        # The lower machine frontend is fixed by its source-owned witness;
        # changing only the joined Mahony API source must fail at the same-event
        # frontend equality rather than creating a second accepted history.
        bad=dict(join); a=list(bad['machine_acc_body']); a[2]=B.rn32(F(a[2])+F(1,100)); bad['machine_acc_body']=tuple(a)
        with self.assertRaises((ValueError,TypeError)):
            X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**bad,**kw)

    def test_detached_guard_runtime_config_is_rejected_before_join(self):
        s,_,_,_,_,_=fixture()
        with self.assertRaisesRegex(ValueError,'guard configuration detached'):
            replace(s,guard_cfg=replace(s.guard_cfg,cutoff_hz=B.rn32(13)))

    def test_MAG_and_HOLD_preserve_joined_machine_source_history(self):
        s,_,_,_,_,_=fixture(); g=s.guard; a=s.separate_source; n=s.source_steps
        pre=X._mtune_state(s.base).base.base.prefix.prefix.live.live_word
        s,_=X.mag_step(s,**LBASE.mag_kwargs(pre)); self.assertIs(s.guard,g); self.assertIs(s.separate_source,a); self.assertEqual(s.source_steps,n)
        s,_=X.set_hold(s,hold=False); self.assertIs(s.guard,g); self.assertIs(s.separate_source,a); self.assertEqual(s.source_steps,n)

    def test_complete_cannot_skip_join_on_600_measurement_edges(self):
        s,_,_,_,_,_=fixture()
        with self.assertRaises((ValueError,TypeError)):
            X.complete(s)

    def test_readiness_closes_parallel_product_gap_only(self):
        r=X.readiness()
        for k in ('same_executed_TuneState_event_supplies_frontend_and_sigma_join',
                  'no_second_physical_or_filter_event_executed_for_frontend_ancestry',
                  'common_machine_guard_and_private_Mahony_history_joined_to_Racc_word',
                  'machine_band_input_bound_to_frontend_consumed_by_Racc_word',
                  'machine_sigma_stillness_bound_to_sigma_target_consumed_by_Racc_word',
                  'machine_guard_runtime_config_bound_in_joined_word',
                  'tracker_LPF_default_cutoff_bound_in_joined_word',
                  'complete_word_requires_join_on_all_600_Racc_measurement_IMU_edges'):
            self.assertTrue(r[k])
        for k in ('tracker_LPF_mutable_setter_ancestry_closed','private_Mahony_mutable_config_setter_ancestry_closed',
                  'startup_joined_machine_history_attached','target_libm_Eigen_and_compiler_profile_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
