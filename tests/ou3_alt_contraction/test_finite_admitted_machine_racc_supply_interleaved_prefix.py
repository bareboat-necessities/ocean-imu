"""Same-event raw applied sigma -> machine Racc -> accelerometer regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_racc_supply_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as M
from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as PRED
from tools.stability.ou3_alt_contraction import finite_admitted_machine_scheduler_interleaved_prefix as SCHED
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_machine_measurement_supply_interleaved_prefix as MBASE
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE
import test_finite_live_magnetic_word as MAG
import test_finite_source_bound_live_word as LBASE
from test_finite_machine_frontend_sigma_source import ready_pair


def state(): return X.begin(MBASE.state())


def event_operands(s):
    # Racc configuration and nominal std are persistent runtime state and must
    # not be supplied as event-local operands.
    return dict(MBASE.event_operands(s.base),
                separate_accel_ldlt=MAG.REJECT,fma_accel_ldlt=MAG.REJECT)


def pending_measurement_state():
    """Build the same admitted wrapper stack from an existing pending MTUNE fixture."""
    mt=TBASE.state(usable=True,pending=True)
    mt=replace(mt,frontends=ready_pair())
    prefix=mt.base.base.prefix; word=prefix.prefix.live.live_word
    bench=B.rn32(F(1,10))
    word=replace(word,runtime=replace(word.runtime,boundary_bench_noise_sigma=bench,bench_noise_sigma=bench))
    admitted=replace(prefix.prefix.live,live_word=word)
    prefix=replace(prefix,prefix=replace(prefix.prefix,live=admitted))
    mt=replace(mt,base=replace(mt.base,base=replace(mt.base.base,prefix=prefix)))
    p=SCHED.LOWER.begin(mt)
    sa=p.base.separate_active; fa=p.base.fma_active
    ss=POST.Scheduler(sa.pseudo_period,F(0),F(0)); fs=POST.Scheduler(fa.pseudo_period,F(0),F(0))
    sched=SCHED.begin(p,ss,fs)
    return M.begin(PRED.begin(sched))


class Tests(unittest.TestCase):
    def test_begin_carries_raw_sigma_separately_from_stationary_covariance(self):
        s=state(); exact=X._exact_live_state(s.base)
        self.assertEqual(s.separate_applied_sigma,exact.tuner.tune.sigma_applied)
        self.assertEqual(s.fma_applied_sigma,exact.tuner.tune.sigma_applied)
        altered=replace(s,separate_applied_sigma=F(7,100))
        self.assertEqual(altered.separate_applied_sigma,F(7,100))
        self.assertEqual(altered.base,s.base)

    def test_nonpending_event_preserves_raw_sigma_and_reexecutes_Racc(self):
        s=state(); kw=event_operands(s)
        out=X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,
            fma_racc_accel_ldlt=MAG.REJECT,**kw)
        live=M._live_result(out.lower.lower)
        self.assertEqual(out.separate.applied_sigma,s.separate_applied_sigma)
        self.assertEqual(out.fma.applied_sigma,s.fma_applied_sigma)
        self.assertEqual(out.state.racc_steps,1)
        self.assertEqual(out.separate.racc.excess_rms,live.guarded.guard.excess_rms)
        self.assertEqual(out.fma.racc.excess_rms,live.guarded.guard.excess_rms)
        self.assertFalse(out.separate.accelerometer.accepted)

    def test_pending_boundary_source_of_raw_sigma_is_mode_stored_snapshot(self):
        ms=pending_measurement_state(); rs=X.begin(ms)
        # This test isolates the boundary-to-Racc state edge.  Prediction roots
        # are a different obligation and, after a changed pending tau, require
        # their own post-boundary exp witnesses.  Execute the same MTUNE event
        # directly so the raw-sigma source cannot be obscured by that later edge.
        mt_state=ms.base.base.base.base
        kwargs,machine=TBASE.event_operands(mt_state)
        mt=TBASE.X.imu_step(mt_state,**machine,**kwargs,
            separate_boundary_noise_sqrt_gain=11,
            fma_boundary_noise_sqrt_gain=12)
        self.assertTrue(mt.machine_boundary.consumed)
        sep,fma=X._applied_sigmas(rs,mt)
        self.assertEqual(sep,mt.machine_boundary.arithmetic.separate.stored_sigma)
        self.assertEqual(fma,mt.machine_boundary.arithmetic.fma.stored_sigma)

    def test_machine_Racc_state_is_persistent_and_MAG_HOLD_preserve_it(self):
        s=state(); sr=s.separate_racc; fr=s.fma_racc; ss=s.separate_applied_sigma
        word=M._preword(s.base)
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertEqual(s.separate_racc,sr); self.assertEqual(s.fma_racc,fr)
        self.assertEqual(s.separate_applied_sigma,ss)
        s,_=X.set_hold(s,hold=False)
        self.assertEqual(s.separate_racc,sr); self.assertEqual(s.fma_racc,fr)
        self.assertEqual(s.separate_applied_sigma,ss)

    def test_wrapper_requires_distinct_machine_Racc_accelerometer_LDLT(self):
        s=state(); kw=event_operands(s)
        with self.assertRaisesRegex(TypeError,'both machine-Racc accelerometer LDLT branches required'):
            X.imu_step(s,separate_racc_accel_ldlt=None,fma_racc_accel_ldlt=MAG.REJECT,**kw)

    def test_explicit_applied_sigma_adapter_is_part_of_Racc_proof_surface(self):
        self.assertTrue(RACC.readiness()['explicit_applied_sigma_adapter_available_for_machine_history'])

    def test_complete_cannot_skip_Racc_recurrence(self):
        with self.assertRaises((ValueError,TypeError)):
            X.complete(state())

    def test_readiness_closes_machine_Racc_coefficient_path_but_not_native_or_storage(self):
        r=X.readiness()
        for k in ('persistent_Racc_config_and_nominal_std_consumed_from_admitted_runtime',
                  'raw_applied_sigma_carried_separately_from_stationary_Sigma_aw',
                  'pending_machine_boundary_updates_same_mode_applied_sigma_before_Racc',
                  'persistent_separate_and_FMA_Racc_states_attached',
                  'same_mode_preupdate_WPE_frequency_drives_machine_Racc',
                  'machine_Racc_coefficient_displacement_attached',
                  'machine_Racc_effect_propagated_through_accelerometer'):
            self.assertTrue(r[k])
        for k in ('machine_Racc_binary32_hypot_sqrt_correspondence_closed',
                  'source_uniform_machine_event_supply_bound_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
