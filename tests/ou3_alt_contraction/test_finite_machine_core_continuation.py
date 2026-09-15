"""Consecutive local supplies must not impersonate a deployed machine word."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_core_continuation as X
from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as DISP
from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as MROOT
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EROOT
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
import test_finite_source_bound_live_word as WORDTEST
import test_finite_admitted_machine_prediction_supply_interleaved_prefix as BASE
import test_finite_source_bound_prediction_word as ROOTBASE


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.before, cls.first = BASE.executed_supply(pending=True)

    def initial(self):
        return X.begin(self.first.separate.predecessor, self.first.fma.predecessor)

    def first_observation(self):
        return X.observe_imu(self.initial(),
            separate_predecessor=self.first.separate.predecessor,
            fma_predecessor=self.first.fma.predecessor,
            separate_after_accel=self.first.separate.machine,
            fma_after_accel=self.first.fma.machine)

    def test_two_local_predictions_do_not_chain_the_machine_covariance(self):
        # This is a conditional algebra fixture, never source-admission evidence.
        # All subsequent measurement branches in the lower fixture reject, so
        # its stored CORE is exactly the shadow prediction.
        word = BASE.X._preword(self.first.state.base)
        core = word.live.live.live.mekf
        self.assertEqual(core, self.first.separate.exact)
        self.assertNotEqual(core.covariance, self.first.separate.machine.covariance)
        residual = max(abs(x-y) for a,b in zip(core.covariance,
                       self.first.separate.machine.covariance) for x,y in zip(a,b))
        self.assertGreater(residual, F(14, 1000))
        self.assertLess(residual, F(15, 1000))

        witness, segment, raw, _ = ROOTBASE.operands(word)
        physical = SOURCE.QualifiedPhysicalSegment(word.source.root, witness, segment)
        args = ROOTBASE.root_args(); args.pop('temperature_c')
        exact = EROOT.build(word, physical, raw, **args)
        relations = []
        for mode in ('separate', 'fma'):
            old = getattr(self.first, mode)
            machine_args = dict(args, ou_alpha=old.machine_roots.ou.alpha,
                                ou_em1=old.machine_roots.ou.em1)
            machine = MROOT.build(word, physical, raw, old.machine_roots.active,
                mode=mode, exact_active=word.live.live.live.active, **machine_args)
            relations.append(DISP.compare(core, segment, raw, exact, machine,
                                           Qbase=word.runtime.Qbase))
        carried = X.observe_imu(self.first_observation(),
            separate_predecessor=relations[0].predecessor,
            fma_predecessor=relations[1].predecessor,
            separate_after_accel=relations[0].machine,
            fma_after_accel=relations[1].machine)
        self.assertFalse(carried.observations[-1].separate_predecessor_connected)
        self.assertFalse(carried.observations[-1].fma_predecessor_connected)
        self.assertIn('separate_machine_CORE_successor_not_consumed', carried.missing)
        self.assertEqual(carried.separate, relations[0].machine)

    def test_MAG_HOLD_do_not_silently_use_the_shadow_successor(self):
        initial = self.first_observation()
        state = X.observe_async(initial, kind='mag')
        state = X.observe_async(state, kind='hold')
        self.assertIs(state.separate, initial.separate)
        self.assertEqual(state.imu_steps, 1)
        self.assertIn('mag_machine_CORE_and_control_successor_not_attached', state.missing)
        self.assertIn('hold_machine_CORE_and_control_successor_not_attached', state.missing)

    def test_600_counters_cannot_replace_event_correspondence(self):
        events = tuple(X.Observation('imu', n, True, True) for n in range(1, 601))
        forged = replace(self.initial(), observations=events)
        with self.assertRaisesRegex(ValueError, 'successor algorithms not attached'):
            X.require_complete(forged)

    def test_covariance_atlas_bias_and_source_are_part_of_predecessor_identity(self):
        state = self.first_observation()
        self.assertEqual(len(state.separate.z), 24)
        self.assertEqual(len(state.separate.covariance), 21)
        self.assertEqual(len(state.separate.covariance[0]), 21)
        self.assertEqual(state.separate.reference.beta, state.separate.z[21:])
        self.assertIn('machine_post_accelerometer_event_and_control_suffix_not_attached',
                      state.missing)

    def test_strongest_word_requires_this_continuation(self):
        r = LIVE.readiness()
        self.assertTrue(r['complete_word_requires_machine_CORE_and_control_continuation'])
        self.assertTrue(r['machine_CORE_local_successors_persist_and_next_predecessors_checked'])
        self.assertTrue(r['machine_CORE_and_control_successor_algorithms_attached'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])

    def test_distinct_predecessor_relation_keeps_prior_covariance_discrepancy(self):
        word=BASE.X._preword(self.first.state.base)
        witness,segment,raw,_=ROOTBASE.operands(word)
        physical=SOURCE.QualifiedPhysicalSegment(word.source.root,witness,segment)
        args=ROOTBASE.root_args(); args.pop('temperature_c')
        exact=EROOT.build(word,physical,raw,**args)
        previous=self.first.separate.machine
        old=self.first.separate.machine_roots
        machine=MROOT.build(word,physical,raw,old.active,mode='separate',
            machine_predecessor=previous,exact_active=word.live.live.live.active,
            **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))
        out=DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
            machine_predecessor=previous,Qbase=word.runtime.Qbase)
        direct=DISP.PRED.prediction_from_raw(previous,segment,raw,
            angular=machine.angular,ou=machine.ou,bias=machine.bias,qaxis=machine.qaxis,
            Qbase=word.runtime.Qbase)
        self.assertEqual(out.machine,direct)
        self.assertIs(out.machine_predecessor,previous)
        self.assertNotEqual(out.machine_predecessor.covariance,out.predecessor.covariance)
        with self.assertRaisesRegex(ValueError,'roots detached from carried predecessor'):
            DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
                         Qbase=word.runtime.Qbase)

    def _actual_runtime(self,*,unlocked=False):
        runtime=BASE.X._preword(self.first.state.base).live
        core=self.first.separate.machine
        magnetic=runtime.magnetic
        if unlocked:
            core=GATE._active(core,magnetic.memory.cfg.gate.initial_ba_std)
            cov=[list(row) for row in core.covariance]
            cov[6][18]=cov[18][6]=F(1,100)
            core=replace(core,covariance=tuple(tuple(row) for row in cov))
            magnetic=replace(magnetic,control=replace(magnetic.control,locked=False,hold=False))
        runtime=replace(runtime,live=replace(runtime.live,
            live=replace(runtime.live.live,mekf=core)),magnetic=magnetic)
        return X.begin_from_interleaved(runtime)

    def test_carried_gyro_bias_error_drives_its_own_angular_roots_and_attitude(self):
        word=BASE.X._preword(self.first.state.base)
        witness,segment,raw,_=ROOTBASE.operands(word)
        physical=SOURCE.QualifiedPhysicalSegment(word.source.root,witness,segment)
        args=ROOTBASE.root_args(); args.pop('temperature_c')
        exact=EROOT.build(word,physical,raw,**args)
        # A nonzero bias error inside the existing small-rate branch needs no
        # caller-selected trigonometric witness, and still changes attitude.
        z=list(self.first.separate.machine.z); z[3]+=F(1,10**8)
        previous=replace(self.first.separate.machine,z=tuple(z))
        old=self.first.separate.machine_roots
        machine=MROOT.build(word,physical,raw,old.active,mode='separate',
            machine_predecessor=previous,exact_active=word.live.live.live.active,
            **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))
        self.assertNotEqual(machine.angular.w,exact.angular.w)
        self.assertEqual(machine.angular.w,raw.required_bias_corrected_relation(previous.z[3:6]))
        out=DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
            machine_predecessor=previous,Qbase=word.runtime.Qbase)
        self.assertNotEqual(out.machine.q_hat,out.exact.q_hat)
        self.assertEqual(out.machine.reference,out.exact.reference)

    def test_concrete_HOLD_release_preserves_machine_mean_and_full_covariance(self):
        state=self._actual_runtime(unlocked=True)
        core=state.separate
        held,events=X.advance_hold(state,hold=True)
        self.assertEqual(held.separate.mode,'H')
        self.assertEqual(held.separate.covariance[6][18],0)
        self.assertEqual(held.separate.z,core.z)
        self.assertEqual(held.separate.q_hat,core.q_hat)
        self.assertEqual(held.separate.attitude_chart,core.attitude_chart)
        released,events=X.advance_hold(held,hold=False)
        self.assertEqual(released.separate.mode,'A')
        self.assertEqual(released.separate.z,core.z)
        self.assertEqual(released.separate_runtime.live.live.mekf,released.separate)
        self.assertEqual(released.missing,())
        self.assertEqual(released.imu_steps,0)

    def test_concrete_MAG_consumes_persisted_machine_covariance_and_control(self):
        state=self._actual_runtime()
        word=BASE.X._preword(self.first.state.base)
        kwargs=WORDTEST.mag_kwargs(word)
        out,events=X.advance_mag(state,**kwargs)
        expected=X.LIVE.mag_step(state.separate_runtime,**kwargs)
        self.assertEqual(out.separate,expected.state.live.live.mekf)
        self.assertEqual(out.separate_runtime.magnetic,expected.state.magnetic)
        self.assertEqual(out.separate.reference,state.separate.reference)
        self.assertEqual(out.imu_steps,0)
        self.assertEqual(out.missing,())
        with self.assertRaisesRegex(TypeError,'override common source operands'):
            X.advance_mag(state,separate={'packet_id':'spliced'},**kwargs)


if __name__ == '__main__':
    unittest.main()
