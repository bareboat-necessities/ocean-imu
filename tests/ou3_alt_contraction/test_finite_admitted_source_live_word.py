"""Admitted BRMM/BIAS restriction composition regressions; not stability evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as A
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as B
from tools.stability.ou3_alt_contraction import finite_bias_history_restriction as BR
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as X
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as DUAL
import test_finite_source_bound_live_word as BASE
import test_finite_live_magnetic_word as MAGBASE


class Tests(unittest.TestCase):
    def root(self):
        lower=BASE.root_state()
        history=A.AdmittedHistory(lower.source.root.history_id)
        bias=B.AdmittedBiasHistory(lower.bias_history_id,lower.source.root.bias_family,BASE.PHI)
        return X.State(lower,history,bias)

    def restrictions(self,s,witness,segment):
        return (A.RestrictedSegment(s.admitted_history,witness.ordinal,segment),
                B.RestrictedBiasStep(s.bias_history,witness.ordinal,segment))

    def dual_startup(self):
        bridge=MAGBASE.FIRST.startup()
        cfg=MAGBASE.X.Config(gate=MAGBASE.X.GATE.Config(mag_delay=0),
            gravity=MAGBASE.X.GRAVITY.Config(mag_delay=0),
            tuner=MAGBASE.X.TUNER.Config(min_samples=1,min_window=0),
            refinement_start=90,refinement_window=0,continuous_enabled=True)
        model=MAGBASE.X.SOURCE.Model((30,0,0),(0,0,0),'fixed-physical-model')
        proxy=bridge.frontend_before.tuner.vertical
        word=MAGBASE.X.START.State(MAGBASE.X.GRAVITY.State(gravity_good=2,aligned_branch=True),
                                   MAGBASE.PREFIX.State(proxy))
        base=MAGBASE.X.begin_startup(word,cfg,model,bridge.state.mekf.reference.history_id,
                                     'startup-reference')
        cert=DUAL.certify_fresh_startup(base)
        out=DUAL.certified_startup_call(cert,bridge.state.mekf.reference,
            residual_body=(0,0,0),packet_id='startup',
            proxy_q_norm=MAGBASE.X.TILT.SqrtWitness(1,1),proxy_yaw_half=MAGBASE.zero_yaw(),
            hi_decay=MAGBASE.HI.Decay(F(1,200),600,1),
            mag_norm=MAGBASE.X.TUNER.SqrtWitness(900,30),
            mean_norm=MAGBASE.X.TUNER.SqrtWitness(900,30),
            horizontal_sqrt=MAGBASE.GAUGE.HorizontalSqrt(900,30),
            ready_yaw_half=MAGBASE.zero_yaw(30))
        return bridge,out.state

    def test_actual_startup_interleave_binds_exact_admitted_tL_origin(self):
        bridge,magnetic=self.dual_startup()
        ref=bridge.state.mekf.reference
        history=A.AdmittedHistory(ref.history_id)
        origin=A.RestrictedOrigin(history,ref)
        bias=B.AdmittedBiasHistory(ref.bias_root,ref.bias_family,BASE.PHI)
        runtime=BASE.root_state().runtime
        s=X.from_startup(
            bridge,magnetic,origin,bias,
            gyro_residual_history_id='gyro-history',
            accel_residual_history_id='accel-history',runtime=runtime,
            proxy_q_norm=MAGBASE.X.TILT.SqrtWitness(1,1),
            proxy_yaw_half=MAGBASE.zero_yaw())
        self.assertEqual(s.live_word.live.live.live.mekf.reference,ref)
        self.assertEqual(s.admitted_history,history)
        self.assertEqual(s.bias_history,bias)
        self.assertEqual(len(s.live_word.source.steps),0)
        self.assertEqual(magnetic.calls,1)

    def test_old_single_clock_startup_state_cannot_enter_admitted_factory(self):
        bridge,old=MAGBASE.startup()
        ref=bridge.state.mekf.reference
        history=A.AdmittedHistory(ref.history_id); origin=A.RestrictedOrigin(history,ref)
        bias=B.AdmittedBiasHistory(ref.bias_root,ref.bias_family,BASE.PHI)
        with self.assertRaisesRegex(TypeError,'certified dual-clock'):
            X.from_startup(bridge,old,origin,bias,
                gyro_residual_history_id='g',accel_residual_history_id='a',
                runtime=BASE.root_state().runtime,
                proxy_q_norm=MAGBASE.X.TILT.SqrtWitness(1,1),proxy_yaw_half=MAGBASE.zero_yaw())
        # Fresh-startup certification now detects prior ready/reference capture
        # directly; this is the stronger form of the old already-used-word guard.
        with self.assertRaisesRegex(ValueError,'must precede ready/reference capture'):
            DUAL.certify_fresh_startup(old)

    def test_next_IMU_consumes_same_kth_BRMM_and_BIAS_restriction(self):
        s=self.root(); witness,segment,raw,dynamic=BASE.next_imu_operands(s.live_word)
        restricted,bias_restricted=self.restrictions(s,witness,segment)
        out=X.imu_step(s,restricted=restricted,bias_restricted=bias_restricted,
                       witness=witness,raw=raw,packet_id='imu-1',**dynamic)
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.bias_history,s.bias_history)
        self.assertEqual(out.state.live_word.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live_word.live.live.live.mekf.reference,segment.after)

    def test_different_admitted_BRMM_history_cannot_supply_next_segment(self):
        s=self.root(); witness,segment,raw,dynamic=BASE.next_imu_operands(s.live_word)
        other=A.AdmittedHistory('different-primary-history')
        # RestrictedSegment itself now proves history ancestry and therefore
        # rejects this detached history before the shipping composer is called.
        with self.assertRaisesRegex(ValueError,'detached from quantified admitted history'):
            A.RestrictedSegment(other,witness.ordinal,segment)
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_different_admitted_BIAS_history_cannot_supply_next_segment(self):
        s=self.root(); witness,segment,raw,dynamic=BASE.next_imu_operands(s.live_word)
        restricted=A.RestrictedSegment(s.admitted_history,witness.ordinal,segment)
        other=B.AdmittedBiasHistory('different-bias-history',s.bias_history.family,s.bias_history.phi_true)
        bias_restricted=B.RestrictedBiasStep(other,witness.ordinal,segment)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted BIAS history'):
            X.imu_step(s,restricted=restricted,bias_restricted=bias_restricted,
                       witness=witness,raw=raw,packet_id='bad',**dynamic)

    def test_BRMM_BIAS_and_source_ordinals_must_match(self):
        s=self.root(); witness,segment,_,_=BASE.next_imu_operands(s.live_word)
        restricted=A.RestrictedSegment(s.admitted_history,1,segment)
        bad=replace(witness,ordinal=2)
        with self.assertRaisesRegex(ValueError,'ordinal detached'):
            A.qualify_step(s.live_word.source.root,restricted,bad)

    def test_bias_restriction_rejects_factor_detached_from_one_history(self):
        s=self.root(); witness,segment,_,_=BASE.next_imu_operands(s.live_word)
        family=next(r for r in BR.restrictions() if r.name==s.bias_history.family)
        other_phi=F.from_float(family.phi_hi)
        if other_phi == s.bias_history.phi_true:
            other_phi=F.from_float(family.phi_lo)
        self.assertNotEqual(other_phi,s.bias_history.phi_true)
        other=B.AdmittedBiasHistory(s.bias_history.history_id,s.bias_history.family,other_phi)
        with self.assertRaisesRegex(ValueError,'factor detached'):
            B.RestrictedBiasStep(other,witness.ordinal,segment)

    def test_history_is_theorem_quantifier_not_runtime_admission_flag(self):
        h=A.AdmittedHistory('H')
        self.assertEqual(h.canonical_source,A.CANONICAL_SOURCE)
        with self.assertRaises(TypeError):
            A.AdmittedHistory('H',admitted=True)
        r=A.readiness()
        self.assertTrue(r['universal_theorem_quantifier_over_admitted_primary_history_explicit'])
        self.assertTrue(r['finite_word_derived_by_restriction_not_runtime_membership_inference'])
        self.assertTrue(r['arbitrary_runtime_tokens_do_not_prove_COMPLETE_BRMM_membership'])
        self.assertFalse(r['complete_600_step_shipping_word_composed_from_restrictions'])

    def test_readiness_does_not_promote_source_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['universal_admitted_COMPLETE_BRMM_history_carried_in_Live_product'])
        self.assertTrue(r['admitted_BIAS_history_carried_in_same_Live_product'])
        self.assertTrue(r['actual_startup_interleave_constructor_consumed_by_admitted_source_factory'])
        self.assertTrue(r['admitted_startup_requires_inductive_dual_clock_magnetic_history'])
        self.assertTrue(r['already_mutated_single_clock_startup_rejected_at_theorem_handoff'])
        self.assertTrue(r['startup_fresh_H18_reference_must_equal_admitted_tL_origin'])
        self.assertTrue(r['synthetic_Live_root_inserted_by_source_bridge'] is False)
        self.assertTrue(r['every_IMU_event_requires_same_BRMM_and_BIAS_restriction_ordinal'])
        self.assertTrue(r['BIAS_generating_history_attached'])
        self.assertTrue(r['startup_source_equality_bridge_available'])
        self.assertFalse(r['sensor_disturbance_admissibility_attached'])
        self.assertFalse(r['complete_600_step_shipping_word_composed_from_restrictions'])
        self.assertFalse(r['startup_deployment_arithmetic_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
