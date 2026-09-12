"""Prediction-root provenance regressions; not storage/contraction evidence."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as QX
import test_finite_source_bound_live_word as BASE

PASS=QX.PSDWitness(True)


def operands(state):
    witness,segment,raw,dynamic=BASE.next_imu_operands(state)
    # Existing stationary fixture has zero bias-corrected rate, so no trig
    # witnesses are consumed by the exact small-rate branch.
    dynamic.pop('angular',None); dynamic.pop('ou',None); dynamic.pop('qaxis',None)
    return witness,segment,raw,dynamic


class Tests(unittest.TestCase):
    def test_prediction_roots_are_derived_from_same_source_and_active_state(self):
        s=BASE.root_state(); witness,segment,raw,_=operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        r=X.build(s,physical,raw,ou_alpha=F(199,200),
                  qaxis_marginal_psd=(PASS,PASS,PASS),
                  qaxis_final_psd=(PASS,PASS,PASS),machine_epsilon=F(1,10**7))
        core=s.live.live.live.mekf; active=s.live.live.live.active
        self.assertEqual(r.angular.w,raw.required_bias_corrected_relation(core.z[3:6]))
        self.assertEqual(r.angular.h,segment.h)
        self.assertEqual(r.ou.h,segment.h)
        self.assertEqual(r.ou.tau,active.tau)
        self.assertEqual(r.qaxis.sigma_aw,active.Sigma_aw)
        self.assertFalse(r.qaxis.correlated)

    def test_source_bound_entry_executes_without_free_prediction_roots(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        out=X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='imu-1',
                       ou_alpha=F(199,200),
                       qaxis_marginal_psd=(PASS,PASS,PASS),
                       qaxis_final_psd=(PASS,PASS,PASS),machine_epsilon=F(1,10**7),
                       **dynamic)
        self.assertEqual(out.state.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live.live.live.mekf.reference,segment.after)

    def test_wrong_source_or_ordinal_fails_before_runtime_execution(self):
        s=BASE.root_state(); _,segment,raw,dynamic=operands(s)
        bad=SOURCE.StepWitness(2,'root','c2','p0','p2')
        with self.assertRaisesRegex(ValueError,'next source ordinal'):
            X.imu_step(s,witness=bad,segment=segment,raw=raw,packet_id='bad',
                       ou_alpha=F(199,200),qaxis_marginal_psd=(PASS,PASS,PASS),
                       qaxis_final_psd=(PASS,PASS,PASS),machine_epsilon=F(1,10**7),
                       **dynamic)

    def test_prediction_roots_cannot_be_overridden_by_theorem_caller(self):
        s=BASE.root_state(); witness,segment,raw,dynamic=operands(s)
        dynamic['ou']='detached'
        with self.assertRaisesRegex(TypeError,'cannot be overridden'):
            X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='bad',
                       ou_alpha=F(199,200),qaxis_marginal_psd=(PASS,PASS,PASS),
                       qaxis_final_psd=(PASS,PASS,PASS),machine_epsilon=F(1,10**7),
                       **dynamic)

    def test_readiness_advances_only_prediction_root_provenance(self):
        r=X.readiness()
        self.assertTrue(r['prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation'])
        self.assertTrue(r['caller_cannot_override_angular_OU_Qaxis_roots'])
        self.assertFalse(r['accelerometer_bias_prediction_root_bound_here'])
        self.assertFalse(r['all_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
