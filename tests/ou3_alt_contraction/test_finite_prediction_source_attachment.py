import unittest
from ou3_interval import Interval
from tools.stability.ou3_alt_contraction import finite_prediction_source_attachment as A
from tools.stability.ou3_alt_contraction import bias_families as B


class FinitePredictionSourceAttachmentTests(unittest.TestCase):
    def test_all_mode_family_rows_attach_correlated_sources(self):
        d=A.build();self.assertEqual(A.validate(d),[])
        self.assertEqual(len(d['rows']),6)
        self.assertTrue(d['all_H_A_BIAS_families_have_one_step_finite_source_attachment'])
        self.assertTrue(d['no_independent_q15_or_bias_boxes'])
        self.assertFalse(d['complete_word_finite_identity'])
        self.assertFalse(d['ALT_LIVE_PASS'])
        for r in d['rows']:
            self.assertEqual(r['q15_coordinate_order'],A.Q15_ORDER)
            self.assertEqual(r['q15_dimension'],15)
            self.assertTrue(r['q15_same_history_relation_attached_to_finite_prediction'])
            self.assertTrue(r['same_physical_bias_driver_used_by_error_and_truth'])
            self.assertTrue(r['bias_parameter_root_token_attached_to_finite_prediction'])
            self.assertFalse(r['multi_transition_BRMM_primitive_chain_closed_here'])
            self.assertFalse(r['covariance_frontend_successor_attached'])

    def test_explicit_family_token_is_preserved(self):
        c=B.contracts()[1]
        r=A.attach(mode='A',tau=Interval.point(1.7),h=Interval.point(.005),bias_contract=c)
        self.assertEqual(r['bias_family'],'BIAS1')
        self.assertEqual(r['bias_parameter_token'],c.parameter_token)
        self.assertEqual(r['bias_phi_true_interval'],c.phi_true)

    def test_invalid_mode_or_detached_bias_contract_fails(self):
        c=B.contracts()[0]
        with self.assertRaises(ValueError):A.attach(mode='X',tau=Interval.point(1),h=Interval.point(.005),bias_contract=c)
        with self.assertRaises(TypeError):A.attach(mode='H',tau=Interval.point(1),h=Interval.point(.005),bias_contract=object())


if __name__=='__main__': unittest.main()
