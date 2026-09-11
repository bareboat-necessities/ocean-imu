import unittest
from tools.stability.ou3_alt_contraction import bias_families as B

class BiasFamilyContractTests(unittest.TestCase):
    def test_all_three_source_families_are_separate_and_analytic(self):
        cs=B.contracts();self.assertEqual(tuple(c.name for c in cs),('BIAS0','BIAS1','BIAS2'))
        self.assertEqual(len({c.parameter_token for c in cs}),3)
        for c in cs:
            self.assertTrue(c.one_history_required);self.assertGreater(c.driver_norm_bound,0);self.assertGreater(c.true_bias_norm_bound,0)
            self.assertLessEqual(c.phi_true.lo,c.phi_true.hi);self.assertLessEqual(c.phi_true.hi,1)
            Q=B.driver_ball_iqc(c);self.assertEqual((len(Q),len(Q[0])),(4,4));self.assertGreater(Q[0][0].lo,0)
            T=B.true_bias_ball_iqc(c);self.assertGreater(T[0][0].lo,0)
        self.assertTrue(next(c for c in cs if c.name=='BIAS2').phi_true.contains(1.0))
    def test_build_does_not_import_broken_fresh_live_aggregate(self):
        d=B.build();self.assertEqual(B.validate(d),[]);self.assertTrue(d['three_families_invoked_separately']);self.assertTrue(d['one_persistent_parameter_token_per_family']);self.assertFalse(d['aggregate_fresh_Live_entry_builder_consumed']);self.assertFalse(d['all_bias_families_attached_to_complete_word']);self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
