import unittest
from tools.stability.ou3_alt_contraction import informative_superword_rho_diagnostic as D

class TestInformativeSuperwordRhoDiagnostic(unittest.TestCase):
    def test_short_point_service_and_nonpromotion(self):
        d=D.diagnose(24)
        self.assertEqual(D.validate(d),[])
        self.assertFalse(d['source_uniform_rho_certified'])
        self.assertFalse(d['storage_search_allowed'])
        self.assertFalse(d['ALT_LIVE_PASS'])
        self.assertTrue(d['identity_metric_rejected_by_this_experiment'])
        self.assertFalse(d['corrected_formulation_falsified_by_this_experiment'])
        d['corrected_formulation_falsified_by_this_experiment'] = True
        self.assertIn('corrected_formulation_falsified_by_this_experiment not false', D.validate(d))

if __name__=='__main__': unittest.main()
