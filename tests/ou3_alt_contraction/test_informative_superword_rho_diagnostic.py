import unittest
from tools.stability.ou3_alt_contraction import informative_superword_rho_diagnostic as D

class TestInformativeSuperwordRhoDiagnostic(unittest.TestCase):
    def test_short_point_service_and_nonpromotion(self):
        d=D.diagnose(24)
        self.assertEqual(D.validate(d),[])
        self.assertFalse(d['source_uniform_rho_certified'])
        self.assertFalse(d['storage_search_allowed'])
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
