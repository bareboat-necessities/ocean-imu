import unittest
from tools.stability.ou3_alt_contraction import inverse_free_universal_local_outer as L

class InverseFreeUniversalLocalOuterTest(unittest.TestCase):
    def test_corrected_local_family_does_not_use_endpoint_covariance(self):
        d=L.build();self.assertEqual(L.validate(d),[])
        print('ALT_INVERSE_FREE_LOCAL_OUTER',L.summary(d))
        self.assertTrue(d['inverse_free_gain_outer_consumed'])
        self.assertFalse(d['endpoint_Pbar_used_for_local_gain'])
        self.assertFalse(d['independent_covariance_box_used'])
        self.assertFalse(d['magnetometer_count_assumption_used'])
        self.assertEqual(d['physical_D_S_bound_m_s'],1100.0)
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
