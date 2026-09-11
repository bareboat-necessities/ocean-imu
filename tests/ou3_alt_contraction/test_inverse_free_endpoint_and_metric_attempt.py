import unittest
from tools.stability.ou3_alt_contraction import inverse_free_endpoint_outer_attempt as E
from tools.stability.ou3_alt_contraction import inverse_free_common_metric_attempt as M

class InverseFreeEndpointAndMetricAttemptTest(unittest.TestCase):
    def test_corrected_600_step_endpoint_attempt(self):
        d=E.build();self.assertEqual(E.validate(d),[])
        print('ALT_INVERSE_FREE_ENDPOINT',E.summary(d))
        self.assertTrue(d['inverse_free_local_family_consumed'])
        self.assertFalse(d['endpoint_Pbar_route_consumed'])
        self.assertFalse(d['common_storage_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])
    def test_corrected_first_common_metric_attempt(self):
        d=M.build();self.assertEqual(M.validate(d),[])
        print('ALT_INVERSE_FREE_COMMON_METRIC',M.summary(d))
        self.assertTrue(d['corrected_inverse_free_endpoint_outer_consumed'])
        self.assertFalse(d['superseded_endpoint_Pbar_local_route_consumed'])
        self.assertTrue(d['same_M_all_modes_bias_families'])
        self.assertFalse(d['async_magnetometer_nonexpansive_same_M_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
