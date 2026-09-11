import unittest
from tools.stability.ou3_alt_contraction import inverse_free_common_metric_attempt as M

class InverseFreeEndpointAndMetricAttemptTest(unittest.TestCase):
    def test_corrected_600_step_endpoint_attempt(self):
        from tools.stability.ou3_alt_contraction import inverse_free_endpoint_outer_attempt as E
        d=E.build();self.assertEqual(E.validate(d),[])
        print('ALT_INVERSE_FREE_ENDPOINT',E.summary(d))
        self.assertTrue(d['inverse_free_local_family_consumed'])
        self.assertFalse(d['endpoint_Pbar_route_consumed'])
        self.assertFalse(d['common_storage_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])
    def test_corrected_first_common_metric_attempt(self):
        with self.assertRaisesRegex(RuntimeError,'finite-state storage blocked'):
            M.build()

if __name__=='__main__':unittest.main()
