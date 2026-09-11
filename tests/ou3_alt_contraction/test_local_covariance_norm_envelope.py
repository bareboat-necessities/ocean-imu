import unittest
from ou3_interval import matrix_point
from tools.stability.ou3_alt_contraction import local_covariance_norm_envelope as C


class LocalCovarianceNormEnvelopeTest(unittest.TestCase):
    def test_interval_norm_upper_contains_point_spectral_norm(self):
        A=matrix_point([[3.0,0.0],[4.0,0.0]])
        # exact spectral norm is 5; sqrt(||A||1 ||A||inf)=sqrt(7*4)>5.
        self.assertGreaterEqual(C.interval_matrix_norm2_upper(A),5.0)

    def test_event_local_covariance_bound_does_not_use_endpoint_Pbar(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        print('ALT_LOCAL_COVARIANCE_NORM',d['reports'])
        self.assertTrue(d['event_local_covariance_norm_envelope_closed'])
        self.assertTrue(d['aw_floor_bound_independent_of_current_P'])
        self.assertFalse(d['endpoint_referenced_Pbar_used_as_local_bound'])
        self.assertFalse(d['trajectory_replay_used'])


if __name__=='__main__':unittest.main()
