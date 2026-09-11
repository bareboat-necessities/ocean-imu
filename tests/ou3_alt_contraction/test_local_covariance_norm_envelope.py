import unittest
from ou3_interval import matrix_point
from tools.stability.ou3_alt_contraction import local_covariance_norm_envelope as C


class LocalCovarianceNormEnvelopeTest(unittest.TestCase):
    def test_interval_norm_upper_contains_point_spectral_norm(self):
        A=matrix_point([[3.0,0.0],[4.0,0.0]])
        # exact spectral norm is 5; sqrt(||A||1 ||A||inf)=sqrt(7*4)>5.
        self.assertGreaterEqual(C.interval_matrix_norm2_upper(A),5.0)

    def test_event_local_covariance_bound_uses_endpoint_Pbar_only_at_word_boundaries(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        print('ALT_LOCAL_COVARIANCE_NORM',d['reports'])
        self.assertTrue(d['event_local_covariance_norm_envelope_closed'])
        self.assertTrue(d['aw_floor_bound_independent_of_current_P'])
        self.assertTrue(d['endpoint_Pbar_used_only_at_certified_sliding_word_endpoints'])
        self.assertFalse(d['endpoint_Pbar_substituted_directly_for_preJoseph_covariance'])
        self.assertTrue(d['H_to_A_release_covariance_ancestry_retained'])
        self.assertTrue(d['A_release_after_H_endpoint_regime_guaranteed'])
        self.assertFalse(d['trajectory_replay_used'])

        H=d['reports']['H']; A=d['reports']['A']
        self.assertGreaterEqual(H['early_word_every_event_upper'], H['live_seed_norm_upper'])
        self.assertGreaterEqual(A['early_active_word_every_event_upper'], A['release_seed_from_H_endpoint_plus_fixed_ba_upper'])
        self.assertGreaterEqual(A['release_seed_from_H_endpoint_plus_fixed_ba_upper'], H['sliding_endpoint_Pbar_norm_upper'])


if __name__=='__main__':unittest.main()
