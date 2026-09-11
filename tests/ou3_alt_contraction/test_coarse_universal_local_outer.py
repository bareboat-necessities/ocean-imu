import unittest
from ou3_interval import Interval
from tools.stability.ou3_alt_contraction import coarse_universal_local_outer as C


class CoarseUniversalLocalOuterTest(unittest.TestCase):
    def test_psd_diagonal_ceiling_implies_cross_entry_outer_box(self):
        P=C.psd_covariance_box_from_diagonal_upper([4.0,9.0])
        self.assertLessEqual(P[0][1].lo,-6.0)
        self.assertGreaterEqual(P[0][1].hi,6.0)
        self.assertLessEqual(P[0][0].lo,0.0)
        self.assertGreaterEqual(P[0][0].hi,4.0)

    def test_first_complete_source_local_outer_attempt(self):
        d=C.build(); self.assertEqual(C.validate(d),[])
        # Finiteness is deliberately not asserted.  A failure here identifies
        # the certified projection that must be split; it is not replaced by a
        # trace or source shrink.  Emit the exact failure ledger into CI.
        print('ALT_COARSE_LOCAL_OUTER',C.summary(d))
        self.assertEqual(d['physical_centered_S_bound_m_s'],1100.0)
        self.assertFalse(d['outer_projection_product_is_history_generator'])
        self.assertFalse(d['replay_or_finite_source_sample_used'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':unittest.main()
