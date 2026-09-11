import unittest
from ou3_interval import matrix_point
from tools.stability.ou3_alt_contraction import inverse_free_gain_outer as G

class InverseFreeGainOuterTest(unittest.TestCase):
    def test_gershgorin_R_lower(self):
        R=matrix_point([[2.0,0.1],[0.1,3.0]])
        self.assertGreater(G.gershgorin_spd_lower(R),1.8)
    def test_source_uniform_gain_bounds_are_finite_without_inverse(self):
        d=G.build();self.assertEqual(G.validate(d),[])
        print('ALT_INVERSE_FREE_GAIN',{m:{k:r['K_norm_upper'] for k,r in rows.items()} for m,rows in d['reports'].items()})
        self.assertTrue(d['all_event_gain_outers_finite'])
        self.assertFalse(d['innovation_inverse_computed'])
        self.assertFalse(d['endpoint_Pbar_used_for_local_gain'])
        self.assertFalse(d['independent_covariance_box_used'])

if __name__=='__main__':unittest.main()
