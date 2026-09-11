import unittest
from ou3_interval import Interval
from tools.stability.ou3_alt_contraction import h18_a21_edge as E

I=Interval.point

def P18():
    P=[[I(0) for _ in range(18)] for _ in range(18)]
    for i in range(18):P[i][i]=I(1+i/10)
    return P

class H18A21EdgeTests(unittest.TestCase):
    def test_held_covariance_lift_is_enable_floor_fixed_point(self):
        p=E.held_covariance_from_H18(P18())
        self.assertEqual((len(p),len(p[0])),(21,21))
        self.assertTrue(E.covariance_fixed_point(P18()))
        for i in range(3):
            self.assertTrue(p[18+i][18+i].contains(E.DEFAULT_SIGMA_BACC0**2))
            for j in range(18):
                self.assertEqual((p[18+i][j].lo,p[18+i][j].hi),(0,0))
                self.assertEqual((p[j][18+i].lo,p[j][18+i].hi),(0,0))
    def test_joint24_mean_edge_is_identity(self):
        A=E.joint24_mean_edge();self.assertEqual((len(A),len(A[0])),(24,24))
        for i in range(24):
            for j in range(24):self.assertTrue(A[i][j].contains(1 if i==j else 0))
    def test_literal_guard_requires_every_shipping_condition_and_no_hold(self):
        kw=dict(startup_live=True,accel_bias_locked=True,mag_updates_applied=250,mag_updates_to_unlock=250,first_mag_finite=True,elapsed_since_first_mag_s=1.0001,external_hold=False)
        self.assertTrue(E.release_guard(**kw)['learning_enabled_after_edge'])
        for key,value in [('startup_live',False),('accel_bias_locked',False),('mag_updates_applied',249),('first_mag_finite',False),('elapsed_since_first_mag_s',1.0),('external_hold',True)]:
            q=dict(kw);q[key]=value;self.assertFalse(E.release_guard(**q)['learning_enabled_after_edge'])
    def test_source_parity_and_conditional_edge_close_but_reachability_does_not(self):
        d=E.build();self.assertEqual(E.validate(d),[]);self.assertTrue(d['conditional_H18_A21_state_and_covariance_edge_closed']);self.assertTrue(d['release_guard_literal_relation_closed']);self.assertFalse(d['release_guard_source_uniform_reachability_closed']);self.assertFalse(d['H18_A21_complete_word_edge_attached']);self.assertFalse(d['storage_search_allowed'])

if __name__=='__main__':unittest.main()
