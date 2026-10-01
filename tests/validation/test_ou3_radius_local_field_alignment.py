"""Exact/rational checks for the radius-local field-alignment exclusion."""
from fractions import Fraction as F
import unittest

G = F(196133, 20000)  # 9.80665
SIGMA_W = F(1, 5)
VMAX = F(11, 2)
JMAX = F(100)
HACC = F(3, 500)
KAW = F(203, 50)  # 4.06 outward-safe source-audited coefficient

def m_phys(T):
    T=F(T)
    return G*SIGMA_W - 2*VMAX/T - JMAX*HACC

def margin(r,T):
    return m_phys(T)-KAW*F(r)

class RadiusLocalFieldAlignment(unittest.TestCase):
    def test_authoritative_seventeen_second_separation(self):
        self.assertEqual(m_phys(17), F(1214261,1700000))
        self.assertGreater(m_phys(17), F('0.71427'))

    def test_positive_separation_threshold(self):
        # m_phys(T)>0 iff T > 11/(g/5-.6)
        tcrit=F(11)/(G*SIGMA_W-JMAX*HACC)
        self.assertGreater(tcrit,F(8))
        self.assertLess(tcrit,F('8.1'))

    def test_outward_safe_radius_max(self):
        rmax=m_phys(17)/KAW
        self.assertGreater(rmax,F('0.1759'))
        self.assertLess(rmax,F('0.1760'))

    def test_conservative_radius_closes(self):
        self.assertEqual(margin(F(3,20),17),F(178961,1700000))
        self.assertGreater(margin(F(3,20),17),F('0.10527'))

    def test_ideal_four_r_is_not_used_for_promotion(self):
        self.assertGreater(KAW,F(4))

if __name__=='__main__': unittest.main()
