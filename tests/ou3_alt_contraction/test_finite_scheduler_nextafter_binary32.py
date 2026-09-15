from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_scheduler_nextafter_binary32 as X

class Tests(unittest.TestCase):
    def test_predecessor_is_exact_for_power_two_and_ordinary_values(self):
        for x in (B.rn32(F(1)),B.rn32(F(3,2)),B.rn32(F(3,200)),B.rn32(F(1,10))):
            p=X.predecessor_positive(x)
            self.assertTrue(B.is_binary32(p)); self.assertGreater(p,0); self.assertLess(p,x)
        self.assertEqual(X.predecessor_positive(B.rn32(F(1))),F(1)-F(1,1<<24))
        self.assertEqual(X.predecessor_positive(B.rn32(F(3,2))),F(3,2)-F(1,1<<23))

    def test_detached_park_witness_rejected(self):
        x=B.rn32(F(3,2)); p=X.predecessor_positive(x)
        self.assertEqual(X.require_witness(x,p),p)
        with self.assertRaisesRegex(ValueError,'exact binary32 predecessor'):
            X.require_witness(x,x-F(1,100))

    def test_readiness_stays_fail_closed_for_storage(self):
        r=X.readiness()
        self.assertTrue(r['machine_scheduler_nextafter_binary32_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
