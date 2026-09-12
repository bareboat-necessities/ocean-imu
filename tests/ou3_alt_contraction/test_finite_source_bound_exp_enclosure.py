"""Same-source exp/expm1 real-enclosure regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as X
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as OU


class Tests(unittest.TestCase):
    def test_ou_fixture_endpoints_are_rigorous_for_x_005(self):
        self.assertEqual(X.enclosure(F(1,200)),
                         (F(199,200),F(79601,80000),F(-1,200),F(-399,80000)))
        d=OU.OUDecay(F(1,200),1,F(199,200),em1=F(-1,200))
        self.assertIs(X.validate_ou(d),d)

    def test_detached_exp_and_expm1_are_rejected_independently(self):
        with self.assertRaisesRegex(ValueError,'exp root detached'):
            X.validate_exp_expm1(F(1,200),F(9,10),F(-1,200))
        with self.assertRaisesRegex(ValueError,'expm1 root detached'):
            X.validate_exp_expm1(F(1,200),F(199,200),F(-1,10))

    def test_distinct_active_BA_arguments_are_checked(self):
        h=F(1,200)
        d=OU.BiasDecay(True,5000,F(999999,1000000),
                       ((0,0,0),(0,0,0),(0,0,0)),em1_2=F(-1,500000))
        self.assertIs(X.validate_bias(d,h),d)
        bad=OU.BiasDecay(True,5000,F(99,100),
                         ((0,0,0),(0,0,0),(0,0,0)),em1_2=F(-1,500000))
        with self.assertRaisesRegex(ValueError,'BA exp root detached'):
            X.validate_bias(bad,h)

    def test_readiness_keeps_binary32_correspondence_open(self):
        r=X.readiness()
        self.assertTrue(r['OU_exp_root_real_enclosed_at_same_h_over_tau'])
        self.assertTrue(r['OU_expm1_root_real_enclosed_at_same_h_over_tau'])
        self.assertTrue(r['BA_exp_and_expm1_distinct_arguments_real_enclosed'])
        self.assertFalse(r['exp_expm1_bit_identity_assumed'])
        self.assertFalse(r['binary32_libm_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
