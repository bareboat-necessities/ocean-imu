"""Source-angle trig enclosure regressions; not deployment/stability evidence."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_source_bound_attitude_trig as X


class Tests(unittest.TestCase):
    def test_taylor_enclosure_contains_rational_unit_circle_witnesses_at_actual_angles(self):
        # Rational unit-circle points from u=6/11 and u=12/47.  They lie in
        # the rigorous alternating-series boxes at theta=1 and theta=1/2.
        full=ATT.TrigWitness(F(1),F(132,157),F(85,157),F(1))
        half=ATT.TrigWitness(F(1,2),F(1128,2353),F(2065,2353),F(1))
        r=ATT.AngularRuntime((1,0,0),F(1),full,half)
        self.assertIs(X.validate(r),r)
        self.assertEqual(X.validate_witness(r.w,full),F(1))
        self.assertEqual(X.validate_witness(r.w,half),F(1,2))

    def test_detached_unit_circle_point_is_rejected(self):
        full=ATT.TrigWitness(F(1),F(0),F(1),F(1))
        half=ATT.TrigWitness(F(1,2),F(0),F(1),F(1))
        r=ATT.AngularRuntime((1,0,0),F(1),full,half)
        with self.assertRaisesRegex(ValueError,'sine witness detached'):
            X.validate(r)

    def test_large_angle_fails_closed_until_retention_is_proved(self):
        full=ATT.TrigWitness(F(1),F(0),F(1),F(1,2))
        half=ATT.TrigWitness(F(1,2),F(0),F(1),F(1,2))
        r=ATT.AngularRuntime((2,0,0),F(1),full,half)
        with self.assertRaisesRegex(ValueError,'outside certified <=1 rad'):
            X.validate(r)

    def test_small_rate_branch_consumes_no_trig_and_passes(self):
        r=ATT.AngularRuntime((0,0,0),F(1,200))
        self.assertIs(X.validate(r),r)

    def test_readiness_exposes_retention_and_binary32_blockers(self):
        r=X.readiness()
        self.assertTrue(r['trig_full_half_angles_derived_from_same_angular_rate_and_step'])
        self.assertTrue(r['detached_unit_circle_points_rejected'])
        self.assertFalse(r['one_radian_guard_retained_for_every_admitted_prefix'])
        self.assertFalse(r['binary32_sqrt_div_sin_cos_correspondence_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
