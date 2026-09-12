"""Source-angle trig enclosure regressions; not deployment/stability evidence."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_source_bound_attitude_trig as X


S1=F(381005220876,452784739765)
C1=F(244640638957,452784739765)
SH=F(452066294520,942933277681)
CH=F(827501801519,942933277681)
S2=F(3994403964948,4392846440677)
C2=F(-1828069149725,4392846440677)


class Tests(unittest.TestCase):
    def test_global_enclosure_contains_close_rational_unit_circle_witnesses(self):
        full=ATT.TrigWitness(F(1),S1,C1,F(1))
        half=ATT.TrigWitness(F(1,2),SH,CH,F(1))
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

    def test_angle_two_radians_no_longer_needs_local_retention_assumption(self):
        full=ATT.TrigWitness(F(1),S2,C2,F(1,2))
        half=ATT.TrigWitness(F(1,2),S1,C1,F(1,2))
        r=ATT.AngularRuntime((2,0,0),F(1),full,half)
        self.assertIs(X.validate(r),r)
        self.assertEqual(X.validate_witness(r.w,full),F(2))
        self.assertEqual(X.validate_witness(r.w,half),F(1))

    def test_global_enclosure_has_rigorous_small_width_at_large_concrete_angle(self):
        slo,shi,clo,chi=X.trig_enclosure(F(10))
        self.assertLessEqual(shi-slo,2*X.TRIG_TOL)
        self.assertLessEqual(chi-clo,2*X.TRIG_TOL)

    def test_small_rate_branch_consumes_no_trig_and_passes(self):
        r=ATT.AngularRuntime((0,0,0),F(1,200))
        self.assertIs(X.validate(r),r)

    def test_readiness_removes_local_guard_but_keeps_binary32_blocker(self):
        r=X.readiness()
        self.assertTrue(r['trig_full_half_angles_derived_from_same_angular_rate_and_step'])
        self.assertTrue(r['detached_unit_circle_points_rejected'])
        self.assertTrue(r['global_finite_rational_angle_enclosure_available'])
        self.assertFalse(r['one_radian_local_angle_guard_required'])
        self.assertTrue(r['one_radian_guard_retained_for_every_admitted_prefix'])
        self.assertFalse(r['binary32_sqrt_div_sin_cos_correspondence_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
