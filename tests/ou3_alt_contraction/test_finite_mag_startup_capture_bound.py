"""Deterministic startup magnetic capture-bound regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_capture_bound as X

class Tests(unittest.TestCase):
    def test_unrestricted_frame_image_bound_does_not_claim_small_angle_capture(self):
        b=X.unrestricted_frame_bound()
        self.assertEqual(b.half_tilt_sin_max,1)
        self.assertEqual(b.earth_rotation_error_max,150)
        self.assertEqual(b.mean_perturbation_max,157)
        self.assertFalse(b.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'):
            X.require_capture(b)
        r=X.readiness()
        self.assertTrue(r['unrestricted_full_frame_chord_bound_closed'])
        self.assertFalse(r['small_frame_accuracy_required_by_finite_atlas_word'])

    def test_current_envelope_exact_open_tilt_margin_is_four_over_75(self):
        self.assertEqual(X.default_tilt_margin_half_sin(),F(4,75))
        at=X.derive(F(4,75))
        self.assertEqual(at.earth_rotation_error_max,8)
        self.assertEqual(at.mean_perturbation_max,15)
        self.assertFalse(at.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'): X.require_capture(at)

    def test_explicit_full_frame_bound_supplies_conditional_capture_algebra(self):
        b=X.derive(F(1,100))
        self.assertEqual(X.DECLARED_STARTUP_TILT_RAD_MAX,F(1,50))
        self.assertEqual(b.half_tilt_sin_max,F(1,100))
        self.assertEqual(b.earth_rotation_error_max,F(3,2))
        self.assertEqual(b.mean_perturbation_max,F(17,2))
        self.assertEqual(b.sin_yaw_error_max,F(17,30))
        self.assertTrue(b.capture_nonzero); self.assertTrue(X.require_capture(b))

    def test_full_attitude_entrance_requires_explicit_frame_and_handoff_bounds(self):
        c=X.conditional_fresh_attitude_certificate(half_frame_sin_max=F(1,100),
                                                 handoff_tilt_rad_max=F(1,50))
        self.assertEqual(c.yaw_rad_upper,F(61,100))
        self.assertEqual(c.tilt_rad_upper,F(1,50))
        self.assertEqual(c.full_attitude_rad_upper,F(63,100))
        self.assertGreater(c.sin_yaw_test_lower,F(17,30))
        self.assertLess(c.full_attitude_rad_upper,F(3,4))
        self.assertTrue(c.below_pi_over_4)

    def test_perfect_tilt_still_retains_full_deterministic_hardiron_and_noise(self):
        b=X.derive(0)
        self.assertEqual(b.mean_perturbation_max,7)
        self.assertEqual(b.sin_yaw_error_max,F(7,15))
        self.assertTrue(b.capture_nonzero)
        self.assertNotEqual(b.mean_perturbation_max,F(7,250))

    def test_no_fake_capture_at_large_tilt(self):
        b=X.derive(F(3,50))
        self.assertEqual(b.earth_rotation_error_max,9)
        self.assertEqual(b.mean_perturbation_max,16)
        self.assertFalse(b.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'): X.require_capture(b)

    def test_tilt_only_bound_cannot_bound_heading_stripped_accumulation_frame(self):
        from tools.stability.ou3_alt_contraction.finite_sensor_source_runtime import q_rotate
        # Perfect tilt at both headings, same norm-25 uT Earth field, no noise.
        field=(F(15),F(0),F(20))
        south=q_rotate((0,0,0,1),field)
        north=q_rotate((1,0,0,0),field)
        mean=tuple((x+y)/2 for x,y in zip(south,north))
        self.assertEqual(mean,(0,0,20))
        self.assertEqual(sum((x-y)**2 for x,y in zip(south,field)),900)
        self.assertGreater(900,X.derive(0).mean_perturbation_max**2)
        # This disproves the frame-bound implication, not a complete shipping
        # capture trace. A temporal source theorem remains necessary.
        with self.assertRaises(TypeError): X.conditional_fresh_attitude_certificate()

    def test_readiness_distinguishes_conditional_algebra_from_source_capture(self):
        r=X.readiness()
        self.assertTrue(r['deterministic_average_does_not_claim_sqrtN_improvement'])
        self.assertTrue(r['current_envelope_half_tilt_sin_margin_is_4_over_75'])
        self.assertTrue(r['conditional_full_frame_0p02_and_handoff_tilt_0p02_imply_pi_over_4_entry'])
        for key in ('full_accumulation_to_handoff_frame_bound_source_qualified',
                    'declared_startup_tilt_0p02_rad_attached_via_sin_x_le_x',
                    'real_arithmetic_yaw_lt_0p61_rad_certified',
                    'real_arithmetic_full_SO3_lt_pi_over_4_certified',
                    'atan2_AngleAxis_binary32_correspondence_attached',
                    'startup_capture_closed','ALT_STARTUP_PASS'):
            self.assertFalse(r[key])

if __name__=='__main__': unittest.main()
