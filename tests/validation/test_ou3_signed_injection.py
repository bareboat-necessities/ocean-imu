"""Signed world-injection transport (Lemma I*) and the rotating-frame floor."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.signed_injection import (  # noqa: E402
    certificate, feasibility, injection_rotation_identity, relaxed_reset_example_exclusion,
    rotating_frame_example, rotation_angle, third_order_reset_remainder)
from tools.stability.ou3_theorem.matrix_certificates import identity  # noqa: E402
from tools.stability.ou3_theorem.world_frame import quaternion_rotation  # noqa: E402


class SignedInjectionTest(unittest.TestCase):
    def test_committed_certificate_matches_exact_reproduction(self):
        path = ROOT/'reports/results/ou3_stability/signed-injection-certificate.json'
        record = json.loads(path.read_text())
        self.assertEqual(record, certificate())
        self.assertFalse(record['half_angle_product_remainder_source_bound'])
        self.assertFalse(record['source_uniform_signed_injection_bound'])
        self.assertFalse(record['theorem_closed'])
        self.assertFalse(record['lemma_I_star']['norm_sum_of_injections_used'])

    def test_ordered_product_identity_on_supplied_words(self):
        quats = [(9, 1, -2, 1), (20, 0, 3, -1), (7, 2, 2, 1), (40, -1, 1, 3), (11, 1, 0, 0)]
        rots = [quaternion_rotation(q) for q in quats]
        for pattern in (('rate', 'inject', 'rate', 'inject', 'inject'),
                        ('inject', 'inject', 'rate', 'rate', 'inject'),
                        ('rate', 'rate', 'rate', 'rate', 'inject')):
            word = list(zip(pattern, rots))
            result = injection_rotation_identity(word, quaternion_rotation((5, 1, 1, -1)))
            self.assertLessEqual(result['injection_angle'], result['bound'])
        with self.assertRaises(ValueError):
            injection_rotation_identity([('inject', [[F(2), 0, 0], [0, 1, 0], [0, 0, 1]])], identity(3))

    def test_single_injection_is_bounded_by_adjacent_attitude_errors(self):
        m0 = quaternion_rotation((3, 1, 0, 1))
        x = quaternion_rotation((2, -1, 1, 0))
        result = injection_rotation_identity([('inject', x)], m0)
        self.assertLessEqual(rotation_angle(x), result['bound'])
        exclusion = relaxed_reset_example_exclusion()
        self.assertGreater(float(exclusion['required_attitude_error_rad_lower']), 1.14)

    def test_third_order_reset_factor(self):
        for v in ([F(1, 2), F(0), F(0)], [F(1, 3), F(1, 3), F(-1, 3)], [F(0), F(-9, 10), F(1, 5)]):
            remainder, bound = third_order_reset_remainder(v)
            self.assertLessEqual(remainder, bound)

    def test_signed_rotations_cancel_and_one_sided_do_not(self):
        example = rotating_frame_example()
        self.assertLess(F(example['one_sided_floor']), F(example['signed_floor']))
        self.assertLessEqual(F(example['signed_floor']), F(example['no_rotation_floor']))

    def test_feasibility_classification(self):
        rows = feasibility()
        retained = [r for r in rows if r['case'].startswith('retained')]
        self.assertFalse(any(r['perturbative_feasible'] for r in retained))
        invariant = [r for r in rows if r['case'].endswith('invariant only')]
        self.assertFalse(any(r['rotating_frame_feasible'] for r in invariant))
        exact_bias = [r for r in rows if r['case'] == 'vanishing radius, bias exact']
        self.assertTrue(all(r['rotating_frame_feasible'] for r in exact_bias))


if __name__ == '__main__':
    unittest.main()
