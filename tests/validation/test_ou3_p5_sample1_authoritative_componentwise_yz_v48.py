import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_authoritative_componentwise_yz_v48 as V48


class Sample1AuthoritativeComponentwiseYzV48Tests(unittest.TestCase):
    def test_authoritative_witness_and_targets_are_frozen(self):
        self.assertEqual(V48.WITNESS, (0, 0, 23))
        self.assertEqual(V48.SCHEMA, 4800)
        self.assertEqual(V48.Q_TARGET, 8.0)
        self.assertAlmostEqual(V48.V45.V41_Q_POST, 8.344528951460543, places=12)

    def test_componentwise_caps_do_not_duplicate_yz_radius(self):
        base = {"sample1_full_residual_norm_upper_mps2": 3.0}
        vr = {"total_residual_perturbation_upper_mps2": 0.2}
        ds = {
            "nominal_theta_y_gain_row_norm_upper": 2.0,
            "nominal_theta_z_gain_row_norm_upper": 1.0,
            "theta_y_gain_perturbation_intersected_upper": 0.5,
            "theta_z_gain_perturbation_intersected_upper": 0.25,
        }
        parent = {
            "x_correction_perturbation_abs_upper_rad": 0.1,
            "yz_correction_perturbation_norm_upper_rad": 10.0,
            "total_correction_perturbation_norm_upper_rad": 10.0,
        }
        d = V48._componentwise_yz_caps(
            base=base, vr=vr, ds_detail=ds, parent_caps=parent)
        self.assertAlmostEqual(d["theta_y_component_abs_upper_rad"], 2.0, places=12)
        self.assertAlmostEqual(d["theta_z_component_abs_upper_rad"], 1.0, places=12)
        self.assertLess(d["theta_z_component_abs_upper_rad"],
                        d["theta_y_component_abs_upper_rad"])
        self.assertGreaterEqual(
            d["componentwise_yz_norm_upper_rad"],
            max(d["theta_y_component_abs_upper_rad"],
                d["theta_z_component_abs_upper_rad"]))
        self.assertLessEqual(d["componentwise_yz_norm_upper_rad"], 10.0)
        self.assertLessEqual(d["componentwise_total_norm_upper_rad"], 10.0)

    def test_componentwise_caps_intersect_existing_v31_parents(self):
        base = {"sample1_full_residual_norm_upper_mps2": 3.0}
        vr = {"total_residual_perturbation_upper_mps2": 0.2}
        ds = {
            "nominal_theta_y_gain_row_norm_upper": 2.0,
            "nominal_theta_z_gain_row_norm_upper": 1.0,
            "theta_y_gain_perturbation_intersected_upper": 0.5,
            "theta_z_gain_perturbation_intersected_upper": 0.25,
        }
        parent = {
            "x_correction_perturbation_abs_upper_rad": 0.1,
            "yz_correction_perturbation_norm_upper_rad": 0.3,
            "total_correction_perturbation_norm_upper_rad": 0.4,
        }
        d = V48._componentwise_yz_caps(
            base=base, vr=vr, ds_detail=ds, parent_caps=parent)
        self.assertLessEqual(d["theta_y_component_abs_upper_rad"], 0.3)
        self.assertLessEqual(d["theta_z_component_abs_upper_rad"], 0.3)
        self.assertLessEqual(d["componentwise_yz_norm_upper_rad"], 0.3)
        self.assertLessEqual(d["componentwise_total_norm_upper_rad"], 0.4)

    def test_box_subset_is_componentwise(self):
        parent = [Interval(-2.0, 2.0), Interval(-3.0, 3.0), Interval(-4.0, 4.0)]
        child = [Interval(-1.0, 1.0), Interval(-2.0, 2.0), Interval(-3.0, 3.0)]
        escaped = [Interval(-1.0, 1.0), Interval(-3.1, 2.0), Interval(-3.0, 3.0)]
        self.assertTrue(V48._box_subset(child, parent))
        self.assertFalse(V48._box_subset(escaped, parent))

    def test_norm_helpers_are_conservative(self):
        self.assertGreaterEqual(V48._norm2_up(3.0, 4.0), 5.0)
        self.assertGreaterEqual(V48._norm3_up(1.0, 2.0, 2.0), 3.0)
        self.assertTrue(math.isfinite(V48._norm3_up(1.0, 2.0, 2.0)))


if __name__ == "__main__":
    unittest.main()
