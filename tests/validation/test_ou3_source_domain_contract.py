import importlib.util
import math
from pathlib import Path
import struct
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_source_domain_contract", ROOT / "tools" / "ou3_source_domain_contract.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

import ou3_sea3_directional_p2_ha_feasibility as sea3  # noqa: E402


def f32(value):
    return struct.unpack(">f", struct.pack(">f", float(value)))[0]


class SourceDomainContractTests(unittest.TestCase):
    def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(d["schema"], 3)
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["source_complete_parameter_domain"])
        self.assertFalse(d["validated_arithmetic"])
        self.assertFalse(d["outward_rounded"])
        self.assertEqual(d["implementation_scalar_semantics"]["type"], "IEEE754_BINARY32")
        self.assertEqual(d["continuous_parameters"]["tau_aw_s"], [f32(0.02), f32(12.0)])
        self.assertEqual(d["continuous_parameters"]["sigma_aw_mps2"], [f32(0.05), f32(4.0)])
        self.assertEqual(set(d["discrete_source_branches"]["mode"]), {"H", "A"})

    def test_constexpr_arithmetic_rounds_as_binary32_after_each_operation(self):
        text = """
        constexpr float A = 0.1f;
        constexpr float B = A + A;
        constexpr float C = B + A;
        constexpr float DT = 1.0f / 200.0f;
        """
        a = f32(0.1)
        b = f32(a + a)
        c = f32(b + a)
        dt = f32(f32(1.0) / f32(200.0))
        self.assertEqual(mod.parse_const(text, "A"), a)
        self.assertEqual(mod.parse_const(text, "B"), b)
        self.assertEqual(mod.parse_const(text, "C"), c)
        self.assertEqual(mod.parse_const(text, "DT"), dt)
        self.assertNotEqual(dt, 0.005)

    def test_validated_parameter_box_outwardly_contains_every_source_endpoint(self):
        d = mod.build(mod.DEFAULT_HEADER)
        box = d["validated_parameter_box"]
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertEqual(box["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(box["continuous_word_enclosed"])
        self.assertFalse(box["nonlinear_word_enclosed"])
        for name, source_bounds in d["continuous_parameters"].items():
            lo, hi = box["continuous_parameters"][name]
            self.assertLess(lo, source_bounds[0])
            self.assertGreater(hi, source_bounds[1])
            self.assertEqual(lo, math.nextafter(source_bounds[0], -math.inf))
            self.assertEqual(hi, math.nextafter(source_bounds[1], math.inf))
        for name, source_value in d["timing_constants_s"].items():
            lo, hi = box["timing_constants_s"][name]
            self.assertLess(lo, source_value)
            self.assertGreater(hi, source_value)
            self.assertEqual(lo, math.nextafter(source_value, -math.inf))
            self.assertEqual(hi, math.nextafter(source_value, math.inf))

    def test_configured_runtime_sampling_assumption_is_explicit_and_source_bound(self):
        d = mod.build(mod.DEFAULT_HEADER)
        runtime = d["configured_runtime_assumption"]
        expected = f32(f32(1.0) / f32(200.0))
        self.assertEqual(runtime["qualification"], "CONFIGURED_VALIDATION_RUNTIME_ASSUMPTION")
        self.assertEqual(runtime["sample_period_contract"], "FIXED_SOURCE_NOMINAL")
        self.assertEqual(runtime["imu_dt_s"], expected)
        self.assertFalse(runtime["api_enforces_this_bound"])
        lo, hi = runtime["imu_dt_outward_interval_s"]
        self.assertLess(lo, expected)
        self.assertGreater(hi, expected)
        self.assertEqual(d["validated_parameter_box"]["configured_runtime"], runtime)

    def test_contract_names_every_hybrid_transition_required_for_deployment(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(
            set(d["hybrid_obligations"]),
            {
                "startup_handoff", "held_to_active", "magnetic_lock",
                "magnetic_regauge_refinement", "tilt_reset", "tilt_relock",
                "cooldown_reentry", "periodic_aw_covariance_sync",
            },
        )
        self.assertEqual(d["periodic_aw_covariance_sync_proof"]["required_mode"], "PSD_NONEXPANSIVE")

    def test_sea3_rao_family_right_inclusion_is_part_of_the_canonical_source_contract(self):
        """`ou3-proof` executes this file, so the RAO bridge has no side workflow."""
        d = sea3.build_inclusion()
        self.assertEqual(sea3.validate(d), [])

        r = d["response_enclosure"]
        self.assertTrue(r["single_worst_envelope_proves_entire_parameter_box_by_monotonicity"])
        self.assertFalse(r["single_nominal_RAO_used"])
        self.assertFalse(r["finite_RAO_grid_used"])
        self.assertTrue(r["uniform_moment_theorem"]["unbanded_acceleration_moment_finite"])
        self.assertGreater(r["acceleration_moment_tightening_vs_flat_6Hz_corner_lower"], 624.0)

        p = d["p2_inclusion"]
        self.assertEqual(p["SEA3_TO_P2_INCLUSION_CERTIFICATE"], "PASS")
        self.assertTrue(p["Lhat_SEA3_subset_L_current_source"])
        self.assertFalse(p["single_RAO_selected_for_inclusion"])
        self.assertFalse(p["P2_pruned_by_SEA3"])
        self.assertEqual(p["P2_physical_source_states"], 800)


if __name__ == "__main__":
    unittest.main()
