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


def f32(value):
    return struct.unpack(">f", struct.pack(">f", float(value)))[0]


class SourceDomainContractTests(unittest.TestCase):
    def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(d["schema"], 2)
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["source_complete_parameter_domain"])
        self.assertFalse(d["validated_arithmetic"])
        self.assertFalse(d["outward_rounded"])
        self.assertEqual(
            d["implementation_scalar_semantics"]["type"], "IEEE754_BINARY32"
        )
        self.assertEqual(
            d["continuous_parameters"]["tau_aw_s"], [f32(0.02), f32(12.0)]
        )
        # apply_ou_tune_ enforces max(0.05f, band_noise_floor_sigma_()) before
        # writing the OU stationary standard deviation.  The contract records
        # the actually deployed binary32 floor, not the decimal source token.
        self.assertEqual(
            d["continuous_parameters"]["sigma_aw_mps2"], [f32(0.05), f32(6.0)]
        )
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

    def test_deployment_step_domain_exposes_missing_finite_upper_guard(self):
        d = mod.build(mod.DEFAULT_HEADER)
        step = d["accepted_update_step_domain_s"]
        self.assertEqual(step["lower_open"], 0.0)
        self.assertIsNone(step["upper"])
        self.assertFalse(step["source_complete_finite_upper_bound"])
        self.assertEqual(
            d["validated_ou_primitive_backend"]["theorem_promotion"],
            "BLOCKED_BY_UNBOUNDED_ACCEPTED_DT",
        )

    def test_rational_taylor_ou_primitive_box_contains_direct_evaluations(self):
        # This is the range the proof backend will use once deployment supplies
        # a finite accepted-step upper guard.  It deliberately spans the full
        # current tau safety range, including the x=h/tau=12.5 corner.
        box = mod.validated_ou_primitives((0.001, 0.25), (0.02, 12.0))
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertLessEqual(box["alpha"][0], box["alpha"][1])
        self.assertLessEqual(box["phi_pa_s2"][0], box["phi_pa_s2"][1])
        self.assertLessEqual(box["phi_Sa_s3"][0], box["phi_Sa_s3"][1])

        for h in (0.001, 0.005, 0.05, 0.25):
            for tau in (0.02, 0.1, 1.0, 12.0):
                x = h / tau
                alpha = math.exp(-x)
                phi_pa = tau * tau * (x + math.expm1(-x))
                phi_sa = tau ** 3 * (0.5 * x * x - x - math.expm1(-x))
                self.assertLessEqual(box["alpha"][0], alpha)
                self.assertGreaterEqual(box["alpha"][1], alpha)
                self.assertLessEqual(box["phi_pa_s2"][0], phi_pa)
                self.assertGreaterEqual(box["phi_pa_s2"][1], phi_pa)
                self.assertLessEqual(box["phi_Sa_s3"][0], phi_sa)
                self.assertGreaterEqual(box["phi_Sa_s3"][1], phi_sa)

    def test_contract_names_every_hybrid_transition_required_for_deployment(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(
            set(d["hybrid_obligations"]),
            {
                "startup_handoff",
                "held_to_active",
                "magnetic_lock",
                "magnetic_regauge_refinement",
                "tilt_reset",
                "tilt_relock",
                "cooldown_reentry",
                "periodic_aw_covariance_sync",
            },
        )
        self.assertEqual(
            d["periodic_aw_covariance_sync_proof"]["required_mode"],
            "PSD_NONEXPANSIVE",
        )


if __name__ == "__main__":
    unittest.main()
