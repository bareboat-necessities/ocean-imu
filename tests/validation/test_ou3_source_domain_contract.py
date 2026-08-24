import importlib.util
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
