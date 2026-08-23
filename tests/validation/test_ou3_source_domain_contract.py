import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_source_domain_contract", ROOT / "tools" / "ou3_source_domain_contract.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class SourceDomainContractTests(unittest.TestCase):
    def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["source_complete_parameter_domain"])
        self.assertFalse(d["validated_arithmetic"])
        self.assertFalse(d["outward_rounded"])
        self.assertEqual(d["continuous_parameters"]["tau_aw_s"], [0.02, 12.0])
        self.assertEqual(d["continuous_parameters"]["sigma_aw_mps2"], [0.0, 6.0])
        self.assertIn("periodic_aw_covariance_sync", d["hybrid_obligations"])
        self.assertEqual(set(d["discrete_source_branches"]["mode"]), {"H", "A"})


if __name__ == "__main__":
    unittest.main()
