import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_hybrid_contract", ROOT / "tools" / "ou3_hybrid_contract.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def row(kind: str) -> dict:
    return {
        "kind": kind,
        "destination_mode": "H",
        "source_level_W_upper": 2.0,
        "jump_gain_upper": 0.5,
        "additive_W_upper": 0.2,
        "new_coordinate_W_upper": 0.0,
        "destination_level_W": 5.0,
        "inward_margin_lower": 3.8,
        "pass": True,
    }


def complete_check() -> dict:
    # Use the legacy names for regauge/cooldown deliberately.  The independent
    # theorem gate must normalize them to the current source-domain contract.
    kinds = [
        "startup_handoff",
        "held_to_active",
        "magnetic_lock",
        "magnetic_regauge",
        "tilt_reset",
        "tilt_relock",
        "cooldown",
    ]
    return {"hybrid": {"pass": False, "bounds": [row(k) for k in kinds]}}


class HybridContractTests(unittest.TestCase):
    def test_source_domain_obligations_are_exact(self):
        self.assertEqual(mod.REQUIRED, {
            "startup_handoff",
            "held_to_active",
            "magnetic_lock",
            "magnetic_regauge_refinement",
            "tilt_reset",
            "tilt_relock",
            "cooldown_reentry",
            "periodic_aw_covariance_sync",
        })

    def test_legacy_names_are_normalized_and_aw_sync_is_analytic(self):
        out = mod.validate(complete_check())
        self.assertTrue(out["pass"], out["failures"])
        self.assertFalse(out["sampled_evidence_used"])
        self.assertEqual(out["missing"], [])
        self.assertTrue(out["periodic_aw_covariance_sync"]["pass"])
        self.assertFalse(out["periodic_aw_covariance_sync"]["strict_inward_margin_required"])
        self.assertEqual(out["periodic_aw_covariance_sync"]["jump_gain_upper"], 1.0)
        self.assertEqual(out["periodic_aw_covariance_sync"]["additive_W_upper"], 0.0)

    def test_missing_source_obligation_fails_closed(self):
        check = complete_check()
        check["hybrid"]["bounds"] = [
            x for x in check["hybrid"]["bounds"] if x["kind"] != "magnetic_lock"
        ]
        out = mod.validate(check)
        self.assertFalse(out["pass"])
        self.assertIn("magnetic_lock", out["missing"])

    def test_aw_sync_source_binding_failure_fails_closed(self):
        bad = {
            "qualification": "SOURCE_BOUND_ANALYTIC_HYBRID_PROOF",
            "sampled_evidence_used": False,
            "source_binding_pass": False,
            "proof_mode": "PSD_NONEXPANSIVE",
            "nonexpansive_information_energy": False,
            "jump_gain_upper": None,
            "additive_W_upper": None,
            "new_coordinate_W_upper": None,
        }
        out = mod.validate(complete_check(), aw_proof=bad)
        self.assertFalse(out["pass"])
        self.assertIn("periodic_aw_covariance_sync", out["missing"])


if __name__ == "__main__":
    unittest.main()
