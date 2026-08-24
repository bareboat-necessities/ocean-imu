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
    out = {
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
    if kind == "tilt_reset":
        out.update({
            "discarded_pre_reset_tilt_excluded_from_multiplicative_gain": True,
            "reset_to_funnel_exact_map": True,
        })
    if kind == "cooldown_reentry":
        out.update({
            "reachable_word_product_used": True,
            "global_worst_word_power_used": False,
        })
    return out


def complete_check() -> dict:
    kinds = [
        "startup_handoff",
        "held_to_active",
        "magnetic_lock",
        "magnetic_regauge_refinement",
        "tilt_reset",
        "tilt_relock",
        "cooldown_reentry",
    ]
    return {"hybrid": {"pass": False, "bounds": [row(k) for k in kinds]}}


class HybridContractTests(unittest.TestCase):
    def test_source_domain_obligations_are_exact(self):
        self.assertEqual(mod.REQUIRED, {
            "startup_handoff", "held_to_active", "magnetic_lock",
            "magnetic_regauge_refinement", "tilt_reset", "tilt_relock",
            "cooldown_reentry", "periodic_aw_covariance_sync",
        })

    def test_current_names_close_and_aw_sync_is_analytic(self):
        out = mod.validate(complete_check())
        self.assertTrue(out["pass"], out["failures"])
        self.assertFalse(out["sampled_evidence_used"])
        self.assertFalse(out["legacy_name_aliases_used"])
        self.assertEqual(out["missing"], [])
        self.assertTrue(out["periodic_aw_covariance_sync"]["pass"])
        self.assertFalse(out["periodic_aw_covariance_sync"]["strict_inward_margin_required"])
        self.assertEqual(out["periodic_aw_covariance_sync"]["jump_gain_upper"], 1.0)
        self.assertEqual(out["periodic_aw_covariance_sync"]["additive_W_upper"], 0.0)

    def test_legacy_names_do_not_satisfy_current_obligations(self):
        check = complete_check()
        for r in check["hybrid"]["bounds"]:
            if r["kind"] == "magnetic_regauge_refinement":
                r["kind"] = "magnetic_regauge"
            if r["kind"] == "cooldown_reentry":
                r["kind"] = "cooldown"
        out = mod.validate(check)
        self.assertFalse(out["pass"])
        self.assertIn("magnetic_regauge_refinement", out["missing"])
        self.assertIn("cooldown_reentry", out["missing"])

    def test_tilt_reset_must_discard_rewritten_tilt_energy(self):
        check = complete_check()
        tilt = next(x for x in check["hybrid"]["bounds"] if x["kind"] == "tilt_reset")
        tilt["discarded_pre_reset_tilt_excluded_from_multiplicative_gain"] = False
        out = mod.validate(check)
        self.assertFalse(out["pass"])
        self.assertIn("tilt_reset", out["missing"])

    def test_cooldown_must_use_reachable_word_products(self):
        check = complete_check()
        cool = next(x for x in check["hybrid"]["bounds"] if x["kind"] == "cooldown_reentry")
        cool["reachable_word_product_used"] = False
        cool["global_worst_word_power_used"] = True
        out = mod.validate(check)
        self.assertFalse(out["pass"])
        self.assertIn("cooldown_reentry", out["missing"])

    def test_missing_source_obligation_fails_closed(self):
        check = complete_check()
        check["hybrid"]["bounds"] = [x for x in check["hybrid"]["bounds"] if x["kind"] != "magnetic_lock"]
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
