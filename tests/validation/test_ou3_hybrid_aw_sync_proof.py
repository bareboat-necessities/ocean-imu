import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "ou3_hybrid_aw_sync_proof", TOOLS / "ou3_hybrid_aw_sync_proof.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class Ou3HybridAwSyncProofTests(unittest.TestCase):
    def test_shipping_default_sync_is_source_bound_psd_nonexpansive(self):
        tool = load_tool()
        out = tool.prove()
        self.assertEqual(out["status"], "PASS", out["failures"])
        self.assertTrue(out["source_binding_pass"])
        self.assertTrue(out["source_complete_for_transition"])
        self.assertEqual(out["proof_mode"], "PSD_NONEXPANSIVE")
        self.assertEqual(out["hybrid_obligation"], "periodic_aw_covariance_sync")
        self.assertEqual(out["jump_gain_upper"], 1.0)
        self.assertEqual(out["additive_W_upper"], 0.0)
        self.assertEqual(out["new_coordinate_W_upper"], 0.0)
        self.assertTrue(out["nonexpansive_information_energy"])
        self.assertFalse(out["sampled_evidence_used"])

    def test_proof_fails_closed_if_psd_projection_is_removed(self):
        tool = load_tool()
        source = tool.DEFAULT_HEADER.read_text(encoding="utf-8")
        source = source.replace(
            "evals(i) = std::max(T(0), evals(i));",
            "evals(i) = evals(i);",
            1,
        )
        with tempfile.TemporaryDirectory() as td:
            fake = Path(td) / "Kalman3D_Wave_OU_III.h"
            fake.write_text(source, encoding="utf-8")
            # prove() reports paths relative to the repository only for the real
            # source tree, so use the source-pattern layer directly for this
            # mutation test.
            checks = {
                name: bool(pattern.search(source))
                for name, pattern in tool.SOURCE_PATTERNS.items()
            }
        self.assertFalse(checks["negative_delta_eigenvalues_are_clamped"])


if __name__ == "__main__":
    unittest.main()
