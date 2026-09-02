import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]

# Current source-generated certificate/status producers. Historical replay/local
# basin and retired information-word producers are intentionally excluded.
PRODUCERS = [
    ROOT / "tools" / "ou3_p2_clock_phase_tuner_graph.py",
    ROOT / "tools" / "ou3_p4_augmented_complete_word_design_v5.py",
    ROOT / "tools" / "ou3_p4_joint_first_accel_cover_v3.py",
    ROOT / "tools" / "ou3_p4_terminal_cluster_p2_reduction.py",
    ROOT / "tools" / "ou3_p4_terminal_source_equivalence.py",
]

FORBIDDEN_RESULT_KEYS = {
    "generated_at", "generated_at_utc", "timestamp", "timestamp_utc",
    "build_id", "run_id", "run_number", "run_attempt", "workflow_id",
    "workflow_run_id", "artifact_id", "ci_run_id",
}
FORBIDDEN_CI_ENV = {
    "GITHUB_RUN_ID", "GITHUB_RUN_NUMBER", "GITHUB_RUN_ATTEMPT",
    "GITHUB_WORKFLOW_REF", "GITHUB_WORKFLOW_SHA",
}


def dotted_name(node: ast.AST) -> str:
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
    return ".".join(reversed(parts))


class Ou3ResultDeterminismTests(unittest.TestCase):
    def test_current_certificate_json_producers_have_no_ephemeral_provenance(self):
        for path in PRODUCERS:
            with self.subTest(path=path.name):
                self.assertTrue(path.is_file(), msg=f"current producer missing: {path}")
                text = path.read_text(encoding="utf-8")
                tree = ast.parse(text, filename=str(path))
                for name in FORBIDDEN_CI_ENV:
                    self.assertNotIn(name, text)
                for node in ast.walk(tree):
                    if isinstance(node, ast.Call):
                        fn = dotted_name(node.func)
                        self.assertNotIn(fn, {
                            "time.time", "time.time_ns", "datetime.now", "datetime.utcnow",
                            "datetime.datetime.now", "datetime.datetime.utcnow",
                        }, msg=f"wall-clock call {fn} in {path}")
                    if isinstance(node, ast.Dict):
                        for key in node.keys:
                            if isinstance(key, ast.Constant) and isinstance(key.value, str):
                                self.assertNotIn(
                                    key.value.lower(), FORBIDDEN_RESULT_KEYS,
                                    msg=f"ephemeral result key {key.value!r} in {path}",
                                )

    def test_readme_build_provenance_is_outside_hashed_results(self):
        inside = ROOT / "reports" / "results" / "readme" / "PROVENANCE.md"
        outside = ROOT / "reports" / "readme-results-provenance.md"
        self.assertFalse(inside.exists())
        self.assertTrue(outside.is_file())
        self.assertIn("intentionally outside `reports/results/`", outside.read_text(encoding="utf-8"))

    def test_result_fingerprint_has_no_wall_clock_dependency(self):
        path = ROOT / "tools" / "ou_replay_fingerprint.py"
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text, filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                self.assertNotIn(dotted_name(node.func), {
                    "time.time", "time.time_ns", "datetime.now", "datetime.utcnow",
                    "datetime.datetime.now", "datetime.datetime.utcnow",
                })


if __name__ == "__main__":
    unittest.main()
