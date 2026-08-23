import ast
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]

PRODUCERS = [
    ROOT / "tools" / "ou3_exact_replay.py",
    ROOT / "tools" / "ou3_numerical_certificate.py",
    ROOT / "tools" / "ou3_information_certificate.py",
    ROOT / "tools" / "ou3_information_completion.py",
    ROOT / "tools" / "ou3_information_enclosure_contract.py",
    ROOT / "tools" / "ou3_certificate_completion.py",
    ROOT / "tools" / "ou3_validate_enclosure.py",
]

# Scientific result JSONs under reports/results must be functions only of the
# repository code/configuration and the simulation inputs. CI/run provenance
# belongs outside reports/results and must never enter these payloads.
FORBIDDEN_RESULT_KEYS = {
    "generated_at",
    "generated_at_utc",
    "timestamp",
    "timestamp_utc",
    "build_id",
    "run_id",
    "run_number",
    "run_attempt",
    "workflow_id",
    "workflow_run_id",
    "artifact_id",
    "ci_run_id",
}
FORBIDDEN_CI_ENV = {
    "GITHUB_RUN_ID",
    "GITHUB_RUN_NUMBER",
    "GITHUB_RUN_ATTEMPT",
    "GITHUB_WORKFLOW_REF",
    "GITHUB_WORKFLOW_SHA",
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
    def test_certificate_json_producers_have_no_ephemeral_provenance(self):
        for path in PRODUCERS:
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8")
                tree = ast.parse(text, filename=str(path))

                # Direct references to CI execution identity are forbidden even
                # if they are not currently written to JSON.
                for name in FORBIDDEN_CI_ENV:
                    self.assertNotIn(name, text)

                for node in ast.walk(tree):
                    # Forbid wall-clock generation in scientific result tools.
                    if isinstance(node, ast.Call):
                        fn = dotted_name(node.func)
                        self.assertNotIn(
                            fn,
                            {
                                "time.time",
                                "time.time_ns",
                                "datetime.now",
                                "datetime.utcnow",
                                "datetime.datetime.now",
                                "datetime.datetime.utcnow",
                            },
                            msg=f"wall-clock call {fn} in {path}",
                        )

                        # JSON output is canonical with stable key ordering.
                        if fn in {"json.dump", "json.dumps"}:
                            kw = {k.arg: k.value for k in node.keywords if k.arg}
                            self.assertIn("sort_keys", kw, msg=f"noncanonical JSON in {path}")
                            self.assertIsInstance(kw["sort_keys"], ast.Constant)
                            self.assertIs(kw["sort_keys"].value, True)

                    # Explicit ephemeral result fields are forbidden. Simulation
                    # time fields such as time_s are intentionally not banned.
                    if isinstance(node, ast.Dict):
                        for key in node.keys:
                            if isinstance(key, ast.Constant) and isinstance(key.value, str):
                                self.assertNotIn(
                                    key.value.lower(),
                                    FORBIDDEN_RESULT_KEYS,
                                    msg=f"ephemeral result key {key.value!r} in {path}",
                                )

    def test_readme_build_provenance_is_outside_hashed_results(self):
        inside = ROOT / "reports" / "results" / "readme" / "PROVENANCE.md"
        outside = ROOT / "reports" / "readme-results-provenance.md"
        self.assertFalse(
            inside.exists(),
            msg="build/run provenance must not be part of the reports/results fingerprint",
        )
        self.assertTrue(outside.is_file())
        text = outside.read_text(encoding="utf-8")
        self.assertIn("intentionally outside `reports/results/`", text)

    def test_result_fingerprint_has_no_wall_clock_dependency(self):
        path = ROOT / "tools" / "ou_replay_fingerprint.py"
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text, filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            self.assertNotIn(
                dotted_name(node.func),
                {
                    "time.time",
                    "time.time_ns",
                    "datetime.now",
                    "datetime.utcnow",
                    "datetime.datetime.now",
                    "datetime.datetime.utcnow",
                },
            )


if __name__ == "__main__":
    unittest.main()
