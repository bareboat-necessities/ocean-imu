from __future__ import annotations

import ast
import importlib
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
STABILITY = TOOLS / "stability"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(STABILITY))

RETIRED_MODULES = {
    "ou3_p3_p2_v1_history_frontier",
    "ou3_p3_p2_v1_stage_phase_translation",
    "ou3_p3_p2_v1_full_state_join",
    "ou3_p3_canonical_gate",
    "ou3_p4_p3_metric_attachment",
    "ou3_p4_accelerometer_corotated_aw",
    "ou3_p4_source_word_timing",
    "ou3_p4_vector_remainder_sector",
    "ou3_p4_signed_joseph_feasibility",
    "ou3_p4_canonical_gate",
    "ou3_p4_source_path_reachability",
    "ou3_p4_source_node_cells",
    "ou3_p4_sample_clock_source_refinement",
}

STALE_RELOCATED_PATH = re.compile(
    r"(?:REPO|ROOT|repo_root)\s*/\s*['\"]tools['\"]\s*/\s*"
    r"['\"]ou3_[^'\"]+\.py['\"]"
)

OPERATING_DOMAIN_BASENAME = "ou3_proof_" + "operating_domain.json"
DIRECTIONAL_RESPONSE_BASENAME = "ou3_sea3_directional_" + "response_domain.json"
STABILITY_JSON_BASENAMES = {
    OPERATING_DOMAIN_BASENAME,
    DIRECTIONAL_RESPONSE_BASENAME,
    "ou3_sea3_spectral_" + "moment_bridge.json",
}
TEXT_SUFFIXES = {".py", ".yml", ".yaml", ".md", ".tex", ".sh", ".txt"}


def _stale_stability_json_references() -> list[str]:
    stale: list[str] = []
    roots = (
        STABILITY,
        ROOT / ".github" / "workflows",
        ROOT / "tests" / "validation",
        ROOT / "doc",
        ROOT / "docs",
    )
    for base in roots:
        if not base.exists():
            continue
        for path in sorted(p for p in base.rglob("*") if p.is_file()):
            if path.suffix not in TEXT_SUFFIXES and path.name != "Makefile":
                continue
            try:
                lines = path.read_text(encoding="utf-8").splitlines()
            except UnicodeDecodeError:
                continue
            for lineno, line in enumerate(lines, 1):
                for basename in STABILITY_JSON_BASENAMES:
                    if basename not in line:
                        continue
                    # Every retained reference must make the new stability
                    # package location explicit.  The basename is constructed
                    # above so this integrity test does not match itself.
                    if "stability" not in line:
                        stale.append(
                            f"{path.relative_to(ROOT)}:{lineno}: {line.strip()}"
                        )
    return stale


class StabilityPackageIntegrityTests(unittest.TestCase):
    def test_all_retained_ou3_imports_resolve_inside_retained_tree(self):
        modules = {p.stem for p in STABILITY.glob("ou3_*.py")}
        top = {p.stem for p in TOOLS.glob("ou3_*.py")}
        unresolved = []
        retired = []
        for path in sorted(STABILITY.glob("ou3_*.py")):
            tree = ast.parse(path.read_text(), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names = [a.name.split(".")[0] for a in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    names = [node.module.split(".")[0]]
                else:
                    names = []
                for name in names:
                    if not name.startswith("ou3_"):
                        continue
                    if name in RETIRED_MODULES:
                        retired.append(f"{path.name}: {name}")
                    if name not in modules and name not in top:
                        unresolved.append(f"{path.name}: {name}")
        self.assertEqual(
            [], retired, "retired stability imports returned:\n" + "\n".join(retired)
        )
        self.assertEqual(
            [], unresolved, "unresolved stability imports:\n" + "\n".join(unresolved)
        )

    def test_every_retained_stability_module_imports(self):
        failures = []
        for path in sorted(STABILITY.glob("ou3_*.py")):
            try:
                importlib.import_module(path.stem)
            except Exception as exc:  # pragma: no cover - failure-reporting path
                failures.append(f"{path.name}: {type(exc).__name__}: {exc}")
        self.assertEqual(
            [], failures, "retained stability module import failures:\n" + "\n".join(failures)
        )

    def test_retained_stability_modules_do_not_use_pre_package_python_paths(self):
        stale = []
        for path in sorted(STABILITY.glob("ou3_*.py")):
            for match in STALE_RELOCATED_PATH.finditer(path.read_text(encoding="utf-8")):
                stale.append(f"{path.name}: {match.group(0)}")
        self.assertEqual(
            [],
            stale,
            "retained stability modules still use pre-package Python paths:\n"
            + "\n".join(stale),
        )

    def test_stability_json_domains_live_only_in_stability_package(self):
        missing = [
            name for name in sorted(STABILITY_JSON_BASENAMES)
            if not (STABILITY / name).is_file()
        ]
        root_duplicates = [
            name for name in sorted(STABILITY_JSON_BASENAMES)
            if (TOOLS / name).exists()
        ]
        stale = _stale_stability_json_references()
        self.assertEqual([], missing, "missing stability JSON domains: " + ", ".join(missing))
        # Report all stale consumers while the temporary root copies still
        # exist; deleting those copies first would turn this into a sequence of
        # opaque FileNotFoundError failures instead of one actionable list.
        self.assertEqual(
            [], stale,
            "stale pre-package stability JSON references:\n" + "\n".join(stale),
        )
        self.assertEqual(
            [], root_duplicates,
            "stability JSON domains still duplicated at tools/: " + ", ".join(root_duplicates),
        )

    def test_all_stability_json_files_have_moved(self):
        self.assertEqual([], sorted(p.name for p in TOOLS.glob("ou3_*.json")))

    def test_no_temporary_relocation_workflows_remain(self):
        workflows = ROOT / ".github" / "workflows"
        self.assertEqual([], sorted(p.name for p in workflows.glob("ou3-stability-*.yml")))

    def test_inline_workflow_imports_use_stability_package(self):
        # Inline Python is not covered by importing retained .py modules.
        stale = []
        pattern = re.compile(r"^\s*(?:import|from)\s+ou3_", re.MULTILINE)
        for path in sorted((ROOT / ".github" / "workflows").glob("*.yml")):
            text = path.read_text(encoding="utf-8")
            if pattern.search(text) and "tools/stability" not in text:
                stale.append(path.name)
            if pattern.search(text) and re.search(
                r"sys\.path\.insert\(0,\s*['\"]tools['\"]\)", text
            ):
                stale.append(path.name + ": root-only inline import path")
        self.assertEqual([], stale)

    def test_canonical_workflow_watches_stability_json(self):
        text = (ROOT / ".github" / "workflows" / "ou3-proof.yml").read_text()
        self.assertEqual(2, text.count('"tools/stability/ou3_*.json"'))

    def test_unrelated_ou3_studies_stay_outside_stability_package(self):
        self.assertEqual(
            {
                "ou3_engine_noise_mitigation.py",
                "ou3_lever_arm_study.py",
                "ou3_lever_arm_tex.py",
            },
            {p.name for p in TOOLS.glob("ou3_*.py")},
        )


if __name__ == "__main__":
    unittest.main()
