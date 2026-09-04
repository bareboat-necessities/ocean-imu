#!/usr/bin/env python3
"""Contracts for the flattened OU-III proof-tool and workflow trees."""
from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
WORKFLOWS = ROOT / ".github" / "workflows"
VALIDATION = ROOT / "tests" / "validation"
OU3_DOC = ROOT / "doc" / "kalman_ou_iii"


class Ou3ProofCleanupTest(unittest.TestCase):
    def test_no_retired_p5_aliases_or_imports_remain(self):
        offenders = []
        for path in sorted(TOOLS.glob("ou3_*.py")):
            text = path.read_text(encoding="utf-8")
            if path.name.startswith("ou3_p5_"):
                offenders.append(f"retired filename: {path.name}")
            if "ou3_p5_" in text:
                offenders.append(f"retired reference: {path.name}")
        self.assertEqual([], offenders, "\n".join(offenders))

    def test_no_superseded_p4_fixed_schedule_routes_remain(self):
        retired_names = {
            "ou3_p4_terminal_cluster_p2_reduction.py",
            "ou3_p4_terminal_source_equivalence.py",
            "ou3_p4_joint_word_dissipation_design.py",
            "ou3_p4_joint_word_postprediction_design.py",
            "ou3_p4_joint_word_gauge_design.py",
            "ou3_p4_joint_word_gauge_design_v2.py",
            "ou3_p4_augmented_complete_word_design.py",
            "ou3_p4_augmented_complete_word_design_v5.py",
        }
        retired_modules = {name[:-3] for name in retired_names}
        offenders = []
        for name in sorted(retired_names):
            if (TOOLS / name).exists():
                offenders.append(f"retired filename: {name}")
        for path in sorted(TOOLS.glob("ou3_*.py")):
            text = path.read_text(encoding="utf-8")
            for module in sorted(retired_modules):
                if module in text:
                    offenders.append(f"retired reference {module}: {path.name}")
        self.assertEqual([], offenders, "\n".join(offenders))

    def test_no_superseded_p3_route_selection_diagnostics_remain(self):
        retired_names = {
            "ou3_p3_four_max_global_label_witness.py",
            "ou3_p3_ordered_witness_covariance_diagnostic.py",
        }
        retired_modules = {name[:-3] for name in retired_names}
        offenders = []
        for name in sorted(retired_names):
            if (TOOLS / name).exists():
                offenders.append(f"retired filename: {name}")
        for path in sorted(TOOLS.glob("ou3_*.py")):
            text = path.read_text(encoding="utf-8")
            for module in sorted(retired_modules):
                if module in text:
                    offenders.append(f"retired reference {module}: {path.name}")
        self.assertEqual([], offenders, "\n".join(offenders))

    def test_no_superseded_sea3_sizing_scaffolding_remains(self):
        retired = {
            TOOLS / "ou3_sea3_initial_estimates.py",
            TOOLS / "ou3_sea3_initial_estimates.json",
            VALIDATION / "test_ou3_sea3_initial_estimates.py",
        }
        offenders = [str(path.relative_to(ROOT)) for path in sorted(retired) if path.exists()]
        self.assertEqual([], offenders, "retired SEA3 sizing scaffolding returned")

    def test_no_retired_extra_proof_workflows_remain(self):
        retired = {
            "ou3-p3-sync-combined-caps-once.yml",
            "ou3-p3-whole-word-probe.yml",
            "ou3-p4-new-clamps-proof.yml",
            "ou3-sea0.yml",
        }
        offenders = [name for name in sorted(retired) if (WORKFLOWS / name).exists()]
        self.assertEqual([], offenders, "retired extra proof workflows returned")

    def test_master_paper_includes_retained_sea3_theorem_parts(self):
        paper = (OU3_DOC / "kalman_ou-w3d.tex").read_text(encoding="utf-8")
        self.assertIn(r"\input{w3d-sea3-stability-theorem.tex-part}", paper)
        self.assertIn(r"\input{w3d-sea3-period-bridge.tex-part}", paper)
        self.assertLess(
            paper.index(r"\input{w3d-sea3-stability-theorem.tex-part}"),
            paper.index(r"\input{w3d-sea3-period-bridge.tex-part}"),
        )


if __name__ == "__main__":
    unittest.main()
