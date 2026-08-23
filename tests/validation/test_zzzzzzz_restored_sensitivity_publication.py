"""Late publication contract for the restored OU-III sensitivity study."""
from __future__ import annotations

import re
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_publication_robustness_sync as sensitivity_sync  # noqa: E402
import test_ou_legacy_law_cleanup as legacy_cleanup  # noqa: E402
import test_ou_robustness as robustness_core  # noqa: E402
import test_zzzz_editable_publication_contract as editable  # noqa: E402


def _ofat_rows(table: str) -> dict[float, tuple[float, ...]]:
    rows: dict[float, tuple[float, ...]] = {}
    body = table.split(r"\midrule", 1)[1].rsplit(r"\bottomrule", 1)[0]
    for line in body.splitlines():
        if "&" not in line or r"\\" not in line:
            continue
        cells = [cell.strip() for cell in line.split("&")]
        try:
            scale = float(cells[0])
        except ValueError:
            continue
        values: list[float] = []
        for cell in cells[1:4]:
            match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", cell)
            if match is None:
                raise AssertionError(f"no numeric mean in OFAT cell: {cell}")
            values.append(float(match.group(0)))
        rows[scale] = tuple(values)
    return rows


def _restored_manuscript_copies_match_generated_evidence(self):
    archived_publication = (
        self.RESULTS / "ou_robustness_publication.tex"
    ).read_text(encoding="utf-8")
    degradation_doc = (
        self.DOC / "w3d-ou-robustness-results-generated.tex-part"
    ).read_text(encoding="utf-8")

    archived_stress = editable._table_block_by_label(
        archived_publication, "tab:ou_robustness_stress"
    )
    doc_stress = editable._table_block_by_label(
        degradation_doc, "tab:ou_robustness_stress"
    )
    self.assertEqual(
        editable._publication_evidence_signature(archived_stress),
        editable._publication_evidence_signature(doc_stress),
        "degradation table payload differs from committed robustness evidence",
    )

    archived_sensitivity = editable._table_block_by_label(
        archived_publication, "tab:ou_robustness_sensitivity"
    )
    ofat_doc = (self.DOC / "w3d-ou-robustness-sensitivity-ofat-generated.tex-part").read_text(
        encoding="utf-8"
    )
    ofat_table = editable._table_block_by_label(
        ofat_doc, "tab:ou_robustness_sensitivity"
    )
    self.assertEqual(_ofat_rows(archived_sensitivity), _ofat_rows(ofat_table))

    archived_macros = editable._macro_definitions(
        (self.RESULTS / "ou_robustness_macros.tex").read_text(encoding="utf-8")
    )
    doc_macros = editable._macro_definitions(
        (self.DOC / "w3d-ou-robustness-macros-generated.tex-part").read_text(
            encoding="utf-8"
        )
    )
    self.assertTrue(doc_macros)
    for name, value in doc_macros.items():
        self.assertEqual(archived_macros[name], value, name)

    self.assertEqual(
        (self.RESULTS / "ou_robustness_stress.svg").read_bytes(),
        (self.DOC / "ou_robustness_stress.svg").read_bytes(),
        "ou_robustness_stress.svg",
    )

    source = (self.DOC / "w3d-ou-robustness.tex-part").read_text(encoding="utf-8")
    for required in (
        "Sensitivity and degradation cases",
        "one-factor-at-a-time",
        "coupled sweeps",
        "deployed SpectralMSE manifold",
        r"c^{6/7}",
        r"c^{24/7}",
        r"c^{41/14}",
        "w3d-ou-robustness-sensitivity-ofat-generated.tex-part",
        "w3d-ou-robustness-sensitivity-current-generated.tex-part",
    ):
        self.assertIn(required, source)
    for retired in (
        r"0.35\,\sigma_{aw}\tau^3",
        r"r_S\to c^{3}r_S",
        "legacy cubic",
    ):
        self.assertNotIn(retired, source)

    current = self.DOC / sensitivity_sync.CURRENT_TABLE
    svg = self.DOC / sensitivity_sync.SENSITIVITY_SVG
    if current.exists():
        current_text = current.read_text(encoding="utf-8")
        for marker in sensitivity_sync.CURRENT_MARKERS:
            self.assertIn(marker, current_text)
        for marker in sensitivity_sync.RETIRED_MARKERS:
            self.assertNotIn(marker, current_text)
        self.assertTrue(svg.exists())
    else:
        self.assertFalse(svg.exists())


def _publication_tree_keeps_studies_but_not_retired_coupling(self):
    source = (DOC / "w3d-ou-robustness.tex-part").read_text(encoding="utf-8")
    ofat = (DOC / "w3d-ou-robustness-sensitivity-ofat-generated.tex-part").read_text(
        encoding="utf-8"
    )
    self.assertIn("tab:ou_robustness_sensitivity", ofat)
    self.assertIn(r"c^{6/7}", source)
    self.assertIn(r"c^{41/14}", source)
    self.assertIn("SpectralMSE", source)
    for token in (
        r"0.35\,\sigma_{aw}\tau^3",
        r"\tau,r_S\propto\tau^3",
        "legacy cubic",
        "historical coupled",
    ):
        self.assertNotIn(token, source + "\n" + ofat)


robustness_core.CommittedRobustnessResultsTests.test_manuscript_copies_match_generated_evidence = (
    _restored_manuscript_copies_match_generated_evidence
)
legacy_cleanup.OULegacyLawCleanupTests.test_publication_tree_has_no_legacy_numerical_law_studies = (
    _publication_tree_keeps_studies_but_not_retired_coupling
)


class AnalyticalSensitivityPublicationSyncTests(unittest.TestCase):
    def test_legacy_bundle_cannot_publish_coupled_table_or_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            robustness = root / "robustness"
            doc = root / "doc"
            robustness.mkdir()
            doc.mkdir()
            (robustness / "ou_robustness_publication.tex").write_text(
                r"\begin{table*}\label{tab:ou_robustness_sensitivity}old\end{table*}",
                encoding="utf-8",
            )
            (doc / sensitivity_sync.CURRENT_TABLE).write_text("stale", encoding="utf-8")
            (doc / sensitivity_sync.SENSITIVITY_SVG).write_text("stale", encoding="utf-8")
            self.assertTrue(sensitivity_sync.sync_current_sensitivity(robustness, doc))
            self.assertFalse((doc / sensitivity_sync.CURRENT_TABLE).exists())
            self.assertFalse((doc / sensitivity_sync.SENSITIVITY_SVG).exists())

    def test_analytical_bundle_publishes_full_table_and_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            robustness = root / "robustness"
            doc = root / "doc"
            robustness.mkdir()
            doc.mkdir()
            table = (
                r"\begin{table*}\caption{deployed SpectralMSE law; $c^{6/7}$; $c^{41/14}$}"
                r"\label{tab:ou_robustness_sensitivity}current\end{table*}"
            )
            (robustness / "ou_robustness_publication.tex").write_text(table, encoding="utf-8")
            (robustness / sensitivity_sync.SENSITIVITY_SVG).write_text("svg", encoding="utf-8")
            self.assertTrue(sensitivity_sync.sync_current_sensitivity(robustness, doc))
            current = (doc / sensitivity_sync.CURRENT_TABLE).read_text(encoding="utf-8")
            self.assertIn("deployed SpectralMSE law", current)
            self.assertIn(r"c^{6/7}", current)
            self.assertIn(r"c^{41/14}", current)
            self.assertEqual((doc / sensitivity_sync.SENSITIVITY_SVG).read_text(), "svg")


if __name__ == "__main__":
    unittest.main()
