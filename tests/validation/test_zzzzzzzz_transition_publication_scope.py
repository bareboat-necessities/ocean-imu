"""Late publication contract: one-way transition evidence is provenance-only."""
from pathlib import Path
import re

import test_ou_robustness as robustness_core
import test_zzzz_editable_publication_contract  # noqa: F401  apply earlier overrides first

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _macro_definitions(text):
    return dict(
        re.findall(
            r"\\providecommand\{\\(OURobustness[A-Za-z]+)\}\{([^}]*)\}", text
        )
    )


def _publication_excludes_one_way_transition(self):
    """Keep archived one-way evidence available, but never render it in the manuscript."""
    archive_publication = (
        self.RESULTS / "ou_robustness_publication.tex"
    ).read_text(encoding="utf-8")
    self.assertIn("Rapid 30 s ramp", archive_publication)
    self.assertIn("Controlled 120 s ramp", archive_publication)

    # Publication synchronization may mirror archive-derived support files and
    # macro definitions into the document tree for provenance/reproducibility.
    # Scope is determined by what the manuscript actually inputs/references,
    # not by whether an unused support file or macro definition exists.
    source = (self.DOC / "w3d-ou-robustness.tex-part").read_text(encoding="utf-8")
    normalized = re.sub(r"\s+", " ", source).lower()
    for required in (
        "spectralmse",
        "one-factor-at-a-time",
        "c^{6/7}",
        "c^{41/14}",
        "low--high--low",
        "provenance only",
    ):
        self.assertIn(required, normalized)
    for forbidden in (
        "rapid 30",
        "controlled 120",
        "rapid transition degradation",
        "w3d-ou-robustness-results-generated.tex-part",
        "ou_robustness_stress.svg",
        "\\ourobustnessrapid",
        "\\ourobustnesscontrolled",
    ):
        self.assertNotIn(forbidden, normalized)

    archived_macros = _macro_definitions(
        (self.RESULTS / "ou_robustness_macros.tex").read_text(encoding="utf-8")
    )
    doc_macros = _macro_definitions(
        (self.DOC / "w3d-ou-robustness-macros-generated.tex-part").read_text(
            encoding="utf-8"
        )
    )
    self.assertTrue(doc_macros)
    for name, value in doc_macros.items():
        self.assertIn(name, archived_macros)
        self.assertEqual(value, archived_macros[name], name)


def _publication_primary_claims_are_low_motion_only(self):
    """Only low-motion degradation macros may be consumed as primary robustness claims."""
    effects = self.read_csv(self.RESULTS / "ou_robustness_paired_effects.csv")
    generated = (
        self.DOC / "w3d-ou-robustness-macros-generated.tex-part"
    ).read_text(encoding="utf-8")
    source = (self.DOC / "w3d-ou-robustness.tex-part").read_text(encoding="utf-8")
    row = next(
        item for item in effects
        if item["comparison"] == "Hs0.05_minus_Hs0.27"
        and item["metric"] == "disp_z_pct_hs"
    )
    for field in (
        "mean_paired_difference",
        "bootstrap_ci95_low",
        "bootstrap_ci95_high",
    ):
        self.assertIn(f"{float(row[field]):+.2f}", generated)

    # Historical one-way macros may be defined by the synchronized archive,
    # but the rendered publication must not consume them.
    self.assertNotIn("\\OURobustnessRapid", source)
    self.assertNotIn("\\OURobustnessControlled", source)


robustness_core.CommittedRobustnessResultsTests.test_manuscript_copies_match_generated_evidence = (
    _publication_excludes_one_way_transition
)
robustness_core.CommittedRobustnessResultsTests.test_generated_primary_claims_match_paired_effects = (
    _publication_primary_claims_are_low_motion_only
)
