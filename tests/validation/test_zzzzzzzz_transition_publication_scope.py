"""Late publication contract for restored OU-III degradation evidence."""
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


def _publication_keeps_current_stress_evidence(self):
    """Publish degradation cases without reviving the retired coupled law."""
    archive_publication = (
        self.RESULTS / "ou_robustness_publication.tex"
    ).read_text(encoding="utf-8")
    self.assertIn("Rapid 30 s ramp", archive_publication)
    self.assertIn("Controlled 120 s ramp", archive_publication)

    source = (self.DOC / "w3d-ou-robustness.tex-part").read_text(encoding="utf-8")
    normalized = re.sub(r"\s+", " ", source).lower()
    for required in (
        "spectralmse",
        "one-factor-at-a-time",
        "c^{6/7}",
        "c^{41/14}",
        "low--high--low",
        "ou_robustness_stress.svg",
        "rapid one-way",
        "controlled",
        "degradation result",
    ):
        self.assertIn(required, normalized)
    for forbidden in (
        r"0.35\,\sigma_{aw}\tau^3",
        r"r_s\to c^{3}r_s",
        "legacy cubic",
    ):
        self.assertNotIn(forbidden, normalized)

    baseline = (self.DOC / "w3d-baseline-comparison.tex-part").read_text(
        encoding="utf-8"
    )
    self.assertIn(r"\input{w3d-ou-robustness.tex-part}", baseline)
    self.assertIn("ou_validation_vertical.svg", baseline)
    self.assertIn("ou_validation_displacement.svg", baseline)
    self.assertIn("ou_validation_attitude.svg", baseline)
    self.assertIn("SpectralMSE", baseline)

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

    self.assertEqual(
        (self.RESULTS / "ou_robustness_stress.svg").read_bytes(),
        (self.DOC / "ou_robustness_stress.svg").read_bytes(),
        "ou_robustness_stress.svg",
    )


def _publication_primary_claims_include_degradation(self):
    """Low-motion and rapid-transition claims must match paired evidence."""
    effects = self.read_csv(self.RESULTS / "ou_robustness_paired_effects.csv")
    generated = (
        self.DOC / "w3d-ou-robustness-macros-generated.tex-part"
    ).read_text(encoding="utf-8")
    checks = (
        ("Hs0.05_minus_Hs0.27", "disp_z_pct_hs"),
        ("rapid_minus_controlled_Adaptive", "disp_z_pct_hs"),
        ("Adaptive_minus_FixedNominal_rapid", "disp_z_pct_hs"),
    )
    for comparison, metric in checks:
        row = next(
            item for item in effects
            if item["comparison"] == comparison and item["metric"] == metric
        )
        for field in (
            "mean_paired_difference",
            "bootstrap_ci95_low",
            "bootstrap_ci95_high",
        ):
            self.assertIn(f"{float(row[field]):+.2f}", generated)


robustness_core.CommittedRobustnessResultsTests.test_manuscript_copies_match_generated_evidence = (
    _publication_keeps_current_stress_evidence
)
robustness_core.CommittedRobustnessResultsTests.test_generated_primary_claims_match_paired_effects = (
    _publication_primary_claims_include_degradation
)
