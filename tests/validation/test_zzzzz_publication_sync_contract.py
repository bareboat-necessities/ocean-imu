"""Late override: restatement is compared after the publication editorial sync.

The validation generator owns evidence values and table structure.  The
publication-sync layer owns editorial claim wording such as the current
SpectralMSE channel-ablation caption.  Restatement therefore remains byte-exact
for numerical artifacts and is byte-exact for the publication fragment after
applying that same deterministic editorial transform.
"""

import contextlib
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_publication_sync  # noqa: E402
import ou_publication_robustness_sync as robustness_sync  # noqa: E402
import test_ou_validation as validation_core  # noqa: E402


def _restating_the_committed_bundle_reproduces_its_derived_files(self):
    source = self.RESULTS / "ou_validation.json"
    if not source.exists():  # pragma: no cover - smoke checkouts
        self.skipTest("no committed validation bundle")

    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp)
        with contextlib.redirect_stdout(io.StringIO()):
            validation_core.validation.restat_bundle(
                source,
                output,
                bootstrap_resamples=10_000,
                stats_seed=20260317,
            )

        # The committed publication fragment has already passed through
        # publication-sync before unittest discovery. Apply the same deterministic
        # caption transform to the restated copy before comparing it.
        publication = output / ou_publication_sync.PUBLICATION_NAME
        text = publication.read_text(encoding="utf-8")
        text, _ = ou_publication_sync._replace_caption_before_label(
            text,
            "tab:ou_mc_channels",
            ou_publication_sync.CURRENT_CHANNEL_CAPTION,
        )
        publication.write_text(text, encoding="utf-8")

        for name in (
            "ou_validation_raw.csv",
            "ou_validation_summary.csv",
            "ou_validation_paired_effects.csv",
            "ou_validation_table.tex",
            "ou_validation_publication.tex",
            "ou_validation_macros.tex",
            "ou_validation_tuning_points.tex",
        ):
            self.assertEqual(
                (output / name).read_bytes(),
                (self.RESULTS / name).read_bytes(),
                name,
            )

        restated = json.loads(
            (output / "ou_validation.json").read_text(encoding="utf-8")
        )
        committed = json.loads(source.read_text(encoding="utf-8"))
        self.assertEqual(restated["raw_runs"], committed["raw_runs"])
        self.assertIn("restated_from", restated["protocol"])


validation_core.RestatTests.test_restating_the_committed_bundle_reproduces_its_derived_files = (
    _restating_the_committed_bundle_reproduces_its_derived_files
)


class PublicationSyncRestatOverrideTests(unittest.TestCase):
    def test_ofat_excerpt_follows_the_new_study(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            results, doc = root / "results", root / "doc"
            results.mkdir()
            doc.mkdir()
            (results / "ou_robustness_publication.tex").write_text(
                r"\begin{table*}\label{tab:ou_robustness_sensitivity}"
                r"\midrule 0.5 & 7.02 & 6.96 & 7.39 & 7.27 & 13.92 \\"
                r"\bottomrule\end{table*}")
            target = doc / robustness_sync.OFAT_TABLE
            target.write_text(r"\midrule old rows\bottomrule")
            self.assertTrue(robustness_sync.sync_ofat_excerpt(results, doc))
            self.assertIn(r"0.5 & 7.02 & 6.96 & 7.39 \\", target.read_text())
            self.assertNotIn("13.92", target.read_text())
            self.assertFalse(robustness_sync.sync_ofat_excerpt(results, doc))

    def test_empty_direction_summary_displays_no_angle(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            with (directory / ou_publication_sync.SUMMARY_NAME).open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("scenario", "family", "mode", "metric", "mean", "std", "n"))
                for scenario, _, _ in ou_publication_sync.DIRECTION_SCENARIOS:
                    for metric in ("dir_axis_rmse_deg", "dir_travel_rmse_deg"):
                        writer.writerow((scenario, "OU_III", "Adaptive", metric, "nan", "nan", 0))
            rows = ou_publication_sync._direction_summary_rows(directory)
            self.assertEqual(sum(" & -- & -- & 0 " in row for row in rows), 8)
            self.assertNotIn("nan", "\n".join(rows))

    def test_unresolved_deterministic_direction_is_preserved(self):
        source = REPO_ROOT / "doc/kalman_ou_iii/w3d-sim-results-generated.tex-part"
        rows = ou_publication_sync._deterministic_direction_rows(source)
        self.assertEqual(sum(" & " in row for row in rows), 8)
        self.assertEqual(sum("{--}" in row for row in rows), 2)
        self.assertEqual(sum("0.0/0.0/100.0" in row for row in rows), 2)

    def test_override_is_installed(self):
        self.assertIs(
            validation_core.RestatTests.
            test_restating_the_committed_bundle_reproduces_its_derived_files,
            _restating_the_committed_bundle_reproduces_its_derived_files,
        )


if __name__ == "__main__":
    unittest.main()
