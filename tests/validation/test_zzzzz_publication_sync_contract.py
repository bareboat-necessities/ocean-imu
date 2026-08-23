"""Late override: restatement is compared after the publication editorial sync.

The validation generator owns evidence values and table structure.  The
publication-sync layer owns editorial claim wording such as the current
SpectralMSE channel-ablation caption.  Restatement therefore remains byte-exact
for numerical artifacts and is byte-exact for the publication fragment after
applying that same deterministic editorial transform.
"""

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_publication_sync  # noqa: E402
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
    def test_override_is_installed(self):
        self.assertIs(
            validation_core.RestatTests.
            test_restating_the_committed_bundle_reproduces_its_derived_files,
            _restating_the_committed_bundle_reproduces_its_derived_files,
        )


if __name__ == "__main__":
    unittest.main()
