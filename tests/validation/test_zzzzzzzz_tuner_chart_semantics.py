"""Late publication contract for the OU-III tuner diagnostic chart."""
from pathlib import Path
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
DRAW = REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh"


class OUIIITunerChartSemanticTests(unittest.TestCase):
    def test_applied_equivalent_frequency_is_not_called_raw_tuner_input(self):
        text = DRAW.read_text(encoding="utf-8")
        self.assertIn("Applied-equivalent wave frequency", text)
        self.assertIn("1.0 / (2.0 * tau_for_plot)", text)
        self.assertIn("deliberately not named", text)
        old_panel = '("freq_tracker_hz", "Frequency (Hz)")'
        new_panel = '("wave_tuning_freq_hz", r"Applied-equivalent wave frequency (Hz)")'
        self.assertIn(old_panel, text)
        self.assertIn(new_panel, text)
        self.assertIn("source.replace(old_panel, new_panel, 1)", text)

    def test_acceleration_and_ou_sigma_are_distinguished(self):
        text = DRAW.read_text(encoding="utf-8")
        self.assertIn("Band accel std (pre noise-floor subtraction)", text)
        self.assertIn(r"$\\sigma_{aw}$ applied", text)

    def test_legacy_shared_column_is_presented_as_ou3_r_s(self):
        text = DRAW.read_text(encoding="utf-8")
        self.assertIn("R_p0_applied is a legacy shared-harness column name", text)

        # The old label is intentionally present as the source anchor for the
        # publication-only rewrite.  Verify the rewrite contract rather than
        # globally banning the anchor string from the wrapper.
        old_panel = '("p0_combo",        r"$R_{p0}$ / $p_{0,S}$ applied")'
        new_panel = '("p0_combo",        r"$r_S$ applied ($m\\\\,s$)")'
        self.assertIn(old_panel, text)
        self.assertIn(new_panel, text)
        self.assertIn(
            "source.replace(old_regularizer_panel, new_regularizer_panel, 1)",
            text,
        )

        # The legacy shared-harness data column is also relabeled in the
        # conditional plotting branch, so no generated OU-III chart calls it p0.
        self.assertIn(
            'label=r"$r_S$ applied")',
            text,
        )
        self.assertIn('label=r"legacy regularizer column")', text)


if __name__ == "__main__":
    unittest.main()
