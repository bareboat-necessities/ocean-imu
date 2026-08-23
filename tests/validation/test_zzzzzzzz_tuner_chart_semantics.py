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
        self.assertNotIn('new_panel = \'        ("wave_tuning_freq_hz", "Applied wave-band frequency (Hz)"),\'', text)

    def test_acceleration_and_ou_sigma_are_distinguished(self):
        text = DRAW.read_text(encoding="utf-8")
        self.assertIn("Band accel std (pre noise-floor subtraction)", text)
        self.assertIn(r"$\\sigma_{aw}$ applied", text)

    def test_legacy_shared_column_is_presented_as_ou3_r_s(self):
        text = DRAW.read_text(encoding="utf-8")
        self.assertIn(r"$r_S$ applied", text)
        self.assertIn("R_p0_applied is a legacy shared-harness column name", text)
        self.assertNotIn(r'new_regularizer_panel = \'        ("p0_combo",        r"$R_{p0}$', text)


if __name__ == "__main__":
    unittest.main()
