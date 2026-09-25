"""Shared-versus-specific configuration contract of the three orchestrators.

src/kalman_common/SeaStateFusionDefaults.h defines, once, every default that
is the same physical or procedural quantity for OU-II, OU-III and TFG, and
tabulates the defaults that differ on purpose.  This test keeps both halves
honest:

* no estimator layer re-defines a shared default as a local literal, which is
  how the three copies drifted apart before they were merged; and
* the table of intentional differences states what the wrappers actually
  ship, so a difference cannot change -- or quietly disappear -- without the
  table being updated in the same change.
"""

from __future__ import annotations

import re
import unittest

from wrapper_sources import ROOT, estimator_layer, resolved_value, wrapper_source

DEFAULTS = (ROOT / "src/kalman_common/SeaStateFusionDefaults.h").read_text(encoding="utf-8")

# Shared defaults and the local names the wrappers expose them under.
SHARED_ALIASES = {
    "ACC_NOISE_FLOOR_SIGMA": ("ACC_NOISE_FLOOR_SIGMA_DEFAULT", "noise_floor_sigma_"),
    "TUNE_FREQ_PRIOR_HZ": ("TUNE_FREQ_PRIOR_HZ", "kTuneFreqPriorHz"),
    "MIN_TUNE_FREQ_HZ": ("MIN_TUNE_FREQ_HZ", "kMinTuneFreqHz"),
    "SIGMA_BAND_LOW_RATIO": ("SIGMA_BAND_LOW_RATIO_DEFAULT",),
    "SIGMA_BAND_HIGH_RATIO": ("SIGMA_BAND_HIGH_RATIO_DEFAULT",),
    "SIGMA_BAND_MIN_HZ": ("SIGMA_BAND_MIN_HZ_DEFAULT",),
    "SIGMA_BAND_MAX_HZ": ("SIGMA_BAND_MAX_HZ_DEFAULT",),
    "ADAPT_TAU_SEC": ("ADAPT_TAU_SEC", "ADAPT_TAU_SEC_DEFAULT"),
    "ADAPT_TAU_SEA_PERIODS": ("ADAPT_TAU_SEA_PERIODS", "ADAPT_TAU_SEA_PERIODS_DEFAULT"),
    "ADAPT_EVERY_SECS": ("ADAPT_EVERY_SECS", "adapt_every_secs_"),
    "MAG_DELAY_SEC": ("MAG_DELAY_SEC",),
    "STARTUP_PROXY_TWO_KP": ("STARTUP_PROXY_TWO_KP_DEFAULT", "proxy_two_kp"),
    "STARTUP_PROXY_TWO_KI": ("STARTUP_PROXY_TWO_KI_DEFAULT", "proxy_two_ki"),
    "ACC_VIBRATION_GUARD_HZ": ("ACC_VIBRATION_GUARD_HZ_DEFAULT",),
    "ACC_VIBRATION_GUARD_POLES": ("ACC_VIBRATION_GUARD_POLES_DEFAULT",),
    "ACC_VIBRATION_RACC_GAIN": ("ACC_VIBRATION_RACC_GAIN_DEFAULT",),
    "PSEUDO_UPDATE_PERIOD_NOMINAL_S": ("PSEUDO_UPDATE_PERIOD_NOMINAL_S", "kPseudoPeriodNominalS"),
    "PSEUDO_UPDATE_TAU_NOMINAL_S": ("PSEUDO_UPDATE_TAU_NOMINAL_S", "kPseudoTauNominalS"),
    "PSEUDO_UPDATE_PERIOD_MIN_S": ("PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT", "kPseudoPeriodMinS"),
    "NOMINAL_IMU_DT_S": ("FREQ_SMOOTHER_DT", "TFG_NOMINAL_DT"),
    "PROXY_STARTUP_TIMEOUT_SEC": ("proxy_startup_timeout_sec",),
    "GRAVITY_GATE_LPF_SEC": ("mag_gravity_align_world_tau_sec", "proxy_gravity_lpf_sec"),
    "MAG_HI_SLEW_TAU_SEC": ("mag_hi_slew_tau_sec",),
}

LITERAL = r"-?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?f?(?:\s*/\s*[0-9.]+f?)?"


def table() -> dict[str, tuple[str, str, str]]:
    """The 'Intentionally different' table, row label -> (OU-II, OU-III, TFG)."""
    block = DEFAULTS.split("Intentionally different", 1)[1]
    rows = {}
    for line in block.splitlines():
        cells = [c for c in re.split(r"\s{2,}", line.lstrip("/ ").strip()) if c]
        if len(cells) == 4 and cells[0] not in ("",) and not cells[0].startswith("OU-"):
            rows[cells[0]] = (cells[1], cells[2], cells[3])
    return rows


def number(cell: str) -> float:
    match = re.search(r"-?[0-9]+(?:\.[0-9]+)?", cell)
    assert match, cell
    value = float(match.group(0))
    if "ms" in cell:
        value /= 1000.0
    return value


class SharedDefaultsTest(unittest.TestCase):
    def test_shared_defaults_are_defined_once(self):
        for family in ("ou2", "ou3", "tfg"):
            layer = estimator_layer(family)
            for shared, aliases in SHARED_ALIASES.items():
                for alias in aliases:
                    with self.subTest(family=family, name=alias):
                        literal = re.search(
                            rf"\b{re.escape(alias)}\s*=\s*{LITERAL}\s*;", layer)
                        self.assertIsNone(
                            literal,
                            f"{family} re-defines shared default {shared} as a literal")

    def test_every_family_resolves_to_the_shared_value(self):
        for family in ("ou2", "ou3", "tfg"):
            source = wrapper_source(family)
            for shared, aliases in SHARED_ALIASES.items():
                expected = resolved_value(DEFAULTS, shared)
                self.assertIsNotNone(expected, shared)
                for alias in aliases:
                    value = resolved_value(source, alias)
                    if value is None:
                        continue  # this family has no such knob
                    with self.subTest(family=family, name=alias):
                        self.assertEqual(value, expected)


class IntentionalDifferencesTest(unittest.TestCase):
    """The documented table matches what each wrapper ships."""

    ROWS = {
        "max tuning frequency": ("MAX_TUNE_FREQ_HZ", "MAX_TUNE_FREQ_HZ", "kMaxTuneFreqHz"),
        "max sigma_aw": ("MAX_SIGMA_A", "MAX_SIGMA_A", "max_sigma_a_"),
        "max pseudo-update period": ("PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
                                     "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
                                     "kPseudoPeriodMaxS"),
        "tau_coeff": ("tau_coeff_", "tau_coeff_", "tau_coeff_"),
        "sigma_coeff": ("sigma_coeff_", "sigma_coeff_", "sigma_coeff_"),
        "mag refine start": ("mag_refine_start_sec", "mag_refine_start_sec",
                             "mag_refine_start_sec"),
    }

    def test_table_matches_sources(self):
        rows = table()
        for label, names in self.ROWS.items():
            self.assertIn(label, rows)
            for family, name, cell in zip(("ou2", "ou3", "tfg"), names, rows[label]):
                with self.subTest(row=label, family=family):
                    self.assertEqual(resolved_value(wrapper_source(family), name),
                                     number(cell))

    def test_documented_differences_are_real(self):
        rows = table()
        for label in self.ROWS:
            values = {number(cell) for cell in rows[label]}
            with self.subTest(row=label):
                self.assertGreater(len(values), 1,
                                   f"{label} is no longer different; move it to the shared defaults")

    def test_still_water_decay_and_laws(self):
        ou2, ou3 = estimator_layer("ou2"), estimator_layer("ou3")
        self.assertIn("float sigma_stillness_decay_sec_ = 5.0f;", ou2)
        self.assertIn("constexpr float STILL_VAR_DECAY_SEC = 1.0f;", ou3)
        self.assertIn("PseudoAdaptationLaw pseudo_law_ = PseudoAdaptationLaw::PhysicalMSE;", ou2)
        self.assertIn("RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;", ou3)
        self.assertIn("RSLaw rs_law_ = RSLaw::SpectralMSE;", estimator_layer("tfg"))
        self.assertIn("ADAPT_R_p0_MULT            = 3.0f", ou2)
        self.assertIn("ADAPT_RS_MULT              = 1.5f", ou3)
        self.assertIn("ADAPT_RS_MULT_DEFAULT          = 1.5f", estimator_layer("tfg"))
        self.assertIn("float P_factor_       = 1.5f;", ou2)
        self.assertIn("float S_factor_      = 1.0f;", ou3)


if __name__ == "__main__":
    unittest.main()
