from __future__ import annotations

import unittest
from pathlib import Path

from wrapper_sources import estimator_layer, wrapper_source

ROOT = Path(__file__).resolve().parents[2]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


class WavePeriodLiveEntryParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.period = text("src/tuner/WavePeriodEstimator.h")
        cls.ou2 = wrapper_source("ou2")
        cls.ou3 = wrapper_source("ou3")
        cls.tfg = wrapper_source("tfg")
        cls.ou2_layer = estimator_layer("ou2")
        cls.ou3_layer = estimator_layer("ou3")

    def test_single_estimator_has_staged_startup_qualification(self) -> None:
        for token in (
            "const float moment_start_sec = 3.0f / lambda_;",
            "const float usable_floor_sec = 4.0f / lambda_;",
            "const float moment_history_sec = elapsed_sec_ - moment_start_sec;",
            "moment_history_sec >= period",
            "const float settled_floor_sec = 6.0f / lambda_;",
            "bool hasUsablePeriod() const { return usable_period_; }",
            "if (usable_period_) return;",
            "usable_period_ = true;",
        ):
            self.assertIn(token, self.period)
        self.assertEqual(self.period.count("class WavePeriodEstimator"), 1)

    def test_ou_ii_and_ou_iii_require_measured_period_for_tuner_ready(self) -> None:
        # The shared stage machine promotes to TunerReady only with a usable
        # measured period, and both OU layers hand it exactly that.
        for source in (self.ou2, self.ou3):
            self.assertIn(
                "if (tuner_.isReady() && wave_period_usable) {",
                source,
            )
            self.assertIn(
                "{ return wave_period_.hasUsablePeriod(); }",
                source,
            )
        for layer in (self.ou2_layer, self.ou3_layer):
            self.assertIn("advanceStartupStage_(wavePeriodUsable())", layer)
        for source in (self.ou2, self.ou3):
            self.assertIn(
                "if (wave_period_.hasUsablePeriod() &&",
                source,
            )
            self.assertIn(
                "inline bool wavePeriodUsable() const noexcept",
                source,
            )

    def test_tfg_uses_same_startup_usable_contract(self) -> None:
        self.assertIn(
            "wave_period_.hasUsablePeriod()) {\n                enterLive_();",
            self.tfg,
        )
        self.assertIn(
            "return wave_period_.hasUsablePeriod() &&",
            self.tfg,
        )
        self.assertIn(
            "wave_period_.hasUsablePeriod() && std::isfinite(f)",
            self.tfg,
        )
        self.assertIn(
            "bool  wavePeriodUsable() const noexcept",
            self.tfg,
        )

    def test_prior_is_only_used_before_startup_usable_period(self) -> None:
        for source in (self.ou2, self.ou3, self.tfg):
            self.assertIn("hasUsablePeriod()", source)
        # Strict readiness remains a separate diagnostic; no startup gate should
        # require it in the three deployed orchestrators.
        self.assertNotIn(
            "tuner_.isReady() && wave_period_.isReady()",
            self.ou2 + self.ou3,
        )
        self.assertNotIn(
            "online_tune_warmup_sec && wave_period_.isReady()",
            self.tfg,
        )


if __name__ == "__main__":
    unittest.main()
