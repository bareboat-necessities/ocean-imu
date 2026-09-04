from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


class WavePeriodLiveEntryParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.period = text("src/tuner/WavePeriodEstimator.h")
        cls.ou2 = text("src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h")
        cls.ou3 = text("src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h")
        cls.tfg = text("src/kalman_tfg/SeaStateFusionFilter_TFG.h")

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
        for source in (self.ou2, self.ou3):
            self.assertIn(
                "if (tuner_.isReady() && wave_period_.hasUsablePeriod()) {",
                source,
            )
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
