"""Publication/source contract for ins-startup.tex."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PAPER = DOC / "ins-startup.tex"
OU3 = REPO_ROOT / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
OU2 = REPO_ROOT / "src" / "kalman_ou_ii" / "SeaStateFusionFilter_OU_II.h"
TFG = REPO_ROOT / "src" / "kalman_tfg" / "SeaStateFusionFilter_TFG.h"
HI = REPO_ROOT / "src" / "tuner" / "ContinuousMagHardIronEstimator.h"
STARTUP_DOC = REPO_ROOT / "docs" / "ou-iii-startup-init.md"
HI_DOC = REPO_ROOT / "docs" / "continuous-mag-hard-iron.md"
OU2_DOC = REPO_ROOT / "docs" / "ou-ii-ou-iii-parity.md"
TFG_DOC = REPO_ROOT / "docs" / "tfg-adaptation-parity.md"


def compact(text: str) -> str:
    return " ".join(text.split())


class InsStartupPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paper = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.paper)
        cls.ou3 = OU3.read_text(encoding="utf-8")
        cls.ou2 = OU2.read_text(encoding="utf-8")
        cls.tfg = TFG.read_text(encoding="utf-8")
        cls.hi = HI.read_text(encoding="utf-8")
        cls.startup_doc = STARTUP_DOC.read_text(encoding="utf-8")
        cls.hi_doc = HI_DOC.read_text(encoding="utf-8")
        cls.ou2_doc = OU2_DOC.read_text(encoding="utf-8")
        cls.tfg_doc = TFG_DOC.read_text(encoding="utf-8")

    def test_standalone_title_and_scope(self):
        self.assertIn(
            "Marine Attitude Initialization and Continuous Magnetic Self-Calibration",
            self.flat,
        )
        for phrase in (
            "MahonyProxy",
            "StagedMekf",
            "ContinuousMagHardIronEstimator",
            "two-stage magnetic acquisition",
            "gauge error",
        ):
            self.assertIn(phrase, self.paper)
        self.assertIn(r"\bibliography{w3d}", self.paper)

    def test_all_three_wrappers_keep_proxy_and_staged_paths(self):
        for source in (self.ou3, self.ou2, self.tfg):
            self.assertIn("enum class StartupInitPolicy", source)
            self.assertIn("StagedMekf", source)
            self.assertIn("MahonyProxy", source)
            self.assertRegex(
                source,
                re.compile(
                    r"startup_init_policy\s*=\s*StartupInitPolicy::MahonyProxy"
                ),
            )
            self.assertIn("ContinuousMagHardIronEstimator", source)

    def test_two_stage_magnetic_policy_matches_sources(self):
        for source in (self.ou3, self.ou2, self.tfg):
            self.assertRegex(source, re.compile(r"mag_refine_enabled\s*=\s*true"))
            self.assertRegex(source, re.compile(r"mag_refine_start_sec\s*=\s*90\.0f"))
            self.assertRegex(source, re.compile(r"mag_refine_window_sec\s*=\s*30\.0f"))
            self.assertRegex(
                source, re.compile(r"mag_continuous_hard_iron\s*=\s*true")
            )

        self.assertIn("90 s, 30 s window", self.ou2_doc)
        self.assertIn("90 s", self.startup_doc)

    def test_proxy_integral_term_and_handoff_numbers_are_current(self):
        # The OU wrappers expose their shared proxy gains through named constants;
        # TFG carries the same values directly in Config.
        for source in (self.ou3, self.ou2):
            self.assertRegex(
                source,
                re.compile(r"STARTUP_PROXY_TWO_KP_DEFAULT\s*=\s*0\.2f"),
            )
            self.assertRegex(
                source,
                re.compile(r"STARTUP_PROXY_TWO_KI_DEFAULT\s*=\s*0\.02f"),
            )
        self.assertRegex(self.tfg, re.compile(r"proxy_two_kp\s*=\s*0\.2f"))
        self.assertRegex(self.tfg, re.compile(r"proxy_two_ki\s*=\s*0\.02f"))

        for source in (self.ou3, self.ou2, self.tfg):
            self.assertRegex(
                source, re.compile(r"proxy_startup_timeout_sec\s*=\s*150\.0f")
            )

        self.assertIn("0.711 deg", self.flat)
        self.assertIn("0.05 deg/s", self.flat)
        self.assertIn("0.71 deg", self.startup_doc)

    def test_hard_iron_derivation_matches_implementation(self):
        for token in (
            "Eigen::Matrix3d::Identity() - A.transpose() * A",
            "M.trace() / 3.0",
            "cfg_.model_ridge_relative",
            "weight_sum_ * lambda_min",
            "max_residual_rms_uT",
            "bias_body_uT",
        ):
            self.assertIn(token, self.hi)

        for token in (
            r"\left(\mat I-\bar{\mat A}\T\bar{\mat A}\right)\vct b",
            r"\mat M=\mat I-\bar{\mat A}\T\bar{\mat A}",
            r"\lambda_{\rm rel}\frac{\operatorname{tr}(\mat M)}{3}",
            r"N_{\rm eff}\lambda_{\min}(\mat M)",
        ):
            self.assertIn(token, self.paper)

    def test_ou3_startup_table_matches_committed_evidence(self):
        expected_rows = (
            r"Staged MEKF & $\sim22$ s & 4.301 & 19.872 & 0.372 & 0.275 & 1.796 & 0.0631",
            r"Mahony proxy, 2-stage & 22--52 s & 4.287 & 19.851 & 0.339 & 0.266 & 1.835 & 0.0574",
            r"Mahony proxy, late single-stage & $\sim105$ s & -- & -- & 0.310 & 0.259 & 1.839 & 0.0517",
        )
        for row in expected_rows:
            self.assertIn(row, self.paper)

        for evidence in (
            "0.372 | 0.275 | 1.796 | 0.0631",
            "0.339 | 0.266 | 1.835 | 0.0574",
            "0.310 | 0.259 | 1.839 | 0.0517",
        ):
            self.assertIn(evidence, self.startup_doc)

    def test_ou2_marginal_startup_attribution_matches_evidence(self):
        for row in (
            r"roll RMS [deg] & 0.348 & 0.286 & $-17.8\%$",
            r"pitch RMS [deg] & 0.308 & 0.300 & $-2.6\%$",
            r"yaw RMS [deg] & 2.360 & 2.411 & $+2.2\%$",
            r"accel-bias RMS [m/s$^2$] & 0.0616 & 0.0520 & $-15.6\%$",
        ):
            self.assertIn(row, self.paper)

        self.assertIn(
            "| roll RMS, deg | 0.354 | 0.348 | 0.348 | 0.286 |", self.ou2_doc
        )
        self.assertIn(
            "| accel-bias 3D RMS, m/s^2 | 0.0595 | 0.0616 | 0.0616 | 0.0520 |",
            self.ou2_doc,
        )
        self.assertIn("0.444 to 0.414 deg", self.ou2_doc)

    def test_hard_iron_results_match_committed_evidence(self):
        for row in (
            r"yaw RMS mean [deg] & 1.835 & 0.822 & $-55.2\%$",
            r"yaw RMS worst [deg] & 2.162 & 1.063 & $-50.8\%$",
            r"yaw RMS mean [deg] & 1.887 & 0.813 & $-56.9\%$",
            r"yaw RMS worst [deg] & 2.161 & 1.089 & $-49.6\%$",
        ):
            self.assertIn(row, self.paper)

        self.assertIn("1.835 deg mean / 2.162 deg worst", self.hi_doc)
        self.assertIn("0.822 deg mean /", self.hi_doc)
        self.assertIn("**1.887** | **0.813**", self.hi_doc)
        self.assertIn("**2.161** | **1.089**", self.hi_doc)

    def test_tfg_feature_attribution_matches_evidence(self):
        for row in (
            "Deployed & 1.077 & 1.528 & 0.308 & 0.340 & 92.7 & 166.7",
            "Hard iron off & 2.349 & 3.303 & 0.316 & 0.350 & 90.2 & 162.2",
            "Mag refinement off & 0.873 & 1.264 & 0.614 & 0.259 & 120.7 & 204.1",
            "Staged startup & 2.202 & 2.906 & 0.696 & 0.308 & 144.8 & 398.1",
        ):
            self.assertIn(row, self.paper)

        for evidence in (
            "| **deployed** | 4.369 | 4.778 | 19.74 | 21.03 | 1.077 | 1.528 | 0.308 | 0.340 | 92.7 | 166.7 |",
            "| `TFG_MAG_HARD_IRON=0` | 4.403 | 4.775 | 19.73 | 21.01 | **2.349** | **3.303** | 0.316 | 0.350 | 90.2 | 162.2 |",
            "| `TFG_STARTUP_INIT=staged` | 4.427 | 4.809 | 20.44 | **24.27** | 2.202 | 2.906 | 0.696 | 0.308 | 144.8 | **398.1** |",
        ):
            self.assertIn(evidence, self.tfg_doc)

    def test_limit_on_soft_iron_is_not_hidden(self):
        for phrase in (
            "5 of 40 pairs became worse",
            "misalignment-dominated",
            "2.98 to 2.74 deg",
            "1.02 deg",
            "cannot identify soft-iron/misalignment error without additional heading excitation",
        ):
            self.assertIn(phrase, self.flat)

        self.assertIn("5 / 40", self.hi_doc)
        self.assertIn("seed=23", self.hi_doc)

    def test_ieee_two_column_tables_stay_single_column(self):
        self.assertNotIn(r"\begin{table*}", self.paper)
        self.assertNotIn(r"\begin{figure*}", self.paper)
        self.assertGreaterEqual(self.paper.count(r"\resizebox{\columnwidth}{!}"), 2)


if __name__ == "__main__":
    unittest.main()
