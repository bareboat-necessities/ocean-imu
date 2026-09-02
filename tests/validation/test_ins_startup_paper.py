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
HI_DOC = REPO_ROOT / "docs" / "continuous-mag-hard-iron.md"


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
        cls.hi_doc = HI_DOC.read_text(encoding="utf-8")

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

    def test_tfg_keeps_the_matched_staged_ablation(self):
        # TFG still carries both startup paths so the published comparison can
        # be re-run.  The OU wrappers dropped the staged path once the proxy
        # became the only deployed one.
        self.assertIn("enum class StartupInitPolicy", self.tfg)
        self.assertIn("StagedMekf", self.tfg)
        self.assertIn("MahonyProxy", self.tfg)
        self.assertRegex(
            self.tfg,
            re.compile(r"startup_init_policy\s*=\s*StartupInitPolicy::MahonyProxy"),
        )

    def test_all_three_wrappers_keep_the_proxy_path(self):
        for source in (self.ou3, self.ou2, self.tfg):
            self.assertIn("proxy_handoff_tilt_sigma_rad", source)
            self.assertIn("proxy_handoff_yaw_sigma_free_rad", source)
            self.assertIn("ContinuousMagHardIronEstimator", source)

    def test_ou_wrappers_dropped_the_staged_startup(self):
        for source in (self.ou3, self.ou2):
            self.assertNotIn("StagedMekf", source)
            self.assertNotIn("enum class StartupInitPolicy", source)
            self.assertNotIn("set_linear_block_enabled", source)

    def test_two_stage_magnetic_policy_matches_sources(self):
        for source in (self.ou3, self.ou2, self.tfg):
            self.assertRegex(source, re.compile(r"mag_refine_enabled\s*=\s*true"))
            self.assertRegex(source, re.compile(r"mag_refine_start_sec\s*=\s*90\.0f"))
            self.assertRegex(source, re.compile(r"mag_refine_window_sec\s*=\s*30\.0f"))
            self.assertRegex(
                source, re.compile(r"mag_continuous_hard_iron\s*=\s*true")
            )

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

    def test_ou2_marginal_startup_attribution_matches_evidence(self):
        for row in (
            r"roll RMS [deg] & 0.348 & 0.286 & $-17.8\%$",
            r"pitch RMS [deg] & 0.308 & 0.300 & $-2.6\%$",
            r"yaw RMS [deg] & 2.360 & 2.411 & $+2.2\%$",
            r"accel-bias RMS [m/s$^2$] & 0.0616 & 0.0520 & $-15.6\%$",
        ):
            self.assertIn(row, self.paper)

    def test_hard_iron_results_match_committed_evidence(self):
        # The paper and the note both quote the ablation as re-measured on the
        # re-cut ridge floor.  The numbers the correction originally landed
        # with are a different baseline and are pinned separately below, so
        # the two cannot be confused for each other again.
        for row in (
            r"yaw RMS mean [deg] & 2.065 & 0.617 & $-70.1\%$",
            r"yaw RMS worst [deg] & 2.132 & 0.878 & $-58.8\%$",
            r"yaw RMS mean [deg] & 1.869 & 0.661 & $-64.6\%$",
            r"yaw RMS worst [deg] & 2.133 & 1.078 & $-49.5\%$",
        ):
            self.assertIn(row, self.paper)

        self.assertIn("2.065 deg mean / 2.133 deg worst", self.hi_doc)
        self.assertIn("0.617 deg mean /", self.hi_doc)
        self.assertIn("| **mean** | | 0.768 | **0.617** | 0.751 | **0.661** |",
                      self.hi_doc)
        self.assertIn("| **worst** | | 1.265 | **0.878** | 1.067 | 1.078 |",
                      self.hi_doc)

    def test_hard_iron_landing_baseline_is_kept_as_history(self):
        # The note keeps what the correction landed with, and says so.  This
        # pins that framing: the older numbers must stay, and must not be the
        # ones the paper quotes.
        self.assertIn("1.835 deg mean / 2.162 deg worst", self.hi_doc)
        self.assertIn("**1.887** | **0.813**", self.hi_doc)
        self.assertIn("**2.161** | **1.089**", self.hi_doc)
        self.assertNotIn(r"yaw RMS mean [deg] & 1.835", self.paper)

    def test_tfg_feature_attribution_matches_evidence(self):
        for row in (
            "Deployed & 0.777 & 1.345 & 0.296 & 0.366 & 97.1 & 164.1",
            "Hard iron off & 2.486 & 3.256 & 0.307 & 0.373 & 94.3 & 159.3",
            "Mag refinement off & 0.688 & 1.534 & 0.600 & 0.241 & 118.0 & 204.3",
            "Staged startup & 2.265 & 2.906 & 0.697 & 0.303 & 145.7 & 398.3",
        ):
            self.assertIn(row, self.paper)

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
