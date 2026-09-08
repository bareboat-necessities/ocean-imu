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
            start = "30" if source == self.tfg else "90"
            self.assertRegex(source, re.compile(rf"mag_refine_start_sec\s*=\s*{start}\.0f"))
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

    def test_current_ablation_tables_are_generated_from_full_replay(self):
        import csv
        import importlib.util
        producer = REPO_ROOT / "tools/startup_publication_sync.py"
        spec = importlib.util.spec_from_file_location("startup_publication_sync", producer)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with (REPO_ROOT / "reports/results/startup_ablation/startup_runs.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 192)
        hard, study = mod.generate(rows)
        self.assertEqual(hard, (DOC / "w3d-hard-iron-results-generated.tex-part").read_text())
        self.assertEqual(study, (DOC / "w3d-startup-results-generated.tex-part").read_text())
        self.assertIn(r"\input{w3d-startup-results-generated.tex-part}", self.paper)
        self.assertIn("misalignment-dominated", self.paper)
        self.assertIn("soft-iron/misalignment error without additional heading excitation", self.flat)

    def test_ieee_two_column_tables_stay_single_column(self):
        self.assertNotIn(r"\begin{table*}", self.paper)
        self.assertNotIn(r"\begin{figure*}", self.paper)
        self.assertNotIn(r"\begin{table*}", (DOC / "w3d-startup-results-generated.tex-part").read_text())


if __name__ == "__main__":
    unittest.main()
