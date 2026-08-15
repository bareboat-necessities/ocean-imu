"""Publication contract for figures and integrated validation in kalman-wave-dir.tex."""

import csv
import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PLOTS = REPO_ROOT / "plots" / "kalman_ou_iii"
VALIDATION_SUMMARY = REPO_ROOT / "reports" / "results" / "ou_validation" / "ou_validation_summary.csv"


class WaveDirectionChartContractTests(unittest.TestCase):
    def test_direction_article_uses_only_direction_pgfs(self):
        article = (DOC / "kalman-wave-dir.tex").read_text(encoding="utf-8")
        results = (DOC / "w3d-wave-direction-results.tex-part").read_text(
            encoding="utf-8"
        )
        frequency = (DOC / "w3d-frequency-tracking.tex-part").read_text(
            encoding="utf-8"
        )
        config = (DOC / "w3d-wave-direction-config.tex-part").read_text(
            encoding="utf-8"
        )
        charts = (DOC / "w3d-wave-direction-charts.tex-part").read_text(
            encoding="utf-8"
        )
        plotter = (PLOTS / "kalman_ou_iii-plots.py").read_text(encoding="utf-8")

        self.assertIn(r"\usepackage{graphicx}", article)
        self.assertIn(r"\usepackage{pgf}", article)
        self.assertIn(r"\usepackage{pgfplots}", article)
        self.assertIn(r"\providecommand{\mathdefault}[1]{#1}", article)
        self.assertIn(r"\input{w3d-wave-direction-charts.tex-part}", results)

        # The shared plot pipeline still creates both direction and OU tuning
        # diagnostics, but this publication consumes only the direction family.
        self.assertIn('finalize_plot(fig, outbase, "_dir")', plotter)
        for name in (
            "w3d_ou3_jonswap_medium_dir.pgf",
            "w3d_ou3_jonswap_low_dir.pgf",
            "w3d_ou3_pmstokes_medium_dir.pgf",
        ):
            self.assertIn(name, charts)

        publication_sources = article + results + frequency + config + charts
        self.assertNotIn("w3d_ou3_jonswap_medium_tuner.pgf", publication_sources)
        for token in ("R_S", "r_S", r"R_{S}", r"r_{S}", "tau_applied", "sigma_a_applied", "sigma_{aw}"):
            self.assertNotIn(token, publication_sources)

        # OU--III is allowed only where the paper explicitly reports the
        # integrated replay that supplied its upstream attitude input.
        non_integration_sources = article + frequency + config + charts
        self.assertNotIn("OU--III", non_integration_sources)
        self.assertIn(r"\subsection{Integrated replay with OU--III attitude input}", results)

        # IEEE conference output is two-column.  Publication plots stay inside
        # one column; no figure* or text-width scaling may creep back in.
        self.assertNotIn(r"\begin{figure*}", charts)
        self.assertNotIn(r"\textwidth", charts)
        self.assertEqual(charts.count(r"\resizebox{\columnwidth}{!}"), 3)

    def test_integrated_ou3_table_is_scoped_and_mirrors_committed_results(self):
        results = (DOC / "w3d-wave-direction-results.tex-part").read_text(
            encoding="utf-8"
        )
        source = (DOC / "w3d-sim-results-generated.tex-part").read_text(
            encoding="utf-8"
        )

        self.assertIn(r"\label{tab:direction-ou3-integration}", results)
        normalized_results = " ".join(results.split())
        self.assertIn(
            "integration evidence, not an OU--III performance claim",
            normalized_results,
        )
        self.assertIn("OUValidationDirectionRMSE", results)
        self.assertNotIn("Z RMS", results)
        self.assertNotIn(r"Z\%$H_s$", results)

        row_pattern = re.compile(
            r"^\s*(JONSWAP|PM--Stokes)\s*&\s*([0-9.]+)\s*&\s*"
            r"[0-9.]+\s*&\s*[0-9.]+\s*&\s*"
            r"([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*"
            r"(-?[0-9.]+)\s*&\s*([0-9.]+/[0-9.]+/[0-9.]+)\s*\\\\$",
            re.MULTILINE,
        )
        rows = row_pattern.findall(source)
        self.assertEqual(len(rows), 8)

        # Copy only the upstream attitude RMS and direction fields.  If the
        # committed deterministic evidence changes, this table must be restated
        # rather than silently drifting out of sync.
        for case, hs, roll, pitch, yaw, theta, tau in rows:
            expected = re.compile(
                rf"{re.escape(case)}\s*&\s*{re.escape(hs)}\s*&\s*"
                rf"{re.escape(roll)}\s*&\s*{re.escape(pitch)}\s*&\s*"
                rf"{re.escape(yaw)}\s*&\s*{re.escape(theta)}\s*&\s*"
                rf"{re.escape(tau)}"
            )
            self.assertRegex(results, expected)

    def test_ou3_direction_rms_table_mirrors_ten_seed_summary(self):
        results = (DOC / "w3d-wave-direction-results.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\label{tab:direction-ou3-rms}", results)

        scenario_labels = {
            "stationary_jonswap_H0_270_L14_047_A30_00_P60_00": ("JONSWAP", "0.27"),
            "stationary_jonswap_H1_500_L50_710_A_30_00_P120_00": ("JONSWAP", "1.50"),
            "stationary_jonswap_H4_000_L112_766_A30_00_P30_00": ("JONSWAP", "4.00"),
            "stationary_jonswap_H8_500_L202_839_A_30_00_P72_00": ("JONSWAP", "8.50"),
            "stationary_pmstokes_H0_270_L14_047_A30_00_P60_00": ("PM--Stokes", "0.27"),
            "stationary_pmstokes_H1_500_L50_710_A_30_00_P120_00": ("PM--Stokes", "1.50"),
            "stationary_pmstokes_H4_000_L112_766_A30_00_P30_00": ("PM--Stokes", "4.00"),
            "stationary_pmstokes_H8_500_L202_839_A_30_00_P72_00": ("PM--Stokes", "8.50"),
        }
        metrics = {}
        with VALIDATION_SUMMARY.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                scenario = row["scenario"]
                if (
                    scenario in scenario_labels
                    and row["family"] == "OU_III"
                    and row["mode"] == "Adaptive"
                    and row["metric"] in {"dir_axis_rmse_deg", "dir_travel_rmse_deg"}
                ):
                    metrics[(scenario, row["metric"])] = (
                        float(row["mean"]),
                        float(row["std"]),
                    )

        self.assertEqual(len(metrics), 16)
        for scenario, (case, hs) in scenario_labels.items():
            axis_mean, axis_std = metrics[(scenario, "dir_axis_rmse_deg")]
            travel_mean, travel_std = metrics[(scenario, "dir_travel_rmse_deg")]
            expected = (
                rf"{re.escape(case)}\s*&\s*{re.escape(hs)}\s*&\s*"
                rf"\${axis_mean:.2f}\\pm{axis_std:.2f}\$\s*&\s*"
                rf"\${travel_mean:.2f}\\pm{travel_std:.2f}\$"
            )
            self.assertRegex(results, expected)


if __name__ == "__main__":
    unittest.main()
