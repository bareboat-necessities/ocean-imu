import csv
import hashlib
import json
import math
import sys
import unittest
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_robustness as robustness  # noqa: E402
import ou_validation as validation  # noqa: E402


class RobustnessDesignTests(unittest.TestCase):
    @staticmethod
    def metric_row(**values):
        row = {
            "experiment": "sensitivity",
            "case": "tau_x1",
            "parameter": "tau",
            "scale_label": "1",
            "scale_multiplier": 1.0,
            "mode": "FixedSensitivity",
            "wave_phase_seed": 1,
            "imu_noise_seed": 2,
            "initialization_seed": 3,
        }
        row.update({metric: 1.0 for metric in validation.METRIC_NAMES})
        row.update(values)
        return row

    def test_sensitivity_scales_exactly_one_parameter(self):
        baseline = validation.TuningPoint(
            tau_s=1.2,
            sigma_a_mps2=0.8,
            RS_ms=3.5,
        )
        tau = robustness.scaled_tuning_point(baseline, "tau", 2.0)
        sigma = robustness.scaled_tuning_point(baseline, "sigma_aw", 0.5)
        r_s = robustness.scaled_tuning_point(baseline, "r_s", 1.5)
        self.assertEqual(tau, validation.TuningPoint(2.4, 0.8, RS_ms=3.5))
        self.assertEqual(sigma, validation.TuningPoint(1.2, 0.4, RS_ms=3.5))
        self.assertEqual(r_s, validation.TuningPoint(1.2, 0.8, RS_ms=5.25))
        robustness.validate_tuning_point(r_s)
        # Just past the implementation ceiling, whatever it currently is.
        with self.assertRaisesRegex(ValueError, "outside implementation bounds"):
            robustness.validate_tuning_point(
                validation.TuningPoint(
                    robustness.TAU_BOUNDS_S[1] * 1.01, 0.8, RS_ms=3.5
                )
            )

    def test_coupled_sweeps_follow_the_deployed_regularization_law(self):
        # r_S = clip(c sigma_aw tau^3, ...), so the coupled sweeps must move
        # r_S linearly with sigma_aw and cubically with tau. A frozen r_S
        # measures only the direct OU process-covariance effect, which the
        # online tuner never produces.
        baseline = validation.TuningPoint(
            tau_s=1.2,
            sigma_a_mps2=0.8,
            RS_ms=3.5,
        )
        coupled_sigma = robustness.scaled_tuning_point(
            baseline, "sigma_aw_rs", 0.5
        )
        coupled_tau = robustness.scaled_tuning_point(baseline, "tau_rs", 1.5)
        self.assertEqual(
            coupled_sigma, validation.TuningPoint(1.2, 0.4, RS_ms=1.75)
        )
        self.assertAlmostEqual(coupled_tau.tau_s, 1.8)
        self.assertAlmostEqual(coupled_tau.sigma_a_mps2, 0.8)
        self.assertAlmostEqual(coupled_tau.RS_ms, 3.5 * 1.5**3)
        self.assertEqual(
            robustness.SENSITIVITY_PARAMETERS,
            robustness.OFAT_PARAMETERS + robustness.COUPLED_PARAMETERS,
        )

    def test_scale_parser_requires_reference(self):
        self.assertEqual(robustness.parse_float_list("0.5,1,2"), [0.5, 1.0, 2.0])
        with self.assertRaisesRegex(Exception, "include the 1.0 reference"):
            robustness.parse_float_list("0.5,2")

    def test_ci_smoke_scales_cover_publication_anchors(self):
        workflow = (
            REPO_ROOT / ".github" / "workflows" / "ou-validation.yml"
        ).read_text(encoding="utf-8")
        # The pull-request gate runs smoke mode and must still exercise the
        # three anchors the publication tables are read at; the dispatch-only
        # regeneration leg runs the full five.
        self.assertIn('--sensitivity-scales "0.5,1.0,1.5"', workflow)
        self.assertIn('--sensitivity-scales "0.5,0.75,1.0,1.25,1.5"', workflow)

    def test_summary_and_sensitivity_effect_are_paired(self):
        rows = []
        for repetition, reference, doubled in ((1, 8.0, 9.0), (2, 10.0, 13.0)):
            seed_fields = {
                "wave_phase_seed": repetition,
                "imu_noise_seed": 100 + repetition,
                "initialization_seed": 200 + repetition,
            }
            rows.append(
                self.metric_row(
                    **seed_fields,
                    disp_z_pct_hs=reference,
                )
            )
            rows.append(
                self.metric_row(
                    **seed_fields,
                    case="tau_x2",
                    scale_label="2",
                    scale_multiplier=2.0,
                    disp_z_pct_hs=doubled,
                )
            )

        summary = robustness.summarize_rows(rows, 500, 7)
        selected = next(
            row for row in summary
            if row["case"] == "tau_x2" and row["metric"] == "disp_z_pct_hs"
        )
        self.assertEqual(selected["n"], 2)
        self.assertAlmostEqual(selected["mean"], 11.0)

        effects = robustness.paired_effect_rows(rows, 500, 7)
        selected_effect = next(
            row for row in effects
            if row["comparison"] == "tau_x2_minus_x1"
            and row["metric"] == "disp_z_pct_hs"
        )
        self.assertEqual(selected_effect["n_pairs"], 2)
        self.assertAlmostEqual(selected_effect["mean_paired_difference"], 2.0)
        self.assertTrue(math.isfinite(selected_effect["cohen_dz"]))


class CommittedRobustnessResultsTests(unittest.TestCase):
    RESULTS = REPO_ROOT / "reports" / "results" / "ou_robustness"
    DOC = REPO_ROOT / "doc" / "kalman_ou_iii"

    @staticmethod
    def read_csv(path):
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    def test_full_bundle_is_paired_bounded_and_self_hashed(self):
        with (self.RESULTS / "ou_robustness_manifest.json").open(
            encoding="utf-8"
        ) as stream:
            manifest = json.load(stream)
        protocol = manifest["protocol"]
        self.assertEqual(protocol["mode"], "full")
        self.assertEqual(protocol["score_window_sec"], 900.0)
        self.assertEqual(len(protocol["seed_triplets"]), 10)
        self.assertEqual(
            protocol["sensitivity_scales"], [0.5, 0.75, 1.0, 1.25, 1.5]
        )
        self.assertEqual(
            [case["end_sec"] - case["start_sec"] for case in protocol["transition_cases"]],
            [360.0, 30.0],
        )
        self.assertIn("analytically derived", protocol["wave_phase_method"])
        self.assertIn(
            "first- and second-derivative", protocol["transition_method"]
        )

        raw = self.read_csv(self.RESULTS / "ou_robustness_raw.csv")
        summary = self.read_csv(self.RESULTS / "ou_robustness_summary.csv")
        effects = self.read_csv(
            self.RESULTS / "ou_robustness_paired_effects.csv"
        )
        scales = protocol["sensitivity_scales"]
        parameters = robustness.SENSITIVITY_PARAMETERS
        self.assertEqual(
            protocol["sensitivity_parameters"], list(parameters)
        )
        sensitivity_cells = len(parameters) * len(scales)
        low_motion_cells = 2
        transition_cells = len(protocol["transition_cases"]) * 2
        cells = sensitivity_cells + low_motion_cells + transition_cells
        self.assertEqual(len(raw), cells * 10)

        # Checked against the metrics the committed bundle carries rather than
        # against the tool's current list, so adding a metric asks for a re-run
        # instead of failing this structural check.  Unknown metrics still fail.
        bundle_metrics = {row["metric"] for row in summary}
        self.assertTrue(
            bundle_metrics <= set(validation.NON_SEGMENT_METRIC_NAMES),
            sorted(bundle_metrics - set(validation.NON_SEGMENT_METRIC_NAMES)),
        )
        self.assertEqual(len(summary), cells * len(bundle_metrics))
        # Sensitivity: every off-reference scale against x1, per direction.
        # Degradation: one low-motion pair, plus rate and adaptation pairs.
        comparisons = len(parameters) * (len(scales) - 1) + 1 + 4
        self.assertEqual(len(effects), comparisons * len(bundle_metrics))
        self.assertEqual(
            Counter(row["experiment"] for row in raw),
            {
                "sensitivity": sensitivity_cells * 10,
                "low_motion": low_motion_cells * 10,
                "transition_rate": transition_cells * 10,
            },
        )
        groups = Counter(
            (row["experiment"], row["case"], row["mode"]) for row in raw
        )
        self.assertEqual(len(groups), cells)
        self.assertEqual(set(groups.values()), {10})
        self.assertEqual({int(row["n"]) for row in summary}, {10})
        self.assertEqual({int(row["n_pairs"]) for row in effects}, {10})
        self.assertEqual({int(row["samples"]) for row in raw}, {180000})
        self.assertEqual({float(row["window_s"]) for row in raw}, {900.0})
        self.assertTrue(
            all(
                int(row["simulator_return_code"])
                == 1 - int(row["historical_60s_gate_pass"])
                for row in raw
            )
        )

        for row in raw:
            if row["experiment"] != "sensitivity":
                continue
            # Bounds come from the module rather than from literals, so that
            # moving the implementation clamps cannot leave this test asserting
            # a range the filter no longer has.
            tau_lo, tau_hi = robustness.TAU_BOUNDS_S
            sigma_lo, sigma_hi = robustness.SIGMA_AW_BOUNDS_MPS2
            rs_lo, rs_hi = robustness.R_S_BOUNDS_MS
            self.assertGreaterEqual(float(row["configured_tau_s"]), tau_lo)
            self.assertLessEqual(float(row["configured_tau_s"]), tau_hi)
            self.assertGreater(float(row["configured_sigma_aw_mps2"]), sigma_lo)
            self.assertLessEqual(float(row["configured_sigma_aw_mps2"]), sigma_hi)
            self.assertGreaterEqual(float(row["configured_r_s_ms"]), rs_lo)
            self.assertLessEqual(float(row["configured_r_s_ms"]), rs_hi)

        for name, metadata in manifest["result_files"].items():
            path = self.RESULTS / name
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), metadata["sha256"], name
            )
            self.assertEqual(path.stat().st_size, metadata["bytes"], name)
        for name, metadata in manifest["source_files"].items():
            path = REPO_ROOT / name
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), metadata["sha256"], name
            )
            self.assertEqual(path.stat().st_size, metadata["bytes"], name)

    def test_manuscript_copies_match_generated_evidence(self):
        for result_name, doc_name in (
            (
                "ou_robustness_publication.tex",
                "w3d-ou-robustness-results-generated.tex-part",
            ),
            ("ou_robustness_sensitivity.svg", "ou_robustness_sensitivity.svg"),
            ("ou_robustness_stress.svg", "ou_robustness_stress.svg"),
        ):
            self.assertEqual(
                (self.RESULTS / result_name).read_bytes(),
                (self.DOC / doc_name).read_bytes(),
                doc_name,
            )

        source = (self.DOC / "w3d-ou-robustness.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn("no sweep point is adopted as a replacement tuning", source)
        self.assertIn("w3d-ou-robustness-results-generated.tex-part", source)
        # The coupled directions and the reason the frozen-companion columns
        # look flat must both be stated, or the sweep invites the wrong
        # conclusion that the filter ignores its OU tuning.
        self.assertIn("coupled", source)
        self.assertIn(r"r_S\to c^{3}r_S", source)
        self.assertIn("near saturation", source)

    def test_generated_primary_claims_match_paired_effects(self):
        effects = self.read_csv(
            self.RESULTS / "ou_robustness_paired_effects.csv"
        )
        generated = (
            self.DOC / "w3d-ou-robustness-results-generated.tex-part"
        ).read_text(encoding="utf-8")
        checks = (
            ("r_s_x0.5_minus_x1", "disp_3d_rms_m", 3),
            ("sigma_aw_x0.5_minus_x1", "pitch_rms_deg", 3),
            ("Hs0.05_minus_Hs0.27", "disp_z_pct_hs", 2),
            ("rapid_minus_controlled_Adaptive", "disp_z_pct_hs", 2),
            ("Adaptive_minus_FixedNominal_rapid", "disp_z_pct_hs", 2),
        )
        for comparison, metric, digits in checks:
            row = next(
                item for item in effects
                if item["comparison"] == comparison and item["metric"] == metric
            )
            for field in (
                "mean_paired_difference",
                "bootstrap_ci95_low",
                "bootstrap_ci95_high",
            ):
                self.assertIn(f"{float(row[field]):+.{digits}f}", generated)


if __name__ == "__main__":
    unittest.main()
