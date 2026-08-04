import csv
import hashlib
import json
import math
import os
import re
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_sim_table as sim_table  # noqa: E402
import ou_validation as validation  # noqa: E402


class WaveSurrogateTests(unittest.TestCase):
    def setUp(self):
        self.columns = [
            "time",
            "disp_x", "disp_y", "disp_z",
            "vel_x", "vel_y", "vel_z",
            "acc_x", "acc_y", "acc_z",
            "acc_bx", "acc_by", "acc_bz",
            "gyro_x", "gyro_y", "gyro_z",
            "roll_deg", "pitch_deg", "yaw_deg",
            "q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z",
        ]
        count = 4096
        time = np.arange(count) * validation.DT_SECONDS
        self.data = np.zeros((count, len(self.columns)), dtype=np.float64)
        self.data[:, 0] = time
        frequencies = (0.25, 0.45, 0.75)
        for offset, name in enumerate(("disp_x", "disp_y", "disp_z")):
            angular_frequency = 2.0 * np.pi * frequencies[offset]
            angle = angular_frequency * time + 0.2 * offset
            self.data[:, self.columns.index(name)] = np.sin(angle)
            axis = "xyz"[offset]
            self.data[:, self.columns.index(f"vel_{axis}")] = (
                angular_frequency * np.cos(angle)
            )
            self.data[:, self.columns.index(f"acc_{axis}")] = (
                -angular_frequency**2 * np.sin(angle)
            )
        self.data[:, self.columns.index("roll_deg")] = 3.0 * np.sin(0.7 * time)
        self.data[:, self.columns.index("pitch_deg")] = 2.0 * np.cos(0.9 * time)
        self.data[:, self.columns.index("yaw_deg")] = 5.0
        self.data = validation.rebuild_body_imu(self.columns, self.data)

    def test_common_phase_preserves_primitive_spectra_and_kinematics(self):
        randomized = validation.phase_randomize_wave(self.columns, self.data, seed=73)
        names = (
            "vel_x", "vel_y", "vel_z",
            "roll_deg", "pitch_deg", "yaw_deg",
        )
        indices = validation._column_indices(self.columns, names)
        original_fft = np.fft.rfft(self.data[:, indices], axis=0)
        randomized_fft = np.fft.rfft(randomized[:, indices], axis=0)
        frequencies = np.fft.rfftfreq(len(self.data), validation.DT_SECONDS)
        retained = (
            (frequencies >= validation.SURROGATE_MIN_FREQ_HZ)
            & (frequencies <= validation.SURROGATE_MAX_FREQ_HZ)
        )
        np.testing.assert_allclose(
            np.abs(randomized_fft[retained]),
            np.abs(original_fft[retained]),
            rtol=2e-11,
            atol=2e-9,
        )
        np.testing.assert_allclose(
            randomized_fft[retained, 0] * np.conj(randomized_fft[retained, 1]),
            original_fft[retained, 0] * np.conj(original_fft[retained, 1]),
            rtol=2e-11,
            atol=2e-8,
        )

        omega = 2.0 * np.pi * frequencies
        for axis in "xyz":
            displacement = np.fft.rfft(
                randomized[:, self.columns.index(f"disp_{axis}")]
            )
            velocity = np.fft.rfft(
                randomized[:, self.columns.index(f"vel_{axis}")]
            )
            acceleration = np.fft.rfft(
                randomized[:, self.columns.index(f"acc_{axis}")]
            )
            np.testing.assert_allclose(
                displacement[retained],
                velocity[retained] / (1j * omega[retained]),
                rtol=2e-11,
                atol=2e-8,
            )
            np.testing.assert_allclose(
                acceleration[retained],
                1j * omega[retained] * velocity[retained],
                rtol=2e-11,
                atol=2e-8,
            )

        quaternion_indices = validation._column_indices(
            self.columns,
            ("q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z"),
        )
        norms = np.linalg.norm(randomized[:, quaternion_indices], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=2e-12)

    def test_phase_seed_changes_realization(self):
        first = validation.phase_randomize_wave(self.columns, self.data, seed=1)
        second = validation.phase_randomize_wave(self.columns, self.data, seed=2)
        index = self.columns.index("disp_z")
        self.assertGreater(np.max(np.abs(first[:, index] - second[:, index])), 0.1)

    def test_smoothstep_has_exact_endpoints(self):
        values = validation.smoothstep_weight(np.array([0.0, 2.0, 3.0, 4.0, 8.0]), 2.0, 4.0)
        np.testing.assert_allclose(values, (0.0, 0.0, 0.5, 1.0, 1.0))

        _, first, second = validation.smoothstep_profile(
            np.array([0.0, 2.0, 3.0, 4.0, 8.0]), 2.0, 4.0
        )
        np.testing.assert_allclose(first[[0, 1, 3, 4]], 0.0)
        np.testing.assert_allclose(second[[0, 1, 3, 4]], 0.0)

    def test_nonstationary_blend_includes_kinematic_cross_terms(self):
        generated = validation.make_nonstationary_wave(
            self.columns,
            self.data,
            self.data,
            seed=19,
            end_scale=1.4,
            transition_start_sec=5.0,
            transition_end_sec=15.0,
        )
        interior = slice(4, -4)
        for axis in "xyz":
            displacement = generated[:, self.columns.index(f"disp_{axis}")]
            velocity = generated[:, self.columns.index(f"vel_{axis}")]
            acceleration = generated[:, self.columns.index(f"acc_{axis}")]
            numerical_velocity = np.gradient(
                displacement, validation.DT_SECONDS, edge_order=2
            )
            numerical_acceleration = np.gradient(
                velocity, validation.DT_SECONDS, edge_order=2
            )
            velocity_scale = np.sqrt(np.mean(velocity[interior] ** 2))
            acceleration_scale = np.sqrt(np.mean(acceleration[interior] ** 2))
            self.assertLess(
                np.sqrt(np.mean((numerical_velocity[interior] - velocity[interior]) ** 2))
                / velocity_scale,
                2e-3,
            )
            self.assertLess(
                np.sqrt(
                    np.mean(
                        (numerical_acceleration[interior] - acceleration[interior]) ** 2
                    )
                )
                / acceleration_scale,
                2e-3,
            )


class StatisticsTests(unittest.TestCase):
    @staticmethod
    def row(family, mode, repetition, value):
        row = {
            "scenario": "sea",
            "family": family,
            "mode": mode,
            "repetition": repetition,
            "wave_phase_seed": repetition,
            "imu_noise_seed": 10 + repetition,
            "initialization_seed": 20 + repetition,
        }
        row.update({metric: value for metric in validation.METRIC_NAMES})
        return row

    def test_seed_broadcasting_is_paired(self):
        seeds = validation.broadcast_seed_triplets([1, 2], [10], [20, 21])
        self.assertEqual(
            seeds,
            [
                validation.SeedTriplet(1, 10, 20),
                validation.SeedTriplet(2, 10, 21),
            ],
        )

    def test_summary_and_paired_effect_use_all_repetitions(self):
        rows = []
        for repetition, ou2, ou3 in ((1, 3.0, 2.0), (2, 5.0, 3.0), (3, 7.0, 4.0)):
            rows.append(self.row("OU_II", "Adaptive", repetition, ou2))
            rows.append(self.row("OU_III", "Adaptive", repetition, ou3))
        summary = validation.summarize_rows(rows, bootstrap_resamples=500, stats_seed=7)
        selected = next(
            row for row in summary
            if row["family"] == "OU_II"
            and row["mode"] == "Adaptive"
            and row["metric"] == "disp_3d_rms_m"
        )
        self.assertEqual(selected["n"], 3)
        self.assertAlmostEqual(selected["mean"], 5.0)
        self.assertAlmostEqual(selected["std"], 2.0)

        effects = validation.paired_effect_rows(rows, 500, 7)
        selected_effect = next(
            row for row in effects
            if row["comparison"] == "OU_III_minus_OU_II"
            and row["metric"] == "disp_3d_rms_m"
        )
        self.assertEqual(selected_effect["n_pairs"], 3)
        self.assertAlmostEqual(selected_effect["mean_paired_difference"], -2.0)
        self.assertTrue(math.isfinite(selected_effect["cohen_dz"]))
        self.assertTrue(math.isfinite(selected_effect["hedges_gz"]))

    def test_machine_readable_metric_parser(self):
        parsed = validation.parse_validation_metrics(
            "VALIDATION_METRICS family=OU_II tuning_mode=adaptive "
            "aw_cov_sync=periodic "
            "input=wave.csv window_s=900 samples=180000 disp_z_rms_m=0.25"
        )
        self.assertEqual(parsed["family"], "OU_II")
        self.assertEqual(parsed["aw_cov_sync"], "periodic")
        self.assertEqual(parsed["samples"], 180000)
        self.assertEqual(parsed["disp_z_rms_m"], 0.25)

    def test_metric_parser_folds_in_segments_and_axial_magnitude(self):
        parsed = validation.parse_validation_metrics(
            "VALIDATION_METRICS family=OU_III tuning_mode=adaptive "
            "aw_cov_sync=periodic input=wave.csv segment=window start_s=300 "
            "window_s=900 samples=180000 disp_z_rms_m=0.25 "
            "disp_z_pct_refrms=12.5 dir_axis_error_deg=-1.5\n"
            "VALIDATION_SEGMENT family=OU_III tuning_mode=adaptive "
            "aw_cov_sync=periodic input=wave.csv segment=blend start_s=420 "
            "window_s=360 samples=72000 disp_z_rms_m=0.31 "
            "disp_z_pct_refrms=17.25\n"
        )
        # The axis is defined mod 180 deg, so signed errors from records of
        # opposite heading cancel; the magnitude is what the claim needs.
        self.assertEqual(parsed["dir_axis_error_deg"], -1.5)
        self.assertEqual(parsed["dir_axis_abs_error_deg"], 1.5)
        self.assertEqual(parsed["seg_blend_disp_z_rms_m"], 0.31)
        self.assertEqual(parsed["seg_blend_disp_z_pct_refrms"], 17.25)
        # The trailing-window row keeps its own values, not the segment's.
        self.assertEqual(parsed["disp_z_pct_refrms"], 12.5)
        self.assertNotIn("segment", parsed)

    def test_generated_macro_names_contain_letters_only(self):
        # TeX ends a control-sequence name at the first non-letter, so a
        # generated name with a digit in it does not fail loudly: the tail of
        # the name is typeset as ordinary text next to the value.  This is what
        # put "27NominalDifference-13.118" into the rendered manuscript, so
        # generation must refuse to emit such a name.
        with self.assertRaises(ValueError) as caught:
            validation.assert_latex_macro_names(
                [r"\providecommand{\OUValidationIIHs027NominalDifference}{-13.1}"]
            )
        self.assertIn("letters only", str(caught.exception))
        validation.assert_latex_macro_names(
            [r"\providecommand{\OUValidationIIHsZeroTwoSevenNominalDifference}{-13.1}"]
        )

    def test_scenario_macro_tokens_are_legal_and_distinct(self):
        scenarios = (
            "stationary_jonswap_H0_270_L14_047_A30_00_P60_00",
            "stationary_jonswap_H1_500_L50_710_A_30_00_P120_00",
            "stationary_jonswap_H4_000_L112_766_A30_00_P30_00",
            "stationary_jonswap_H8_500_L202_839_A_30_00_P72_00",
            "stationary_pmstokes_H1_500_L50_710_A_30_00_P120_00",
            "nonstationary_H1_5_to_H4_0_Tp5_7_to_11_4",
        )
        tokens = [validation.scenario_macro_token(name) for name in scenarios]
        for token in tokens:
            self.assertTrue(token.isalpha(), token)
        self.assertEqual(len(set(tokens)), len(tokens))

    def test_partial_adaptation_modes_form_a_factorial(self):
        # AdaptiveRSOnly and AdaptiveOUOnly complete a 2x2 with Adaptive and
        # FixedNominal: each freezes exactly one channel, and all four share
        # the deployed covariance policy so the only factor that moves is
        # which parameters are estimated online.
        self.assertEqual(
            validation.MODE_SETTINGS["AdaptiveRSOnly"],
            ("adaptive_rs_only", "periodic"),
        )
        self.assertEqual(
            validation.MODE_SETTINGS["AdaptiveOUOnly"],
            ("adaptive_ou_only", "periodic"),
        )
        for mode in validation.PARTIAL_MODES:
            self.assertEqual(validation.MODE_SETTINGS[mode][1], "periodic")
            # OU-II regularizes with (R_p0, R_v0) rather than a single
            # integral scale, so the channel ablation is not defined for it.
            self.assertNotIn(mode, validation.FAMILY_MODES["OU_II"])
            self.assertIn(mode, validation.FAMILY_MODES["OU_III"])

    def test_pmstokes_is_not_pooled_into_the_primary_aggregate(self):
        jonswap = "stationary_jonswap_H1_500_L50_710_A_30_00_P120_00"
        pmstokes = "stationary_pmstokes_H1_500_L50_710_A_30_00_P120_00"
        self.assertEqual(validation.scenario_spectrum(jonswap), "jonswap")
        self.assertEqual(validation.scenario_spectrum(pmstokes), "pmstokes")

        rows = []
        for repetition in (1, 2, 3):
            for family, offset in (("OU_II", 0.0), ("OU_III", -1.0)):
                for scenario, value in ((jonswap, 8.0), (pmstokes, 40.0)):
                    row = self.row(family, "Adaptive", repetition, value)
                    row["scenario"] = scenario
                    row["disp_z_pct_hs"] = value + repetition + offset
                    rows.append(row)

        primary = validation.stationary_normalized_aggregate(rows, 200, 7)
        secondary = validation.stationary_normalized_aggregate(
            rows, 200, 7, spectrum="pmstokes"
        )
        # Pooling would drag the primary mean toward the PM-Stokes level.
        self.assertAlmostEqual(primary["OU_II"]["mean"], 10.0)
        self.assertAlmostEqual(secondary["OU_II"]["mean"], 42.0)

    def test_segment_metrics_are_optional(self):
        # Segment columns exist only for the non-stationary scenario, so the
        # summary must treat a missing key as an absence rather than fail.
        rows = [self.row("OU_III", "Adaptive", repetition, 1.0) for repetition in (1, 2)]
        for row in rows:
            for metric in list(row):
                if metric.startswith("seg_"):
                    del row[metric]
        summary = validation.summarize_rows(rows, 100, 3)
        segment = [
            row for row in summary
            if row["metric"] == "seg_blend_disp_z_pct_refrms"
        ]
        self.assertTrue(segment)
        self.assertEqual(segment[0]["n"], 0)
        self.assertNotIn(
            "seg_blend_disp_z_pct_refrms", validation.NON_SEGMENT_METRIC_NAMES
        )

    def test_covariance_sync_ablation_is_a_matched_pair(self):
        # Each *PeriodicSync* mode must differ from its partner in exactly one
        # factor: the covariance-reset policy. Otherwise the ablation cannot
        # separate that factor from online adaptation.
        for left, right in validation.COVARIANCE_SYNC_PAIRS:
            left_tuning, left_sync = validation.MODE_SETTINGS[left]
            right_tuning, right_sync = validation.MODE_SETTINGS[right]
            self.assertEqual(left_tuning, right_tuning)
            self.assertEqual(left_sync, "periodic")
            self.assertEqual(right_sync, "reconfigure")
        # All three primary modes must share the deployed covariance policy,
        # otherwise the adaptation ablation is confounded by it.
        for mode in validation.PRIMARY_MODES:
            self.assertEqual(validation.MODE_SETTINGS[mode][1], "periodic")

    def test_paired_effects_separate_the_two_families(self):
        rows = []
        for repetition, low, high in ((1, 1.0, 2.0), (2, 1.5, 3.0), (3, 2.0, 4.0)):
            rows.append(self.row("OU_II", "Adaptive", repetition, low))
            rows.append(self.row("OU_II", "FixedNominal", repetition, low - 0.25))
            rows.append(self.row("OU_III", "Adaptive", repetition, high))
            rows.append(self.row("OU_III", "FixedNominal", repetition, high - 1.0))
        effects = validation.paired_effect_rows(rows, 200, 7)
        by_family = {
            str(row["family"]): row
            for row in effects
            if row["comparison"] == "Adaptive_minus_FixedNominal"
            and row["metric"] == "disp_z_pct_hs"
        }
        self.assertAlmostEqual(
            by_family["OU_II"]["mean_paired_difference"], 0.25
        )
        self.assertAlmostEqual(
            by_family["OU_III"]["mean_paired_difference"], 1.0
        )

    def test_transition_window_composition_is_reported(self):
        composition = validation.transition_window_composition(
            transition_start_sec=420.0,
            transition_end_sec=780.0,
            duration_sec=1200.0,
            window_sec=900.0,
        )
        self.assertAlmostEqual(composition["window_start_sec"], 300.0)
        self.assertAlmostEqual(composition["pure_start_sea_sec"], 120.0)
        self.assertAlmostEqual(composition["blended_sec"], 360.0)
        self.assertAlmostEqual(composition["pure_end_sea_sec"], 420.0)

    def test_crossfade_midpoint_is_below_a_linear_height_ramp(self):
        # The two blended records are independent, so their variances add and
        # the crossfade never reaches the height a linear H_s ramp would.
        midpoint = validation.mixture_significant_height_m(1.5, 4.0, 0.5)
        self.assertAlmostEqual(midpoint, math.hypot(0.75, 2.0), places=6)
        self.assertLess(midpoint, 0.5 * (1.5 + 4.0))
        self.assertAlmostEqual(
            validation.mixture_significant_height_m(1.5, 4.0, 0.0), 1.5
        )
        self.assertAlmostEqual(
            validation.mixture_significant_height_m(1.5, 4.0, 1.0), 4.0
        )

    def test_publication_labels_and_stationary_aggregate(self):
        small = "stationary_jonswap_H0_270_L14_047_A30_00_P60_00"
        large = "stationary_jonswap_H4_000_L112_766_A30_00_P30_00"
        transition = "nonstationary_H1_5_to_H4_0_Tp5_7_to_11_4"
        self.assertEqual(validation.scenario_display_label(small), "$H_s=0.27$ m")
        self.assertEqual(validation.scenario_display_label(transition), "Transition")
        self.assertLess(
            validation.scenario_sort_key(small),
            validation.scenario_sort_key(large),
        )
        self.assertLess(
            validation.scenario_sort_key(large),
            validation.scenario_sort_key(transition),
        )

        rows = []
        for repetition in (1, 2, 3):
            for family, offset in (("OU_II", 0.0), ("OU_III", -1.0)):
                for scenario, value in (
                    (small, 8.0 + repetition),
                    (large, 10.0 + repetition),
                ):
                    row = self.row(family, "Adaptive", repetition, value)
                    row["scenario"] = scenario
                    row["disp_z_pct_hs"] = value + offset
                    rows.append(row)

        aggregate = validation.stationary_normalized_aggregate(rows, 500, 7)
        self.assertEqual(aggregate["OU_III_minus_OU_II"]["n_pairs"], 3)
        self.assertAlmostEqual(
            aggregate["OU_III_minus_OU_II"]["mean_paired_difference"],
            -1.0,
        )


class CommittedFullResultsTests(unittest.TestCase):
    RESULTS = REPO_ROOT / "reports" / "results" / "ou_validation"

    @staticmethod
    def read_csv(path):
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    def test_full_result_bundle_is_complete_and_self_consistent(self):
        manifest_path = self.RESULTS / "ou_validation_manifest.json"
        with manifest_path.open(encoding="utf-8") as stream:
            manifest = json.load(stream)
        self.assertEqual(manifest["protocol"]["mode"], "full")
        self.assertEqual(len(manifest["protocol"]["seed_triplets"]), 10)
        self.assertEqual(manifest["protocol"]["score_window_sec"], 900.0)
        self.assertIn(
            "analytically derived", manifest["protocol"]["wave_phase_method"]
        )
        self.assertIn(
            "first- and second-derivative",
            manifest["protocol"]["transition_method"],
        )

        raw = self.read_csv(self.RESULTS / "ou_validation_raw.csv")
        summary = self.read_csv(self.RESULTS / "ou_validation_summary.csv")
        effects = self.read_csv(
            self.RESULTS / "ou_validation_paired_effects.csv"
        )
        protocol = manifest["protocol"]
        run_modes = set(protocol["adaptation_modes"])
        pmstokes_modes = set(protocol["pmstokes_modes"])
        scenarios = {row["scenario"] for row in raw}
        families = {row["family"] for row in raw}
        self.assertEqual(families, {"OU_II", "OU_III"})
        # Four JONSWAP seas, four PM-Stokes seas, one controlled transition.
        self.assertEqual(
            Counter(validation.scenario_spectrum(s) for s in scenarios),
            Counter({"jonswap": 5, "pmstokes": 4}),
        )

        def expected_modes(scenario, family):
            allowed = set(validation.FAMILY_MODES[family]) & run_modes
            if validation.scenario_spectrum(scenario) == "pmstokes":
                allowed &= pmstokes_modes
            return allowed

        # The mode grid is family- and spectrum-dependent by design: the
        # channel ablation exists only for OU-III, and PM-Stokes is a
        # secondary ensemble scored on the primary modes only.
        cells = 0
        for scenario in sorted(scenarios):
            for family in sorted(families):
                expected = expected_modes(scenario, family)
                observed = {
                    row["mode"] for row in raw
                    if row["scenario"] == scenario and row["family"] == family
                }
                self.assertEqual(observed, expected, (scenario, family))
                cells += len(expected)
        self.assertGreater(cells, 0)
        self.assertEqual(len(raw), cells * 10)

        # The bundle carries one summary row per cell per metric. It is checked
        # against the metrics the bundle itself contains rather than against the
        # tool's current list, so that adding a metric asks for a re-run instead
        # of failing this structural check. Unknown metrics still fail.
        bundle_metrics = {row["metric"] for row in summary}
        self.assertTrue(
            bundle_metrics <= set(validation.METRIC_NAMES),
            sorted(bundle_metrics - set(validation.METRIC_NAMES)),
        )
        self.assertEqual(len(summary), cells * len(bundle_metrics))

        # Comparisons are emitted only where both arms were actually run.
        expected_comparisons = set()
        for scenario in scenarios:
            for family in families:
                modes = expected_modes(scenario, family)
                if "Adaptive" in modes:
                    for baseline in ("FixedNominal", "FixedOracle"):
                        if baseline in modes:
                            expected_comparisons.add(
                                (scenario, family, f"Adaptive_minus_{baseline}")
                            )
                if "FixedNominal" in modes:
                    for partial in validation.PARTIAL_MODES:
                        if partial in modes:
                            expected_comparisons.add(
                                (
                                    scenario,
                                    family,
                                    f"{partial}_minus_FixedNominal",
                                )
                            )
                for left, right in validation.COVARIANCE_SYNC_PAIRS:
                    if left in modes and right in modes:
                        expected_comparisons.add(
                            (scenario, family, f"{left}_minus_{right}")
                        )
            expected_comparisons.add(
                (scenario, "OU_III_vs_OU_II", "OU_III_minus_OU_II")
            )
        self.assertEqual(
            {
                (row["scenario"], row["family"], row["comparison"])
                for row in effects
            },
            expected_comparisons,
        )

        groups = Counter(
            (row["scenario"], row["family"], row["mode"])
            for row in raw
        )
        self.assertEqual(len(groups), cells)
        self.assertEqual(set(groups.values()), {10})
        # Segment metrics exist only where extra scoring intervals were
        # requested, which is the non-stationary scenario.
        for row in summary:
            segmented = row["metric"].startswith("seg_")
            transition = row["scenario"].startswith("nonstationary_")
            self.assertEqual(
                int(row["n"]), 0 if segmented and not transition else 10,
                (row["scenario"], row["metric"]),
            )
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

        for name, metadata in manifest["result_files"].items():
            path = self.RESULTS / name
            self.assertTrue(path.is_file(), name)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(digest, metadata["sha256"], name)
            self.assertEqual(path.stat().st_size, metadata["bytes"], name)

        for result_name, doc_name in (
            (
                "ou_validation_publication.tex",
                "w3d-ou-validation-results-generated.tex-part",
            ),
            (
                "ou_validation_macros.tex",
                "w3d-ou-validation-macros-generated.tex-part",
            ),
            (
                "ou_validation_tuning_points.tex",
                "w3d-ou-validation-tuning-points-generated.tex-part",
            ),
            ("ou_validation_vertical.svg", "ou_validation_vertical.svg"),
            ("ou_validation_transition.svg", "ou_validation_transition.svg"),
        ):
            self.assertEqual(
                (self.RESULTS / result_name).read_bytes(),
                (REPO_ROOT / "doc" / "kalman_ou_iii" / doc_name).read_bytes(),
                doc_name,
            )

    def test_every_generated_macro_used_in_the_manuscript_resolves(self):
        # A generated macro that is quoted but never defined does not stop the
        # build; it renders as an undefined-control-sequence artifact or, if
        # its name contains a digit, as the value with a stray text tail.  That
        # is how "27NominalDifference-13.118" reached the rendered paper, so
        # the reference and definition sets are checked against each other.
        doc_dir = REPO_ROOT / "doc" / "kalman_ou_iii"
        sources = sorted(doc_dir.glob("*.tex-part")) + [
            doc_dir / "kalman_ou-w3d.tex"
        ]
        defined = set()
        referenced = set()
        definition = re.compile(
            r"\\(?:provide|new|renew)command\s*\{\\([A-Za-z]+)\}"
        )
        reference = re.compile(r"\\((?:OUValidation|OURobustness)[A-Za-z]*)")
        for source in sources:
            text = source.read_text(encoding="utf-8")
            defined.update(definition.findall(text))
            for name in reference.findall(text):
                referenced.add(name)
        # Names quoted in prose but defined nowhere.
        self.assertEqual(sorted(referenced - defined), [])
        # And no definition may carry a character TeX cannot put in a name.
        for source in sources:
            validation.assert_latex_macro_names(
                source.read_text(encoding="utf-8").splitlines()
            )

    def test_transition_protocol_records_its_own_confounds(self):
        with (self.RESULTS / "ou_validation_manifest.json").open(
            encoding="utf-8"
        ) as stream:
            manifest = json.load(stream)
        protocol = manifest["protocol"]

        # The crossfade is not a spectral ramp, and the scoring window is not a
        # uniform sample of transitioning conditions. Both must stay recorded
        # with the evidence, not only asserted in prose.
        self.assertIn("crossfade", protocol["transition_method"])
        self.assertIn(
            "not a continuously evolving", protocol["transition_method"]
        )
        composition = protocol["transition_window_composition_sec"]
        self.assertGreater(composition["pure_end_sea_sec"], 0.0)
        self.assertAlmostEqual(
            composition["pure_start_sea_sec"]
            + composition["blended_sec"]
            + composition["pure_end_sea_sec"],
            protocol["score_window_sec"],
        )
        self.assertLess(
            protocol["transition_midpoint_mixture_hs_m"],
            protocol["transition_midpoint_linear_hs_m"],
        )

        # FixedOracle must not be described as an optimum.
        self.assertIn("not optimized", protocol["fixed_oracle_definition"])

        # The transition is also scored per interval, because one number
        # normalized by the final Hs is not comparable with a stationary score.
        self.assertEqual(
            set(protocol["transition_segments_sec"]),
            set(validation.TRANSITION_SEGMENTS),
        )
        bounds = protocol["transition_segments_sec"]
        self.assertEqual(bounds["blend"][0], protocol["transition_start_sec"])
        self.assertEqual(bounds["blend"][1], protocol["transition_end_sec"])
        self.assertEqual(bounds["start"][1], bounds["blend"][0])
        self.assertEqual(bounds["end"][0], bounds["blend"][1])

        # Interval construction and the multiplicity position must travel with
        # the evidence, not only with the prose.
        self.assertIn("percentile bootstrap", protocol["bootstrap_method"])
        self.assertIn("resample", protocol["bootstrap_method"])
        self.assertIn("PCG64", protocol["bootstrap_rng"])
        self.assertGreaterEqual(protocol["bootstrap_resamples"], 10_000)
        self.assertIn("not adjusted", protocol["multiplicity"])
        self.assertIn("JONSWAP", protocol["primary_endpoint"])
        self.assertIn("not pooled", protocol["pmstokes_pooling"])

        # The channel ablation is only meaningful if r_S keeps tracking the
        # live estimates while the OU channel is frozen.
        self.assertIn("live tau", protocol["partial_adaptation_definition"])

        policies = protocol["aw_covariance_sync_policies"]
        for left, right in validation.COVARIANCE_SYNC_PAIRS:
            self.assertEqual(policies[left], "periodic")
            self.assertEqual(policies[right], "reconfigure")

    def test_abstract_reports_committed_stationary_aggregate(self):
        with (self.RESULTS / "ou_validation_manifest.json").open(
            encoding="utf-8"
        ) as stream:
            aggregate = json.load(stream)["stationary_normalized_aggregate"]

        manuscript = (
            REPO_ROOT / "doc/kalman_ou_iii/kalman_ou-w3d.tex"
        ).read_text(encoding="utf-8")
        abstract = manuscript.split("\\begin{abstract}", 1)[1].split(
            "\\end{abstract}", 1
        )[0]
        ou2 = aggregate["OU_II"]
        ou3 = aggregate["OU_III"]
        difference = aggregate["OU_III_minus_OU_II"]

        self.assertIn(f"${ou3['mean']:.2f}\\pm{ou3['std']:.2f}\\%$", abstract)
        self.assertIn(f"${ou2['mean']:.2f}\\pm{ou2['std']:.2f}\\%$", abstract)
        self.assertIn(
            f"${difference['mean_paired_difference']:.3f}$ percentage point",
            abstract,
        )
        self.assertIn(
            "$["
            f"{difference['bootstrap_ci95_low']:.3f},"
            f"{difference['bootstrap_ci95_high']:.3f}"
            "]$",
            abstract,
        )

        # A reader must not be able to mistake the quoted spread for a
        # confidence interval, and the interval must name its construction.
        self.assertIn("sample standard deviation", abstract)
        self.assertIn("percentile-bootstrap", abstract)
        # The negative three-dimensional result and the channel ablation both
        # bound the contribution, so neither may be left to the body text.
        self.assertIn("full-3D displacement error", abstract)
        self.assertIn("channel ablation", abstract)
        self.assertIn("tolerance", abstract)


class ManuscriptMethodologyTests(unittest.TestCase):
    DOC = REPO_ROOT / "doc" / "kalman_ou_iii"

    @classmethod
    def read(cls, name):
        return (cls.DOC / name).read_text(encoding="utf-8")

    @classmethod
    def read_flat(cls, name):
        """Source with line wrapping collapsed, for phrase assertions."""

        return " ".join(cls.read(name).split())

    @staticmethod
    def paragraph(flat_source, opening):
        """One \\paragraph block of a flattened source, by its opening words."""

        start = flat_source.index(opening)
        rest = flat_source[start:]
        end = rest.find("\\paragraph{", 1)
        return rest if end < 0 else rest[:end]

    @staticmethod
    def three_d_verdict_counts():
        """Per-scenario 3D verdicts over the primary set, from the bundle.

        The primary set is the four stationary JONSWAP seas plus the
        transition; PM--Stokes is reported separately and is not pooled into
        it.  A scenario counts as decided only when its bootstrap interval
        excludes zero.
        """
        path = (
            REPO_ROOT / "reports" / "results" / "ou_validation"
            / "ou_validation_paired_effects.csv"
        )
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))

        higher = lower = unresolved = 0
        for row in rows:
            if row["comparison"] != "OU_III_minus_OU_II":
                continue
            if row["metric"] != "disp_3d_pct_refmax":
                continue
            if "pmstokes" in row["scenario"]:
                continue
            low = float(row["bootstrap_ci95_low"])
            high = float(row["bootstrap_ci95_high"])
            if low > 0.0:
                higher += 1
            elif high < 0.0:
                lower += 1
            else:
                unresolved += 1
        return higher, lower, unresolved

    def test_fixed_reference_and_transition_limits_are_stated(self):
        protocol = self.read("w3d-sim-charts.tex-part")
        results = self.read("w3d-baseline-comparison.tex-part")

        # "Oracle" must be defined as a scenario-calibrated fixed reference,
        # not as an optimum, and its exact values must be tabulated.
        self.assertIn("scenario-calibrated fixed reference", protocol)
        self.assertIn("searching against reference displacement error", protocol)
        self.assertIn("tab:ou_fixed_points", protocol)
        self.assertIn(
            "w3d-ou-validation-tuning-points-generated.tex-part", results
        )

        # The crossfade and its consequences must be stated, not implied.
        self.assertIn("crossfade", protocol)
        self.assertIn("2.14", protocol)
        self.assertIn("bimodal", protocol)
        self.assertIn("does not isolate adaptation rate", protocol)

        # The covariance-reset factor must be described as a controlled
        # ablation and used to qualify the adaptation result.
        self.assertIn("AdaptiveHeldCovariance", protocol)
        self.assertIn("tab:ou_mc_covsync", protocol)
        self.assertIn("OUValidationCovSyncWorstDifference", results)
        self.assertIn("fig:ou_transition", results)

        # The unsupported lag attribution must be gone.
        self.assertNotIn("transition lag", protocol + results)

    def test_inference_is_qualified_rather_than_asserted(self):
        protocol = self.read_flat("w3d-sim-charts.tex-part")
        results = self.read_flat("w3d-baseline-comparison.tex-part")
        robustness = self.read_flat("w3d-ou-robustness.tex-part")

        # No power or precision calculation is supplied, so the protocol may
        # not describe itself as powered.
        self.assertNotIn("statistically powered", protocol + results)

        # Interval construction, the pairing level, and the multiplicity
        # position must all be stated.  The method string, resample count, and
        # RNG seed come from the generated macros so they cannot drift from
        # the settings that actually produced the intervals.
        self.assertIn("OUValidationBootstrapMethod", protocol)
        self.assertIn("OUValidationBootstrapRNG", protocol)
        self.assertIn("OUValidationBootstrapResamples", protocol)
        self.assertIn("OUValidationStatsSeed", protocol)
        self.assertIn("resample the seed-level differences", protocol)
        self.assertIn("multiplicity", protocol)
        self.assertIn("primary endpoint", protocol)

        # A quoted a +- b is a spread, and the standardized effect size must
        # stay subordinate to the raw paired difference at n = 10.
        self.assertIn("not an interval estimate for the mean", protocol)
        self.assertIn("deliberately secondary", protocol)
        self.assertIn("primary statement of effect size", results)

        # Wave energy scales with Hs^2, so 1.5 m against 0.27 m or 8.5 m is a
        # factor of about 31, not two orders of magnitude.
        self.assertNotIn("two orders of magnitude", results)
        self.assertIn("factor of about $31$", results)

        # Compressing the ramp changes what is scored as well as how fast the
        # sea moves, so it cannot be read as a pure response-limit measurement.
        self.assertIn("confounds the two", robustness)
        self.assertNotIn(
            "These results quantify, rather than conceal, the response limit",
            robustness,
        )

    def test_the_sixty_second_table_is_generated_not_typed(self):
        """The deterministic 60 s table has to come from the simulator.

        It was hand-typed for a long time, and hand-typed numbers do not follow
        the code: its roll column still reported values from before the startup
        tilt reference was fixed, so the manuscript quoted an error the filter
        no longer made.  This asserts the protocol section inputs the generated
        file and carries no table body of its own.
        """
        source = self.read("w3d-sim-charts.tex-part")

        self.assertIn("w3d-sim-results-generated.tex-part", source)

        generated = self.read("w3d-sim-results-generated.tex-part")
        self.assertIn("tab:sim_results_rms_dir_basic", generated)
        self.assertIn("Generated by tools/ou_sim_table.py", generated)

        # The label must live in exactly one place, or the fallback and the
        # generated table would both define it.
        body = re.compile(r"^\s*(JONSWAP|PM--Stokes)\s*&", re.M)
        self.assertFalse(
            body.search(source),
            "w3d-sim-charts.tex-part still carries hand-typed table rows",
        )
        self.assertEqual(len(body.findall(generated)), 8)

    def test_three_dimensional_and_channel_results_are_reported(self):
        protocol = self.read_flat("w3d-sim-charts.tex-part")
        results = self.read_flat("w3d-baseline-comparison.tex-part")
        conclusion = self.read_flat("w3d-conclusion-summary.tex-part")

        # Absolute per-axis and 3D error for both filters, not only the
        # difference, so the magnitude of the negative result is readable.
        self.assertIn("tab:ou_mc_axes", results)
        self.assertIn("higher", results)

        # The 3D result must be stated with its remaining exception named, in
        # whichever direction the evidence puts it, and must not be rounded up
        # into a claim in either direction.
        #
        # The direction is read from the committed bundle rather than pinned as
        # a phrase.  It has already flipped once -- the fifth scenario was
        # unresolved and is now resolvably in OU--III's favour -- and while it
        # was pinned the assertion kept passing against prose that said the
        # opposite of the evidence beside it.
        higher, lower, unresolved = self.three_d_verdict_counts()
        self.assertEqual(higher + lower + unresolved, 5)

        self.assertNotIn("consistently higher OU--III 3D error", conclusion)
        self.assertNotIn("uniformly better in three dimensions", conclusion)
        self.assertIn("OUValidationThreeDHigherCount", conclusion)

        # Scoped to the 3D paragraph: "spanning zero" is a correct description
        # of other intervals elsewhere in the same section.
        three_d = self.paragraph(results, "The vertical gain does not extend")

        if unresolved:
            self.assertIn("spanning zero", three_d)
        else:
            # Nothing may be described as straddling zero when nothing does.
            self.assertNotIn("spanning zero", three_d)
            self.assertNotIn("inconclusive in the fifth", conclusion)

        if lower:
            # A scenario where OU--III wins in 3D has to be admitted as such,
            # not left implied by a count.
            self.assertIn("lower", conclusion)
            self.assertNotIn(
                "no scenario in which OU--III is resolvably better in three "
                "dimensions",
                three_d,
            )

        # The channel ablation and its no-implicit-freeze construction.
        self.assertIn("tab:ou_mc_channels", results)
        self.assertIn("par:channel-ablation", protocol)
        self.assertIn("live", protocol)
        self.assertIn("integral-state regularization", results)

    def test_transition_and_secondary_ensembles_are_rescored(self):
        protocol = self.read_flat("w3d-sim-charts.tex-part")
        results = self.read_flat("w3d-baseline-comparison.tex-part")

        # Normalizing the whole transition window by the final Hs is not
        # comparable with a stationary score; the interval breakdown and the
        # scale-free normalization must both be present.
        self.assertIn("tab:ou_transition_segments", protocol)
        self.assertIn("Z_{\\mathrm{ref}}", results)
        self.assertIn("not comparable with the stationary", results)

        # PM-Stokes and direction must be ensemble results, not single runs.
        self.assertIn("tab:ou_mc_pmstokes", results)
        self.assertIn("not pooled into the primary", results)
        self.assertIn("tab:ou_mc_direction", results)
        self.assertIn("OUValidationDirectionAbsError", results)
        # Travel sense is now a correctness rate, but only because it is scored
        # against the record's physical propagation direction.  The manuscript
        # must report that rate and must keep saying that the FORWARD/BACKWARD
        # classes themselves are a gauge quantity, since they invert under a
        # 180-degree heading change while the directed estimate does not.
        self.assertIn("OUValidationTravelCorrect", results)
        self.assertIn("generator azimuth plus", results)
        self.assertIn("are not scored", results)
        self.assertIn("heading rotated by", results)

    def test_singer_relationship_and_contribution_wording(self):
        intro = self.read("w3d-intro.tex-part")
        bibliography = self.read("w3d.bib")
        self.assertIn("Singer1970_ManeuverModel", intro)
        self.assertIn("We do not claim the OU", intro)
        self.assertIn("doi     = {10.1109/TAES.1970.310128}", bibliography)
        self.assertNotIn("novel", intro.lower())

    def test_nomenclature_is_included_and_disambiguates_covariances(self):
        manuscript = self.read("kalman_ou-w3d.tex")
        nomenclature = self.read("w3d-nomenclature.tex-part")
        self.assertIn("\\input{w3d-nomenclature.tex-part}", manuscript)
        self.assertIn("\\Sigma}_{aw}^{\\mathrm{stat}}", nomenclature)
        self.assertIn("\\mat{P}_{a_wa_w}", nomenclature)
        self.assertIn("Kalman gain", nomenclature)

    def test_analytic_completion_and_small_step_regime_are_explicit(self):
        analytic = self.read("w3d-analytic-coeff.tex-part")
        lti = self.read("w3d-lti-discrete.tex-part")
        manuscript = self.read("kalman_ou-w3d.tex")
        self.assertIn("\\scomp(\\mat{U})", analytic)
        self.assertIn("it is not the\naveraging operator", analytic)
        self.assertNotIn("\\sym", manuscript + analytic)
        self.assertIn("5.30\\times10^{-3}", analytic)
        self.assertIn("h-s", lti)
        self.assertNotIn("h-\\tau", lti)

    def test_tracker_and_heel_material_are_scoped_as_unevaluated_appendices(self):
        manuscript = self.read("kalman_ou-w3d.tex")
        fusion = self.read("w3d-fus-methods.tex-part")
        trackers = self.read("w3d-tracker-alternatives.tex-part")
        heel = self.read("w3d-wind-heel.tex-part")
        appendix_sources = "".join(
            self.read(name)
            for name in (
                "w3d-tracker-alternatives.tex-part",
                "w3d-wind-heel.tex-part",
                "w3d-gps-fusion.tex-part",
                "w3d-iss-stability.tex-part",
            )
        )
        self.assertLess(
            manuscript.index("\\appendices"),
            manuscript.index("\\input{w3d-tracker-alternatives.tex-part}"),
        )
        self.assertEqual(manuscript.count("\\appendices"), 1)
        self.assertNotIn("\\appendices", appendix_sources)
        self.assertNotIn("Aranovskiy Frequency Tracker", fusion)
        self.assertIn("not evaluated as OU--III tracker ablations", trackers)
        self.assertIn("disabled in every deterministic simulation", heel)

    def test_baseline_fairness_thresholds_and_hardware_limits_are_recorded(self):
        baseline = self.read("w3d-baseline-comparison.tex-part")
        fusion = self.read("w3d-fus-methods.tex-part")
        simulation = self.read("w3d-sim-charts.tex-part")
        self.assertIn("frozen before the comparison", baseline)
        self.assertIn("performs no parameter search", baseline)
        self.assertIn("never uses\nreference errors to choose gains", baseline)
        self.assertIn("tab:baseline-tuning-policy", baseline)
        self.assertIn("tab:implementation-gates", fusion)
        # The implementation bounds the table records have to be the ones the
        # filter actually enforces, so they are read from the header rather
        # than pinned as literals here.
        header = (
            REPO_ROOT / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
        ).read_text(encoding="utf-8")

        def clamp(name: str) -> float:
            match = re.search(rf"{name}\s*=\s*([0-9.]+)f", header)
            self.assertIsNotNone(match, name)
            return float(match.group(1))

        def tabulated(label: str) -> tuple[float, float]:
            match = re.search(
                re.escape(label) + r".*?\$\[([0-9.]+),([0-9.]+)\]\$",
                fusion,
                re.S,
            )
            self.assertIsNotNone(match, label)
            return float(match.group(1)), float(match.group(2))

        # Anchored on the symbol pair the row defines, which is unique, and
        # compared numerically so that 12 and 12.0 are the same bound.
        self.assertEqual(
            tabulated(r"[\tau_{\min},\tau_{\max}]$"),
            (clamp("MIN_TAU_S"), clamp("MAX_TAU_S")),
        )
        self.assertEqual(
            tabulated(r"[r_{S,\min},r_{S,\max}]$"),
            (clamp("MIN_R_S"), clamp("MAX_R_S")),
        )
        for value in ("[0.2,6.0]", "70^\\circ"):
            self.assertIn(value, fusion)
        self.assertIn("did not\ninstrument update latency", simulation)
        self.assertIn("does not establish a\nquantitative real-time-performance claim", simulation)


STUB_SIMULATOR = """#!/bin/sh
echo "XYZ RMS (m): X=0.0372370929 Y=0.0284884889 Z=0.0186269842"
echo "XYZ RMS (%Hs): X=13.7915154% Y=10.5512915% Z=6.89888287% (Hs=0.27)"
echo "Angles RMS (deg): Roll=0.562675 Pitch=0.157529 Yaw=1.60768"
echo "dir_deg_gen ([-90,90], 0=+Y CW): mean_circ=-28.275 circ_std 9.21715 deg"
echo "sense: +AXIS=11199 (93.325%) -AXIS=251 (2.09167%) UNCERTAIN=550 (4.58333%)"
"""


class SixtySecondTableToolTests(unittest.TestCase):
    def test_a_relative_binary_path_survives_the_child_directory_change(self):
        """The table generator has to accept the path the workflow passes.

        `build.yml` calls it from the repository root with
        `--binary tests/kalman_ou_iii/kalman_ou_iii-sim`, and the child is run
        from the simulator's own directory.  An unresolved relative path is then
        looked up under that directory instead, one level too deep, and the run
        dies reporting a missing simulator that is sitting exactly where the
        argument said it was.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            relative_binary = Path("tests") / "kalman_ou_iii" / "kalman_ou_iii-sim"
            binary = root / relative_binary
            binary.parent.mkdir(parents=True)
            binary.write_text(STUB_SIMULATOR, encoding="utf-8")
            binary.chmod(0o755)
            (root / "record.csv").write_text("", encoding="utf-8")

            previous = os.getcwd()
            os.chdir(root)
            try:
                row = sim_table.run_record(relative_binary, Path("record.csv"))
            finally:
                os.chdir(previous)

        self.assertEqual(row["z_rms_m"], (0.0186269842,))
        self.assertEqual(row["angles"], (0.562675, 0.157529, 1.60768))


if __name__ == "__main__":
    unittest.main()
