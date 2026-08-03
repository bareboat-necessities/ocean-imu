import csv
import hashlib
import json
import math
import sys
import unittest
from collections import Counter
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

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
        scenarios = {row["scenario"] for row in raw}
        families = {row["family"] for row in raw}
        modes = {row["mode"] for row in raw}
        self.assertEqual(len(scenarios), 5)
        self.assertEqual(families, {"OU_II", "OU_III"})
        self.assertEqual(modes, set(validation.MODE_SETTINGS))
        cells = len(scenarios) * len(families) * len(modes)
        self.assertEqual(len(raw), cells * 10)
        self.assertEqual(len(summary), cells * len(validation.METRIC_NAMES))
        # Per scenario: one cross-family comparison, plus, for each family,
        # two adaptation baselines and two covariance-reset controls.
        comparisons_per_scenario = 1 + len(families) * (
            2 + len(validation.COVARIANCE_SYNC_PAIRS)
        )
        self.assertEqual(
            len(effects),
            len(scenarios)
            * comparisons_per_scenario
            * len(validation.METRIC_NAMES),
        )
        groups = Counter(
            (row["scenario"], row["family"], row["mode"])
            for row in raw
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


class ManuscriptMethodologyTests(unittest.TestCase):
    DOC = REPO_ROOT / "doc" / "kalman_ou_iii"

    @classmethod
    def read(cls, name):
        return (cls.DOC / name).read_text(encoding="utf-8")

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
        for value in ("[0.2,6.0]", "[0.02,3.0]", "[0.4,35]", "70^\\circ"):
            self.assertIn(value, fusion)
        self.assertIn("did not\ninstrument update latency", simulation)
        self.assertIn("does not establish a\nquantitative real-time-performance claim", simulation)


if __name__ == "__main__":
    unittest.main()
