import contextlib
import csv
import hashlib
import io
import json
import math
import sys
import tempfile
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



class ShardMergeTests(unittest.TestCase):
    """The merge is what lets regeneration run on several runners at once.

    It has to put the rows back in exactly the order one process would have
    produced, because the raw CSV is hashed in the manifest and the bundle is
    the published evidence.  A merge that silently dropped or reordered rows
    would still write a plausible-looking bundle, so the failure modes are
    checked here rather than left to a reviewer noticing.
    """

    @staticmethod
    def shard_rows(count, per_shard):
        """Rows as the shards would emit them: interleaved by repetition."""
        shards = [[] for _ in range(count)]
        order = 0
        for repetition in range(count * 2):
            for _ in range(per_shard):
                shards[repetition % count].append(
                    {"run_id": f"r{repetition}", robustness.SHARD_ORDER_KEY: order}
                )
                order += 1
        return shards

    def write_shards(self, directory, shards):
        for index, rows in enumerate(shards):
            robustness.write_shard(
                robustness.shard_path(Path(directory), index), rows
            )

    def test_merge_restores_single_process_order(self):
        shards = self.shard_rows(3, 4)
        with tempfile.TemporaryDirectory() as directory:
            self.write_shards(directory, shards)
            merged = robustness.read_shards(Path(directory))

        # Position is scaffolding and must not reach the bundle.
        self.assertNotIn(robustness.SHARD_ORDER_KEY, merged[0])
        expected = [
            row["run_id"]
            for row in sorted(
                (row for shard in shards for row in shard),
                key=lambda row: row[robustness.SHARD_ORDER_KEY],
            )
        ]
        self.assertEqual([row["run_id"] for row in merged], expected)

    def test_a_missing_shard_is_an_error_rather_than_a_short_bundle(self):
        shards = self.shard_rows(3, 4)
        with tempfile.TemporaryDirectory() as directory:
            self.write_shards(directory, shards[:-1])
            with self.assertRaisesRegex(ValueError, "not a complete set"):
                robustness.read_shards(Path(directory))

    def test_a_repeated_shard_is_an_error(self):
        shards = self.shard_rows(2, 3)
        with tempfile.TemporaryDirectory() as directory:
            self.write_shards(directory, shards)
            robustness.write_shard(
                robustness.shard_path(Path(directory), 7), shards[0]
            )
            with self.assertRaisesRegex(ValueError, "duplicate rows"):
                robustness.read_shards(Path(directory))

    def test_no_shards_at_all_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                robustness.read_shards(Path(directory))

    def test_shard_round_trip_preserves_non_finite_metrics(self):
        # Several metrics are legitimately NaN, and the JSON writer used for
        # the bundle maps those to null.  If shard files went through it, a
        # sharded run would disagree with an unsharded one on exactly the rows
        # a degenerate case produces.
        rows = [
            {
                "disp_3d_pct_refmax": float("nan"),
                "dir_axis_error_deg": float("inf"),
                "disp_z_pct_hs": 4.25,
                "mode": "Adaptive",
                robustness.SHARD_ORDER_KEY: 0,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            robustness.write_shard(
                robustness.shard_path(Path(directory), 0), rows
            )
            merged = robustness.read_shards(Path(directory))

        self.assertTrue(math.isnan(merged[0]["disp_3d_pct_refmax"]))
        self.assertEqual(merged[0]["dir_axis_error_deg"], float("inf"))
        self.assertEqual(merged[0]["disp_z_pct_hs"], 4.25)
        self.assertEqual(merged[0]["mode"], "Adaptive")


class RestatTests(unittest.TestCase):
    """`--restat-from` has to restate a bundle, not quietly rewrite it.

    The replays are the expensive half of this study, and this study borrows
    its statistics from ``ou_validation.py``, so an edit over there strands the
    bundle's source pin without changing a single number.  Recomputing the
    tables from the archived rows is the honest way out of that, but only if a
    restat of the committed bundle reproduces the committed derived files byte
    for byte, which is what this checks.
    """

    RESULTS = REPO_ROOT / "reports" / "results" / "ou_robustness"

    def test_restating_the_committed_bundle_reproduces_its_derived_files(self):
        source = self.RESULTS / "ou_robustness.json"
        if not source.exists():  # pragma: no cover - smoke checkouts
            self.skipTest("no committed robustness bundle")

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            with contextlib.redirect_stdout(io.StringIO()):
                robustness.restat_bundle(
                    source,
                    output,
                    bootstrap_resamples=10_000,
                    stats_seed=20260329,
                )
            for name in (
                "ou_robustness_raw.csv",
                "ou_robustness_summary.csv",
                "ou_robustness_paired_effects.csv",
                "ou_robustness_publication.tex",
            ):
                self.assertEqual(
                    (output / name).read_bytes(),
                    (self.RESULTS / name).read_bytes(),
                    name,
                )

            restated = json.loads(
                (output / "ou_robustness.json").read_text(encoding="utf-8")
            )
            committed = json.loads(source.read_text(encoding="utf-8"))
            # Same rows, same operating point, and the restat must say what it
            # restated.
            self.assertEqual(restated["raw_runs"], committed["raw_runs"])
            self.assertEqual(
                restated["nominal_tuning_point"], committed["nominal_tuning_point"]
            )
            self.assertIn("restated_from", restated["protocol"])

    def test_a_restat_keeps_covering_the_files_it_did_not_rewrite(self):
        """A restat rewrites the tables; the figures stay as the run left them.

        Matplotlib's SVG output is not a function of the rows alone, so the
        figures are carried rather than redrawn.  The manifest is the bundle's
        inventory, and dropping the files a restat happens not to touch would
        quietly shrink what the hash check covers.
        """
        manifest = json.loads(
            (self.RESULTS / "ou_robustness_manifest.json").read_text(encoding="utf-8")
        )
        covered = set(manifest["result_files"])
        for name in (
            "ou_robustness_sensitivity.svg",
            "ou_robustness_stress.svg",
        ):
            self.assertIn(name, covered, name)
        # The code sources are pinned by whichever path wrote the manifest, so
        # a restat must not be a way to lose them.
        pinned = set(manifest["source_files"])
        for path in robustness.CODE_SOURCE_PATHS:
            self.assertIn(str(path.relative_to(REPO_ROOT)), pinned)

    def test_a_restat_records_the_sources_that_moved_under_it(self):
        """The rows come from the earlier run; the tables come from this tree.

        When those two disagree about a source file the manifest has to say so,
        because a moved simulator source means the replays are stale and no
        amount of restating fixes that.
        """
        previous = {
            "tools/ou_validation.py": {"sha256": "0" * 64, "bytes": 1},
            "plots/kalman_ou_ii/absent.csv": {"sha256": "1" * 64, "bytes": 2},
        }
        restated, moved = robustness._restated_source_files(previous)
        # Present and changed: re-pinned, and reported as moved.
        self.assertEqual(set(moved), {"tools/ou_validation.py"})
        self.assertEqual(moved["tools/ou_validation.py"], previous["tools/ou_validation.py"])
        self.assertNotEqual(restated["tools/ou_validation.py"]["sha256"], "0" * 64)
        # Absent from the checkout: carried across, never dropped.
        self.assertEqual(
            restated["plots/kalman_ou_ii/absent.csv"],
            previous["plots/kalman_ou_ii/absent.csv"],
        )

    def test_a_restat_cannot_manufacture_an_ensemble(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.json"
            empty.write_text(json.dumps({"raw_runs": []}), encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(ValueError):
                    robustness.restat_bundle(empty, Path(tmp), 100, 1)


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
        # Source pins are provenance, and the two kinds of source drift for
        # different reasons.
        #
        # The wave records are versioned release inputs.  They do not change
        # under editing, and a hash that moved would mean the bundle was scored
        # against different seas than the ones it names, so they are still
        # checked byte for byte.
        #
        # The code sources are the tree that produced the bundle, and that tree
        # moves whenever anyone touches the tools -- a comment, a docstring, a
        # crash fix on a path the full study never takes.  Asserting they still
        # match the checkout made every such edit fail this test, and the only
        # way to satisfy it was a multi-hour regeneration dispatch, so a
        # one-line fix could not be landed without re-running the study.  It
        # also contradicts the tools themselves: `_restated_source_files`
        # already treats a moved code hash as something to *record* rather than
        # something to reject.  What has to hold is that the manifest keeps
        # naming them with well-formed entries, so provenance cannot quietly
        # shrink to nothing.
        code_sources = {
            str(path.relative_to(REPO_ROOT))
            for path in robustness.CODE_SOURCE_PATHS
        }
        self.assertTrue(
            code_sources <= set(manifest["source_files"]),
            sorted(code_sources - set(manifest["source_files"])),
        )
        for name, metadata in manifest["source_files"].items():
            self.assertRegex(metadata["sha256"], r"\A[0-9a-f]{64}\Z", name)
            self.assertGreater(metadata["bytes"], 0, name)
            if name in code_sources:
                continue
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
