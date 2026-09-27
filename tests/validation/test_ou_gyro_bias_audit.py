"""Paired engineering evidence must reject corrupt or incomparable inputs."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou_gyro_bias_audit import compare_audit_directories, compare_metrics  # noqa: E402


class GyroAuditComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.before, self.after = [Path(self.temporary.name)/name for name in ("before", "after")]
        observation = dict(family=3, predictions=10, corrections=8, nonfinite=0,
            projection_count=0, max_bias_rad_s=.02, max_prediction_angle_rad=.003,
            max_projection_rad_s=0, max_dt_s=.005, candidate_exceedances=[0, 0, 0, 0])
        for directory, metric in ((self.before, 1.0), (self.after, 1.125)):
            (directory/"validation").mkdir(parents=True)
            audit = dict(observer_sha256="same-observer", source_sha256={"header": directory.name},
                candidate_radii_rad_s=[.1, .2, .5, 1.], runs=[], replay_study="validation",
                validation_exit_code=0, validation_replays=[dict(observations=[observation],
                    input_immutable_verified=True, input_identity="wave.csv", input_sha256="same-input",
                    input_rows=10, environment={"W3D_SEED": "1"})])
            study = dict(protocol={"mode": "full"}, raw_runs=[dict(run_id="case", family="OU_III",
                samples=10, disp_z_ref_rms_m=.5, disp_3d_rms_m=metric)])
            (directory/"audit.json").write_text(json.dumps(audit))
            (directory/"validation/ou_validation.json").write_text(json.dumps(study))

    def mutate(self, path, operation):
        data = json.loads(path.read_text())
        operation(data)
        path.write_text(json.dumps(data))

    def test_pairs_metrics_and_keeps_projection_observations(self):
        self.mutate(self.after/"audit.json", lambda d: d["validation_replays"][0]["observations"][0].update(
            projection_count=1, max_projection_rad_s=.1))
        result = compare_audit_directories(self.before, self.after)["replay"]
        self.assertEqual(result["paired_replay_count"], 1)
        self.assertEqual(result["metrics"]["max_absolute_metric_differences"]["OU_III"]["disp_3d_rms_m"], .125)
        self.assertEqual(result["modified"]["3"]["projection_count"], 1)
        self.assertEqual(result["modified"]["3"]["max_projection_rad_s"], .1)

    def test_rejects_changed_physical_input(self):
        self.mutate(self.after/"audit.json", lambda d: d["validation_replays"][0].update(input_sha256="changed"))
        with self.assertRaisesRegex(ValueError, "physical inputs"):
            compare_audit_directories(self.before, self.after)

    def test_rejects_unverified_or_incomplete_replays(self):
        self.mutate(self.after/"audit.json", lambda d: d["validation_replays"][0].update(input_immutable_verified=False))
        with self.assertRaisesRegex(ValueError, "immutability"):
            compare_audit_directories(self.before, self.after)
        self.mutate(self.after/"audit.json", lambda d: d.update(validation_exit_code=1))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            compare_audit_directories(self.before, self.after)

    def test_rejects_changed_reference_motion(self):
        self.mutate(self.after/"validation/ou_validation.json", lambda d: d["raw_runs"][0].update(disp_z_ref_rms_m=.6))
        with self.assertRaisesRegex(ValueError, "reference motion"):
            compare_audit_directories(self.before, self.after)

    def test_rejects_duplicate_metric_identities(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            compare_metrics([dict(run_id="same"), dict(run_id="same")], [dict(run_id="same")])

    def test_standalone_uses_and_verifies_retained_stdout(self):
        for directory, metric in ((self.before, 1.), (self.after, 1.125)):
            stdout = f"VALIDATION_METRICS input=wave family=OU_III rms={metric}\n"
            (directory/"sim.stdout").write_text(stdout)
            row = dict(target="sim", exit_code=0, test_source_sha256="same-test", observations=[],
                       stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest())
            self.mutate(directory/"audit.json", lambda d, row=row: d.update(runs=[row]))
        result = compare_audit_directories(self.before, self.after, part="standalone")["standalone"]
        self.assertEqual(result["metrics"][0]["max_absolute_metric_differences"]["OU_III"]["rms"], .125)
        (self.after/"sim.stdout").write_text("changed")
        with self.assertRaisesRegex(ValueError, "stdout changed"):
            compare_audit_directories(self.before, self.after, part="standalone")


if __name__ == "__main__":
    unittest.main()
