"""Startup evidence must retain fixed-horizon trials when shipping times out."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import report_startup_learning as report
import startup_learning as experiment
import compare_startup_slope as slope
from startup_source_snapshot import SourceSnapshot


class StartupAccuracyTests(unittest.TestCase):
    def summary(self):
        cells = []
        for case in ('following', 'head'):
            for profile, method in (('calibrated', 'joint'), ('residual', 'joint_cal')):
                for name, schedule, yaw, tilt, bg, n in (
                    ('shipping', 'handoff', 3, 1.2, .05, 1),
                    (method, 'matched', .8, .3, .02, 1),
                    # A second trial has no shipping handover but MUST be
                    # included in the independent deadline accuracy checks.
                    (method, 'fixed', .9, .35, .02, 2),
                ):
                    cells.append(dict(case=case, profile=profile, method=name, schedule=schedule,
                                      horizon='35' if schedule == 'fixed' else 'handoff', n=n,
                                      yaw_rms_deg=yaw, tilt_rms_deg=tilt, bg_rms_dps=bg,
                                      worst_yaw_deg=1.5, worst_tilt_deg=.6))
        return cells

    def test_deadline_miss_cannot_hide_a_candidate_regression(self):
        summary = self.summary()
        self.assertTrue(all(c['passed'] for c in experiment.check_accuracy(summary)))
        # Matched-time RMS stays good, but the unmatched trial is inaccurate.
        endpoint = next(s for s in summary if s['case'] == 'following' and
                        s['profile'] == 'calibrated' and s['schedule'] == 'fixed')
        for metric, value in (('worst_yaw_deg', 2.1), ('worst_tilt_deg', 1.1), ('bg_rms_dps', .041)):
            with self.subTest(metric=metric):
                original = endpoint[metric]
                endpoint[metric] = value
                failed = [c for c in experiment.check_accuracy(summary) if not c['passed']]
                self.assertEqual([(c['schedule'], c['metric']) for c in failed], [('fixed', metric)])
                endpoint[metric] = original

    def test_missing_endpoint_is_an_error(self):
        summary = [s for s in self.summary() if s['schedule'] != 'fixed']
        with self.assertRaisesRegex(ValueError, 'missing 35-second accuracy cell'):
            experiment.check_accuracy(summary)


class StartupComparisonTests(unittest.TestCase):
    def test_source_or_revision_changes_reject_the_receipt(self):
        with tempfile.TemporaryDirectory() as temp, patch('startup_source_snapshot.subprocess.check_output',
                                                        return_value='original-commit\n') as revision:
            root = Path(temp)
            source = root / 'filter.h'
            source.write_text('original estimator')
            snapshot = SourceSnapshot(root, [source])
            snapshot.verify()
            source.write_text('changed estimator')
            with self.assertRaisesRegex(ValueError, 'filter.h'):
                snapshot.verify()
            source.write_text('original estimator')
            revision.return_value = 'different-commit\n'
            with self.assertRaisesRegex(ValueError, 'HEAD'):
                snapshot.verify()
            revision.return_value = 'original-commit\n'
            source.unlink()
            with self.assertRaisesRegex(ValueError, 'filter.h'):
                snapshot.verify()

    def test_different_or_missing_executable_cannot_be_attributed_to_weighting(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folders = [root / name for name in ('fixed', 'scaled')]
            for folder in folders:
                folder.mkdir()
                metrics = folder / 'metrics.csv.gz'
                with gzip.open(metrics, 'wt') as stream:
                    stream.write('method,profile,schedule,elapsed_s,case,input,seed,offset,epoch_s,valid,'
                                 'yaw_error_deg,tilt_error_deg,bg_error_dps\n'
                                 'joint,calibrated,fixed,35,following,wave.csv,7,0,35,1,.8,.3,.02\n')
                protocol = dict.fromkeys(('datasets', 'seeds', 'offsets_s', 'profiles', 'imu_hz',
                                          'mag_hz', 'deadline_s'), 'matched')
                protocol.update(force_slope_model=folder.name, metrics_file=metrics.name,
                                metrics_sha256=slope.digest(metrics), executable_sha256='same-binary')
                (folder / 'provenance.json').write_text(json.dumps(protocol))
            output = root / 'comparison.csv'
            slope.compare(*folders, output)
            self.assertEqual(len(output.read_text().splitlines()), 2)
            receipt = folders[1] / 'provenance.json'
            for replacement in ('different-binary', None):
                changed = json.loads(receipt.read_text())
                changed['executable_sha256'] = replacement
                receipt.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, 'same recorded executable'):
                    slope.compare(*folders, output)


class StartupReportTests(unittest.TestCase):
    def test_missing_shipping_handovers_still_render_candidates_and_deadlines(self):
        # Exercise both a missing first-row baseline (the legend must still
        # include the second-row marker) and an entirely absent baseline.
        for missing in ({'following'}, {'following', 'head'}):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                summary, coverage = [], []
                for case in ('following', 'head'):
                    for method in ('joint', 'joint_cal'):
                        for horizon in ('20', '35'):
                            summary.append(dict(
                                case=case, profile='calibrated', method=method, schedule='fixed',
                                horizon=horizon, n=1, yaw_rms_deg=.7, tilt_rms_deg=.3, bg_rms_dps=.02,
                                worst_yaw_deg=.7, worst_tilt_deg=.3,
                                ready_fraction=1, false_ready_fraction=0))
                    if case not in missing:
                        summary.append(dict(
                            case=case, profile='calibrated', method='shipping', schedule='handoff',
                            horizon='handoff', n=1, mean_elapsed_s=31.84,
                            yaw_rms_deg=3, tilt_rms_deg=1, bg_rms_dps=.05,
                            field_direction_rms_deg=.5, field_magnitude_rms_uT=.03))
                    entry = dict(case=case, profile='calibrated', shipping_deadline_missed=case in missing,
                                 shipping_handover_s=None if case in missing else 31.84)
                    for method in ('joint', 'joint_cal'):
                        entry.update({f'{method}_first_ready_s': 20,
                                      f'{method}_accuracy_by_s': 20,
                                      f'{method}_false_first_ready': False})
                    coverage.append(entry)
                for name, data in (
                    ('summary', summary), ('coverage', coverage),
                    ('provenance', dict(base_commit='fixture', seeds='1', offsets_s='0')),
                ):
                    (root / f'{name}.json').write_text(json.dumps(data))
                with gzip.open(root / 'metrics.csv.gz', 'wt') as stream:
                    stream.write('method,schedule,elapsed_s,solve_ms\njoint,fixed,35,10\n')
                with patch.object(sys, 'argv', ['report_startup_learning.py', str(root)]):
                    report.main()
                text = (root / 'REPORT.md').read_text()
                svg = (root / 'convergence.svg').read_text()
                self.assertIn('| following | calibrated | joint | 1 | 1 | — | 1/1 |', text)
                self.assertIn('joint', svg)
                self.assertIn('Head seas', svg)
                self.assertEqual('shipping at its handover' in svg, len(missing) == 1)
                self.assertTrue((root / 'convergence.png').read_bytes().startswith(b'\x89PNG'))


if __name__ == '__main__':
    unittest.main()
