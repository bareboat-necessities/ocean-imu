"""Fail-closed attachment checks for the sampled runtime audit."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tests/kalman_ou_iii'))
import ou3_brmm_runtime_audit as audit  # noqa: E402
import ou3_brmm_runtime_overlay as overlay  # noqa: E402


def fixture():
    n = 2200
    source = pd.DataFrame(np.zeros((n, 6)), columns=['acc_'+a for a in 'xyz']+['gyro_'+a for a in 'xyz'])
    flags = ['outer_pre','outer_post','inner_pre','inner_post','bias_pre','bias_post',
             'predict','acc_calls','acc_applied','guard_calls','guard_identical','mag_lock','mag_refined']
    samples = pd.DataFrame({key: np.ones(n, dtype=int) for key in flags})
    samples['index'], samples['dt'] = np.arange(n), float(np.float32(.005))
    for key in ['mag_calls','mag_applied','s_calls','s_applied','tilt_resets','reference_writes','engagement','excess']:
        samples[key] = 0.0 if key in ('engagement', 'excess') else 0
    samples['resets'] = 1
    domain = json.loads(audit.DOMAIN.read_text())
    events = []
    for i in range(n):
        kinds = ['acc'] + (['S'] if i % 10 == 0 else []) + (['mag'] if i % 8 == 0 else [])
        for kind in kinds:
            vector = [0, 0, -9.8] if kind == 'acc' else [40, 0, 20]
            r = np.eye(3)*({'acc':.04, 'mag':.09, 'S':1}[kind])
            events.append([i,len(events),kind,*vector,*vector,*r.ravel()])
            if kind != 'acc':
                column = 's' if kind == 'S' else 'mag'
                samples.loc[i, column+'_calls'] += 1
                samples.loc[i, column+'_applied'] += 1
                samples.loc[i, 'resets'] += 1
    events = pd.DataFrame(events, columns=['index','sequence','kind','bx','by','bz','wx','wy','wz']+audit.R_FIELDS)
    return source, samples, events, domain


class RuntimeAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = fixture()

    def test_all_half_open_intervals_and_window_boundaries(self):
        np.testing.assert_array_equal(audit.runs([1,1,0,1,0,1]), [[0,2],[3,4],[5,6]])
        np.testing.assert_array_equal(audit.all_window_samples([1,1,0,1,1], 2), [1,0,0,1])

    def test_live_does_not_hide_cap_or_config_failures(self):
        source, samples, events, domain = copy.deepcopy(self.data)
        source.loc[300:309, 'acc_x'] = 5
        source.loc[600:609, 'gyro_y'] = 1
        events.loc[(events.kind == 'acc') & (events['index'] == 700), 'R00'] = .001
        r, masks = audit.analyze(source, samples, events, domain)
        self.assertEqual(r['regimes']['A21']['samples'], 2200)
        self.assertEqual(r['violations']['physical_acceleration_above_cap']['A21']['samples'], 10)
        self.assertEqual(r['violations']['physical_body_rate_above_cap']['A21']['samples'], 10)
        self.assertTrue(masks['acc_configured_R_mismatch'][700])
        self.assertEqual(r['source_samples_pruned'], 0)

    def test_dropped_S_or_shifted_event_fails_closed(self):
        for defect in ['dropped', 'shifted']:
            source, samples, events, domain = copy.deepcopy(self.data)
            index = events.index[events.kind == 'S'][0]
            if defect == 'dropped':
                events = events.drop(index)
                events.sequence = np.arange(len(events))
            else:
                events.loc[index, 'index'] += 1
            with self.assertRaises(ValueError):
                audit.analyze(source, samples, events, domain)

    def test_phase_caps_keep_startup_and_later_live_exit_distinct(self):
        source, samples, events, domain = copy.deepcopy(self.data)
        samples.loc[:99, ['outer_pre', 'outer_post']] = 0
        samples.loc[500:509, ['inner_pre', 'inner_post']] = 0
        source.loc[50, 'acc_z'] = 12
        source.loc[200, 'acc_z'] = 15
        source.loc[505, 'acc_z'] = 20
        r, _ = audit.analyze(source, samples, events, domain)
        phases = r['physical_phase_extrema']
        self.assertEqual([phases[k]['samples'] for k in phases], [100, 2090, 10])
        for key, value in [('before_Live', 12), ('Live', 15), ('non_Live_after_entry', 20)]:
            self.assertAlmostEqual(phases[key]['acceleration_norm_max_upper_mps2'], value)
            self.assertEqual(phases[key]['acceleration_norm_above_g_samples'], 1)

    def test_hybrid_and_guard_are_visible_without_pruning(self):
        source, samples, events, domain = copy.deepcopy(self.data)
        samples.loc[500, 'bias_pre'] = 0
        samples.loc[600, 'tilt_resets'] = 1
        samples.loc[900:999, 'engagement'] = .1
        r, _ = audit.analyze(source, samples, events, domain)
        self.assertEqual(r['violations']['bias_mode_transition']['all']['samples'], 1)
        self.assertEqual(r['violations']['guard_not_dormant_transparent']['A21']['samples'], 100)
        self.assertEqual(r['violations']['hard_tilt_reset']['A21']['samples'], 1)

    def test_passing_observations_do_not_promote_proof(self):
        r, _ = audit.analyze(*self.data)
        self.assertEqual(r['passes_observed_sample_checks_only']['samples'], 2200)
        self.assertEqual(r['S_applied'], 220)
        self.assertFalse(r['P3_certified_by_this_audit'])
        self.assertFalse(r['continuous_BRMM_admission_certified'])
        self.assertFalse(r['transported_vector_PE_certified'])
        self.assertFalse(r['P4_or_P5_promoted'])

    def test_continuous_reference_write_is_not_automatic_hard_regauge(self):
        source, samples, events, domain = copy.deepcopy(self.data)
        samples.loc[100, 'reference_writes'] = 1
        r, _ = audit.analyze(source, samples, events, domain)
        self.assertEqual(r['passes_observed_sample_checks_only']['samples'], 2200)
        self.assertEqual(r['observations']['magnetic_reference_write']['samples'], 1)

    def test_overlay_is_insertions_only_and_bound_to_source(self):
        with tempfile.TemporaryDirectory() as temp:
            manifest = overlay.build(ROOT, Path(temp))
            audit.verify_manifest(manifest, ROOT)
            for path, entry in manifest['files'].items():
                generated = (Path(temp)/path).read_text()
                for insertion in entry['insertions']:
                    generated = generated.replace(insertion, '')
                self.assertEqual(generated, (ROOT/path).read_text())
            manifest['files'][overlay.MEKF]['shipping_sha256'] = '0'*64
            with self.assertRaises(ValueError):
                audit.verify_manifest(manifest, ROOT)
        with self.assertRaises(ValueError):
            overlay.instrument((ROOT/overlay.MEKF).read_text().replace('last_acc_diag_ = MeasDiag3{};', ''), overlay.MEKF)


if __name__ == '__main__':
    unittest.main()
