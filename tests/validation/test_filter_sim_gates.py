"""Accuracy reports must not silently omit failed, missing or invalid records."""
import importlib.util
import json
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
spec = importlib.util.spec_from_file_location('filter_sim_gates', ROOT / 'tools/filter_sim_gates.py')
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)


def fixture():
    lines = []
    for name in G.REFERENCE_CSV_SHA256:
        lines.append('VALIDATION_METRICS input=' + name +
                     ' start_s=300 window_s=900 samples=180000 ' +
                     ' '.join(key + '=0.1' for key in G.METRICS))
        for kind in ('acc', 'gyro'):
            lines.append(f'Bias error RMS (% of max TRUE bias) ({kind}): X=1% Y=2% Z=3% |3D|=4%')
    return '\n'.join(lines)


class ReplayEvidenceTests(unittest.TestCase):
    def test_rejects_partial_duplicate_nonfinite_and_wrong_window(self):
        text = fixture()
        self.assertEqual(len(G.parse_w3d(text)), 8)
        for invalid in (text.split('\n', 3)[-1], text + '\n' + text,
                        text.replace('=0.1', '=nan', 1),
                        text.replace('samples=180000', 'samples=12000', 1),
                        text.replace('Z=3%', 'Z=nan%', 1),
                        text.replace('start_s=300', 'start_s=0', 1)):
            with self.subTest(invalid=invalid[:80]):
                with self.assertRaises(ValueError):
                    G.parse_w3d(invalid)

    def test_normalized_bias_and_family_maxima_are_independent(self):
        rows = G.parse_w3d(fixture())
        rows[0]['disp_z_pct_hs'] = 7.0
        rows[-1]['gyro_bias_3d_percent'] = 9.0
        worst = G.measured_worst(rows)
        self.assertEqual(worst['err_limit_percent_z_jonswap'], 7.0)
        self.assertEqual(worst['err_limit_percent_z_pmstokes'], 0.1)
        self.assertEqual(worst['bias_3d_percent'], 4.0)
        self.assertEqual(worst['gyro_bias_3d_percent'], 9.0)

    def test_every_simulator_has_an_audit_owner(self):
        owned = {f'tests/{folder}/{binary}.cpp' for folder, binary, _ in G.REPLAYS.values()}
        owned.update(('tests/kalman_carma/carma_process-sim.cpp',
                      'tests/kalman_carma/carma_lever_arm-sim.cpp'))
        self.assertEqual(owned, {str(p.relative_to(ROOT)) for p in (ROOT / 'tests').rglob('*-sim.cpp')})

    def test_current_limits_cover_measurements_without_loosening_existing_gates(self):
        evidence = json.loads((ROOT / 'reports/results/filter_sim_gates/limits.json').read_text())
        header = (ROOT / 'tests/common/CurrentReplayLimits.h').read_text()
        for name, family in evidence['families'].items():
            block = re.search(r'W3dFailureLimits ' + name + r'\{(.*?)\};', header, re.S)[1]
            compiled = {key: float(value) for key, value in
                        re.findall(r'\.(\w+)\s*=\s*([0-9.]+)f', block)}
            self.assertEqual(compiled, family['current_limits'], name)
            self.assertEqual(set(compiled), set(G.GATE_METRICS), name)
            for gate, value in compiled.items():
                self.assertGreater(value, family['measured_worst'][gate], (name, gate))
                previous = family['previous_production_limits'].get(gate)
                if previous:
                    self.assertLessEqual(value, previous, (name, gate))

    def test_make_log_keeps_sampling_profiles_separate_and_rejects_failed_runs(self):
        families = ('OU_II', 'OU_III', 'TFG', 'CARMA3_A41', 'CARMA_SHARED_MODES_CSM',
                    'TFG_CARMA3_T41', 'TFG_CARMA_SHARED_MODES_CSM', 'TFG_LEGACY')
        blocks = []
        for family in families:
            for legacy in (False, True) if family != 'TFG_LEGACY' else (False,):
                mode = 'legacy_interval' if legacy else 'instantaneous_default'
                marker = ('TFG legacy comparison simulation starting' if family == 'TFG_LEGACY'
                          else 'REGRESSION_PROTOCOL: SAMPLING=' + mode if 'SHARED' in family
                          else 'GYRO_REPLAY_PROFILE ' + mode)
                data = fixture().replace('VALIDATION_METRICS input=',
                                         f'VALIDATION_METRICS family={family} input=')
                if legacy:
                    data = data.replace('disp_z_pct_hs=0.1', 'disp_z_pct_hs=9.0')
                blocks.append(marker + '\n' + data)
        log = '\n'.join(blocks)
        parts = G.split_make_log(log)
        self.assertEqual(G.parse_w3d(parts['ou2', 'legacy'])[0]['disp_z_pct_hs'], 9.0)
        self.assertEqual(G.parse_w3d(parts['ou2', 'default'])[0]['disp_z_pct_hs'], 0.1)
        for invalid in (log + '\n' + blocks[0], log + '\nQUALITY_GATE: PASS=0',
                        log + '\nERROR: raw Z RMS above limit'):
            with self.assertRaises(ValueError):
                G.split_make_log(invalid)


if __name__ == '__main__':
    unittest.main()
