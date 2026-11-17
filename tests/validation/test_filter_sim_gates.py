"""Accuracy reports must not silently omit failed, missing or invalid records."""
import importlib.util
import json
from pathlib import Path
import re
import sys
import subprocess
import tempfile
import unittest
from unittest import mock

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


def pii_fixture():
    return '\n'.join('PII_METRICS input=' + name +
        ' start_s=300 window_s=900 samples=180000 disp_z_rms_m=0.1 disp_z_pct_hs=1'
        ' roll_rms_deg=0.2 pitch_rms_deg=0.3 yaw_rms_deg=0.4\nQUALITY_GATE: PASS=1 REASON=ok'
        for name in G.REFERENCE_CSV_SHA256)


class ReplayEvidenceTests(unittest.TestCase):
    def test_cli_rejects_failed_gate_even_when_process_returns_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            logs = Path(tmp)
            runs = []
            for name in G.REPLAYS:
                profiles = []
                for profile, _ in G.profiles_for(name):
                    if name == 'pii':
                        text = pii_fixture()
                    elif name == 'ahrs':
                        text = 'Angles RMS (deg): Roll=1 Pitch=1 Yaw=1\n' * 8
                    elif name == 'nlo':
                        text = 'Z RMS raw (%Hs): 1%\n' * 8
                    elif name == 'heave':
                        methods = G.HEAVE_METHODS | ({'pii'} if profile != 'truth' else set())
                        frontend = 'mahony' if profile == 'default' else profile
                        text = '\n'.join(f'HEAVE_ONLY method={method} frontend={frontend} record={record} z_pct_hs=1'
                                         for method in methods for record in G.REFERENCE_CSV_SHA256)
                    else:
                        lines = fixture().splitlines()
                        text = '\n'.join(line + ('\nQUALITY_GATE: PASS=1 REASON=ok' if line.startswith('Bias error RMS (% of max TRUE bias) (gyro)') else '')
                                         for line in lines)
                    if (name, profile) == ('ou3', 'default'):
                        text = text.replace('QUALITY_GATE: PASS=1 REASON=ok',
                                            'QUALITY_GATE: PASS=0 REASON=accuracy', 1)
                    path = logs / f'{name}-{profile}.log'
                    path.write_text(text)
                    profiles.append(dict(profile=profile, exit_code=0, log_sha256=G.digest(path)))
                runs.append(dict(name=name, binary_sha256='fixture', profiles=profiles))
            (logs / 'runs.json').write_text(json.dumps(runs))
            output = logs / 'measurements.json'
            with mock.patch.object(G, 'input_provenance', return_value={}), mock.patch.object(
                    sys, 'argv', ['filter_sim_gates.py', '--log-dir', str(logs), '--output', str(output)]):
                with self.assertRaisesRegex(SystemExit, 'replay gates failed'):
                    G.main()
            report = json.loads(output.read_text())
            self.assertNotEqual(report['acceptance'], 'all replay gates passed')
            self.assertIn('ou3/default', report['gate_error'])
            self.assertEqual(report['replay_exit_codes']['ou3']['default'], 0)

    def test_per_record_inventory_caps_and_generated_cpp(self):
        from generate_record_replay_limits import render
        from lock_simulation_baseline import cap
        data = json.loads((ROOT / 'reports/results/filter_sim_gates/per_record_limits.json').read_text())
        for case, name in data['measurement_files'].items():
            path = ROOT/'reports/results/filter_sim_gates'/name
            self.assertEqual(G.digest(path), data['measurement_sha256'][case])
            self.assertEqual(json.loads(path.read_text())['case'], case)
        keys = {(e['case'], e['family'], e['profile'], e['input']) for e in data['entries']}
        self.assertEqual(len(keys), 544)
        self.assertEqual(len(keys), len(data['entries']))
        for case in ('following', 'head'):
            self.assertEqual({e['family'] for e in data['entries'] if e['case'] == case},
                set(G.REPLAYS) - {'heave'} | {'heave_' + m for m in G.HEAVE_METHODS | {'pii'}})
        for entry in data['entries']:
            for metric, value in entry['limits'].items():
                self.assertEqual(value, cap(entry['measured_max'][metric]))
                self.assertGreater(value, entry['measured_max'][metric])
        self.assertEqual(render(data), (ROOT / 'tests/common/RecordReplayLimits.inc').read_text())
        horizontal = set(G.REPLAYS) - {'ahrs', 'heave', 'nlo', 'pii'}
        y_entries = [e for e in data['entries'] if 'displacement_y' in e['limits']]
        self.assertEqual(len(y_entries), 272)
        self.assertEqual({e['family'] for e in y_entries}, horizontal)
        self.assertTrue(all(('displacement_y' in e['limits']) == (e['family'] in horizontal)
                            for e in data['entries']))

    def test_compiled_record_lookup_rejects_hidden_regressions_and_missing_channels(self):
        # Exercise the actual C++ lookup, including path handling and failure modes.
        source = r'''
#include "tests/common/RecordW3dLimits.h"
#include "tests/common/CurrentReplayLimits.h"
#include <cassert>
using namespace ocean_imu::sim::record_limits;
namespace rl = ocean_imu::sim::record_limits;
int main() {
    const auto& e = rl::find("ou2", "default", "/tmp/wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv");
    const float global = 8.52521f;
    const float local = tighten(e, Metric::heave, global);
    const float regression = (local + global) / 2;
    assert(regression <= global && !(regression <= local));
    assert(tighten(e, Metric::heave, 0.01f) == 0.01f);
    const auto bars = ocean_imu::sim::current_limits::ou2;
    const auto capped = w3d("ou2", "default", e.input, bars);
    assert(capped.err_limit_y_rms_m == e.limits[static_cast<std::size_t>(Metric::displacement_y)]);
    const float y_regression = (capped.err_limit_y_rms_m + bars.err_limit_y_rms_m) / 2;
    assert(y_regression <= bars.err_limit_y_rms_m && y_regression > capped.err_limit_y_rms_m);
    auto diagnostic = bars;
    diagnostic.err_limit_y_rms_m = 0.0f;
    assert(w3d("ou2", "default", e.input, diagnostic).err_limit_y_rms_m == 0.0f);
    int caught = 0;
    for (const auto* unknown : {"missing.csv", "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv.bak"}) {
        try { (void)rl::find("ou2", "default", unknown); }
        catch (const std::invalid_argument&) { ++caught; }
    }
    try { (void)rl::find("ou2", "typo", e.input); }
    catch (const std::invalid_argument&) { ++caught; }
    try { (void)rl::find("typo", "default", e.input); }
    catch (const std::invalid_argument&) { ++caught; }
    try { (void)tighten(rl::find("nlo", "default", e.input), Metric::yaw, 10); }
    catch (const std::invalid_argument&) { ++caught; }
    assert(caught == 5);
    setenv("W3D_IMU_SEED", "1234", 1);
    assert(w3d("ou2", "default", e.input, bars).err_limit_y_rms_m == bars.err_limit_y_rms_m);
    assert(limit("ou3", "default", "diagnostic_crossfade.csv", Metric::heave, 8.0f) == 8.0f);
    unsetenv("W3D_IMU_SEED");
    const char* regular = "/tmp/wave_data_gerstner_H0.270_L14.047_A0.00_P60.00.csv";
    assert(uses_record_caps(regular));
    setenv("W3D_ALL_WAVE_TYPES", "1", 1);
    assert(!uses_record_caps(regular));
    assert(limit("nlo", "default", regular, Metric::heave, 8.19f) == 8.19f);
    assert(uses_record_caps(e.input));
    assert(uses_record_caps("wave_data_jonswap_typo.csv"));
    assert(limit("ou2", "default", e.input, Metric::heave, global) == e.limits[0]);
    setenv("W3D_ALL_WAVE_TYPES", "0", 1);
    assert(uses_record_caps(regular));
    unsetenv("W3D_ALL_WAVE_TYPES");
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp, binary = Path(tmp) / 'test.cpp', Path(tmp) / 'test'
            cpp.write_text(source)
            eigen = ROOT / 'third_party/eigen'
            if not (eigen / 'Eigen/Dense').exists():
                eigen = Path('/usr/include/eigen3')
            subprocess.run(['g++', '-std=c++20', '-Wall', '-Wextra', '-Werror', '-I', str(ROOT), '-I', str(ROOT / 'src'),
                            '-isystem', str(eigen),
                            str(cpp), '-o', str(binary)], check=True, capture_output=True)
            subprocess.run([str(binary)], check=True, capture_output=True)

    def test_heave_ablation_flags_do_not_use_default_record_bars(self):
        # Compile the simulator's actual protocol predicate without its Eigen
        # replay machinery; changes to the experiment must retain old bars.
        sim = (ROOT / 'tests/heave_only/heave_only-sim.cpp').read_text()
        enums = '\n'.join(re.findall(r'enum class (?:Method|Levelling) \{[^}]+\};', sim))
        helpers = sim[sim.index('bool env_flag('):sim.index('// Regression sentinels')]
        source = '#include <cstdlib>\n#include <string>\n#include <cassert>\n' + enums
        source += '\nLevelling g_frontend = Levelling::Mahony;\n' + helpers
        source += r'''
int main() {
    for (auto m : {Method::PII, Method::Godhavn, Method::GodhavnExt, Method::RichterZd, Method::Kuchler})
        assert(!diagnostic_configuration(m));
    for (const auto* flag : {"HO_KUCHLER_BODY_Z", "HO_KUCHLER_OFFSET_FIT"}) {
        setenv(flag, "1", 1);
        assert(diagnostic_configuration(Method::Kuchler));
        assert(!diagnostic_configuration(Method::Godhavn));
        setenv(flag, "0", 1);
        assert(!diagnostic_configuration(Method::Kuchler));
        unsetenv(flag);
    }
    for (const auto* value : {"HO_KUCHLER_ZETA_Q", "HO_KUCHLER_OMEGA_RW"}) {
        setenv(value, "0", 1);
        assert(diagnostic_configuration(Method::Kuchler));
        assert(!diagnostic_configuration(Method::PII));
        unsetenv(value);
    }
    setenv("HO_FIXED_WC", "0.2", 1);
    assert(diagnostic_configuration(Method::Godhavn));
    assert(diagnostic_configuration(Method::GodhavnExt));
    assert(diagnostic_configuration(Method::RichterZd));
    assert(!diagnostic_configuration(Method::Kuchler));
    unsetenv("HO_FIXED_WC");
    for (const auto* flag : {"HO_BIAS_TERM", "HO_SUBTRACT_MEAN"}) {
        setenv(flag, "1", 1);
        assert(diagnostic_configuration(Method::Godhavn));
        assert(!diagnostic_configuration(Method::GodhavnExt));
        assert(!diagnostic_configuration(Method::RichterZd));
        unsetenv(flag);
    }
    setenv("HO_MAHONY_TWOKP", "1", 1);
    assert(diagnostic_configuration(Method::PII));
    g_frontend = Levelling::Proxy;
    assert(diagnostic_configuration(Method::Kuchler));
    g_frontend = Levelling::Truth;
    assert(!diagnostic_configuration(Method::Kuchler));
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp, binary = Path(tmp) / 'test.cpp', Path(tmp) / 'test'
            cpp.write_text(source)
            subprocess.run(['g++', '-std=c++20', '-Wall', '-Wextra', '-Werror',
                            str(cpp), '-o', str(binary)], check=True, capture_output=True)
            subprocess.run([str(binary)], check=True, capture_output=True)

    def test_heave_requires_all_methods_inputs_and_the_selected_frontend(self):
        methods = sorted(G.HEAVE_METHODS | {'pii'})
        text = '\n'.join(f'HEAVE_ONLY method={method} frontend=mahony record={name} z_pct_hs=1'
                         for method in methods for name in G.REFERENCE_CSV_SHA256)
        self.assertEqual(len(G.parse_heave(text, 'default')), 40)
        invalid = (text.replace('frontend=mahony', 'frontend=proxy', 1),
                   text.replace('method=godhavn ', 'method=unknown ', 1),
                   text.replace('z_pct_hs=1', 'z_pct_hs=nan', 1),
                   '\n'.join(text.splitlines()[1:]),
                   text.replace(next(iter(G.REFERENCE_CSV_SHA256)),
                                list(G.REFERENCE_CSV_SHA256)[1], 1))
        for bad in invalid:
            with self.assertRaises(ValueError):
                G.parse_heave(bad, 'default')

    def test_pii_rejects_missing_duplicate_and_nonfinite_records(self):
        text = pii_fixture()
        self.assertEqual(len(G.parse_pii(text)), 8)
        for bad in (text.split('\n', 2)[-1], text + '\n' + text,
                    text.replace('yaw_rms_deg=0.4', 'yaw_rms_deg=nan', 1),
                    text.replace('samples=180000', 'samples=12000', 1)):
            with self.assertRaises(ValueError):
                G.parse_pii(bad)

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
        rows[-1]['disp_y_rms_m'] = 0.003
        rows[0]['disp_y_rms_m'] = 0.2
        worst = G.measured_worst(rows)
        self.assertEqual(worst['err_limit_percent_z_jonswap'], 7.0)
        self.assertEqual(worst['err_limit_percent_z_pmstokes'], 0.1)
        self.assertEqual(worst['bias_3d_percent'], 4.0)
        self.assertEqual(worst['gyro_bias_3d_percent'], 9.0)
        self.assertEqual(worst['err_limit_y_rms_m'], 0.2)

    def test_every_simulator_has_an_audit_owner(self):
        owned = {f'tests/{folder}/{binary}.cpp' for folder, binary, _ in G.REPLAYS.values()}
        owned.update(('tests/kalman_carma/carma_process-sim.cpp',
                      'tests/kalman_carma/carma_lever_arm-sim.cpp'))
        expected = {str(p.relative_to(ROOT)) for p in (ROOT / 'tests').rglob('*-sim.cpp')}
        expected.add('tests/pii_observer/pii_observer-adaptive.cpp')
        expected.add('tests/kalman_tangent_shared_modes/tangent-replay.cpp')
        self.assertEqual(owned, expected)

    def test_locked_limits_match_qualifying_receipts_and_current_scores_pass(self):
        from lock_simulation_baseline import make_baseline, render_global
        actual = json.loads((ROOT / 'reports/results/filter_sim_gates/per_record_limits.json').read_text())
        baseline = {c: json.loads((ROOT/'reports/results/filter_sim_gates'/name).read_text())
                    for c, name in actual['measurement_files'].items()}
        expected, global_caps = make_baseline(baseline)
        self.assertEqual(actual['entries'], expected['entries'])
        self.assertEqual(actual['inputs'], expected['inputs'])
        self.assertEqual(global_caps, json.loads((ROOT / 'reports/results/filter_sim_gates/limits.json').read_text()))
        self.assertEqual(render_global(global_caps), (ROOT / 'tests/common/CurrentReplayLimits.h').read_text())
        reports = {c: json.loads((ROOT/f'reports/results/sim_v1.4.0/{c}/measurements.json').read_text())
                   for c in ('following', 'head')}
        for report in reports.values():
            self.assertEqual(report['acceptance'], 'all replay gates passed')
            self.assertFalse(any(code for profiles in report['replay_exit_codes'].values()
                                 for code in profiles.values()))
        current, _ = make_baseline(reports)
        def key(row):
            return tuple(row[k] for k in ('case', 'family', 'profile', 'input'))
        caps = {key(row): row['limits'] for row in actual['entries']}
        self.assertEqual({key(row) for row in current['entries']}, set(caps))
        for row in current['entries']:
            for metric, score in row['measured_max'].items():
                self.assertLessEqual(score, caps[key(row)][metric], (key(row), metric))

    def test_make_log_keeps_sampling_profiles_separate_and_rejects_failed_runs(self):
        families = ('OU_II', 'OU_III', 'TFG', 'CARMA3_A41', 'CARMA_SHARED_MODES_CSM',
                    'TFG_CARMA3_T41', 'TFG_CARMA_SHARED_MODES_CSM', 'TFG_LEGACY', 'TANGENT_CARMA_SHARED_MODES_CSM')
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
        log = '\n'.join(blocks) + '\nAdaptiveVerticalPIIMahony simulation starting\n' + pii_fixture()
        parts = G.split_make_log(log)
        self.assertEqual(G.parse_w3d(parts['ou2', 'legacy'])[0]['disp_z_pct_hs'], 9.0)
        self.assertEqual(G.parse_w3d(parts['ou2', 'default'])[0]['disp_z_pct_hs'], 0.1)
        for invalid in (log + '\n' + blocks[0], log + '\nQUALITY_GATE: PASS=0',
                        log + '\nERROR: raw Z RMS above limit'):
            with self.assertRaises(ValueError):
                G.split_make_log(invalid)


if __name__ == '__main__':
    unittest.main()
