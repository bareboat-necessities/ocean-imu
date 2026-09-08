"""Attach observable Normal-Live premises to unchanged reference replays.

This is a non-promoting necessary-condition audit. It keeps out-of-domain
samples, actual accepted Jacobian vectors and every applied anisotropic R_S.
Passing sampled checks does not establish continuous BRMM, transported PE,
the covariance theorem, nonlinear storage attachment or prefix retention.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile

import numpy as np
import pandas as pd

import ou3_brmm_audit as SOURCE

ROOT = Path(__file__).resolve().parents[2]
DOMAIN = ROOT/"tools/stability/ou3_proof_operating_domain.json"
R_FIELDS = [f"R{i}{j}" for i in range(3) for j in range(3)]
ROUNDING = float(8*np.finfo(np.float32).eps)


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def runs(mask):
    """All maximal half-open sample-index intervals, including endpoint runs."""
    edges = np.diff(np.r_[False, np.asarray(mask, dtype=bool), False].astype(int))
    return np.column_stack((np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def summary(mask):
    intervals = runs(mask)
    return {"samples": int(np.count_nonzero(mask)), "intervals": len(intervals),
            "first_interval": intervals[0].tolist() if len(intervals) else None,
            "longest_interval_samples": int(np.max(intervals[:, 1]-intervals[:, 0])) if len(intervals) else 0}


def all_window_samples(mask, width):
    c = np.r_[0, np.cumsum(np.asarray(mask, dtype=np.int64))]
    return c[width:]-c[:-width] == width


def analyze(source, samples, events, domain):
    n = len(source)
    if n == 0 or len(samples) != n or not np.array_equal(samples['index'], np.arange(n)):
        raise ValueError("runtime/source sample index or length mismatch")
    if not np.all(np.isfinite(samples.to_numpy(dtype=float))):
        raise ValueError("nonfinite runtime sample")
    if not np.allclose(samples.dt, float(np.float32(.005)), rtol=1e-8, atol=0):
        raise ValueError("runtime is not the unchanged 200 Hz reference clock")
    if not np.array_equal(events.sequence, np.arange(len(events))):
        raise ValueError("missing/reordered runtime measurement event")
    if not np.all(np.isfinite(events.drop(columns='kind').to_numpy(dtype=float))):
        raise ValueError("nonfinite runtime measurement operand")
    event_indices = events['index'].to_numpy()
    if (np.any(event_indices != event_indices.astype(int)) or np.any(event_indices < 0)
            or np.any(event_indices >= n) or np.any(np.diff(event_indices) < 0)):
        raise ValueError("measurement event index is detached from sample clock")
    if not set(events.kind) <= {'acc', 'mag', 'S'}:
        raise ValueError("unknown measurement event")
    for kind, column in [('acc', 'acc_applied'), ('mag', 'mag_applied'), ('S', 's_applied')]:
        counts = np.bincount(events.loc[events.kind == kind, 'index'].to_numpy(dtype=int), minlength=n)
        if not np.array_equal(counts, samples[column]):
            raise ValueError(f"{kind} event capture count mismatch")
    if not np.array_equal(samples.resets, samples.acc_applied+samples.mag_applied+samples.s_applied):
        raise ValueError("Joseph/reset count mismatch")
    if np.any(samples.guard_calls != 1):
        raise ValueError("guard capture does not cover every IMU sample")

    live = (samples.outer_pre == 1) & (samples.outer_post == 1) & (samples.inner_pre == 1) & (samples.inner_post == 1)
    stable_mode = samples.bias_pre == samples.bias_post
    a = source[['acc_'+axis for axis in 'xyz']].to_numpy()
    g = source[['gyro_'+axis for axis in 'xyz']].to_numpy()
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(g)):
        raise ValueError("nonfinite physical source")
    an = SOURCE.sqrt(SOURCE.norm2(SOURCE.point(a)))
    gn = SOURCE.sqrt(SOURCE.norm2(SOURCE.point(g)))
    cfg, nl = domain['configured_runtime'], domain['normal_live']
    rate_cap = SOURCE.mul(SOURCE.point(nl['body_rate_norm_upper_deg_s']),
                          (SOURCE.down(np.pi/180), SOURCE.up(np.pi/180)))
    masks = {
        'outside_outer_and_inner_Live': ~live.to_numpy(),
        'bias_mode_transition': ~stable_mode.to_numpy(),
        'physical_acceleration_above_cap': an[0] > nl['non_gravitational_cog_acceleration_norm_upper_mps2'],
        'physical_body_rate_above_cap': gn[0] > rate_cap[1],
        'guard_not_dormant_transparent': ((samples.guard_identical != 1) | (samples.engagement != 0) | (samples.excess != 0)).to_numpy(),
        'prediction_or_accelerometer_missing': ((samples.predict != 1) | (samples.acc_calls != 1) | (samples.acc_applied != 1)).to_numpy(),
        'S_factorization_failure': (samples.s_calls != samples.s_applied).to_numpy(),
        'mag_factorization_failure': (samples.mag_calls != samples.mag_applied).to_numpy(),
        'hard_tilt_reset': (samples.tilt_resets != 0).to_numpy(),
        'magnetic_lock_or_refinement_transition': ((samples.mag_lock.diff().fillna(0) != 0) |
                                                  (samples.mag_refined.diff().fillna(0) != 0)).to_numpy(),
    }
    masks['source_cap_rounding_unresolved'] = (
        ((an[0] <= nl['non_gravitational_cog_acceleration_norm_upper_mps2']) &
         (an[1] > nl['non_gravitational_cog_acceleration_norm_upper_mps2'])) |
        ((gn[0] <= rate_cap[1]) & (gn[1] > rate_cap[0])))
    vector_report = {}
    for kind, sensor, lower, upper in [
        ('acc', 'accelerometer_mps2', nl['specific_force_norm_lower_mps2'], nl['specific_force_norm_upper_mps2']),
        ('mag', 'magnetometer_uT', nl['magnetic_vector_norm_lower_uT'], nl['magnetic_vector_norm_upper_uT'])]:
        e = events[events.kind == kind]
        index = e['index'].to_numpy(dtype=int)
        matrices = e[R_FIELDS].to_numpy().reshape(-1, 3, 3)
        expected = np.diag(np.asarray(cfg['measurement_noise_std'][sensor], dtype=float)**2)
        mismatch = np.any(np.abs(matrices-expected) > ROUNDING*np.maximum(np.abs(expected), 1e-30), axis=(1, 2))
        mask = np.zeros(n, dtype=bool)
        mask[index[mismatch]] = True
        masks[kind+'_configured_R_mismatch'] = mask
        norms = np.linalg.norm(e[['bx', 'by', 'bz']].to_numpy(), axis=1)
        mask = np.zeros(n, dtype=bool)
        mask[index[(norms < lower) | (norms > upper)]] = True
        masks[kind+'_Jacobian_vector_norm_outside_bounds'] = mask
        vector_report[kind] = {
            'accepted_events': len(e), 'configured_R': expected.tolist(),
            'R_diagonal_min': np.diagonal(matrices, axis1=1, axis2=2).min(axis=0).tolist() if len(e) else None,
            'R_diagonal_max': np.diagonal(matrices, axis1=1, axis2=2).max(axis=0).tolist() if len(e) else None,
            'norm_min': float(norms.min()) if len(e) else None,
            'norm_max': float(norms.max()) if len(e) else None,
            'comparison_rounding_relative_tolerance': ROUNDING}
    # Necessary packet recurrence only. These diagnostics do not assert the
    # transported two-vector matrix lower bound from the P3 premise manifest.
    width = round(nl['vector_pe_recurrence_window_s']*200)
    mag_windows = all_window_samples(samples.mag_applied == 0, width)
    valid_live_windows = all_window_samples(live & stable_mode, width)
    paired = events[events.kind == 'acc'][['index', 'wx', 'wy', 'wz']].merge(
        events[events.kind == 'mag'][['index', 'wx', 'wy', 'wz']], on='index', suffixes=('_a', '_m'))
    va = paired[[x+'_a' for x in ['wx', 'wy', 'wz']]].to_numpy()
    vm = paired[[x+'_m' for x in ['wx', 'wy', 'wz']]].to_numpy()
    denominator = np.linalg.norm(va, axis=1)*np.linalg.norm(vm, axis=1)
    sine = np.divide(np.linalg.norm(np.cross(va, vm), axis=1), denominator,
                     out=np.zeros(len(paired)), where=denominator > 0)
    small_sine = np.zeros(n, dtype=bool)
    small_sine[paired['index'].to_numpy(dtype=int)[sine < nl['vector_sine_separation_lower']]] = True
    # The body and world vectors come from the accepted Jacobian, not the
    # measured noisy vector. Their common-frame sine is still a point diagnostic.
    regimes = {'all': np.ones(n, dtype=bool), 'outer_inner_Live': live.to_numpy(),
               'H18': (live & stable_mode & (samples.bias_post == 0)).to_numpy(),
               'A21': (live & stable_mode & (samples.bias_post == 1)).to_numpy()}
    violations = {name: {mode: summary(mask & select) for mode, select in regimes.items()}
                  for name, mask in masks.items()}
    physical_ok = ~(masks['physical_acceleration_above_cap'] | masks['physical_body_rate_above_cap'] |
                    masks['source_cap_rounding_unresolved'])
    necessary = ~np.logical_or.reduce(list(masks.values()))
    # Phase partition depends only on the runtime flags, never on a cap.
    # acc_xyz is geographic Z-up CoG motion, without gravity or sensor bias.
    force = SOURCE.add(SOURCE.point(a), SOURCE.point(np.array([0., 0., domain['startup']['gravity_mps2']])))
    fn = SOURCE.sqrt(SOURCE.norm2(force))
    phase_extrema = {}
    entered = np.maximum.accumulate(live.to_numpy())
    for phase, select in [('before_Live', ~entered), ('Live', live.to_numpy()),
                          ('non_Live_after_entry', entered & ~live.to_numpy())]:
        phase_extrema[phase] = {
            'samples': int(select.sum()),
            'acceleration_norm_max_upper_mps2': float(an[1][select].max()) if select.any() else None,
            'body_rate_norm_max_upper_rad_s': float(gn[1][select].max()) if select.any() else None,
            'specific_force_norm_min_lower_mps2': float(fn[0][select].min()) if select.any() else None,
            'specific_force_norm_max_upper_mps2': float(fn[1][select].max()) if select.any() else None,
            'vertical_acceleration_min_mps2': float(a[select, 2].min()) if select.any() else None,
            'vertical_acceleration_max_mps2': float(a[select, 2].max()) if select.any() else None,
            'acceleration_norm_above_g_samples': int((an[0][select] > domain['startup']['gravity_mps2']).sum())}
    # Continuous reference adaptation is recorded but is not automatically a
    # forbidden hard regauge: the domain explicitly retains that variation.
    observations = {'magnetic_reference_write': summary(samples.reference_writes != 0),
                    'same_sample_vector_sine_below_bound': summary(small_sine)}
    impulse_report = {}
    if n >= 2000:
        primitive = SOURCE.cumulative(SOURCE.mul(SOURCE.point(a), (SOURCE.down(.005), SOURCE.up(.005))))
        impulse = SOURCE.sqrt(SOURCE.norm2(SOURCE.intervals(primitive, 0, 2000, n-1999)))
        for mode, select in regimes.items():
            whole_window = all_window_samples(select, 2000)
            impulse_report[mode] = {
                'complete_10s_windows': int(whole_window.sum()),
                'candidate_2_mps_definite_violations': summary(whole_window & (impulse[0] > 2.)),
                'max_upper_mps': float(impulse[1][whole_window].max()) if whole_window.any() else None}
        # These intervals index window STARTS, not individual violating samples.
        masks['candidate_10s_impulse_above_2_mps_window_starts'] = impulse[0] > 2.
    se = events[events.kind == 'S']
    return {
        'samples': n, 'regimes': {k: summary(v) for k, v in regimes.items()},
        'violations': violations, 'vectors': vector_report,
        'physical_phase_extrema': phase_extrema,
        'observations': observations, 'candidate_impulse_windows': impulse_report,
        'physical_caps_inside_runtime_Live': summary(physical_ok & live.to_numpy()),
        'passes_observed_sample_checks_only': summary(necessary),
        'live_1s_windows_without_accepted_mag': int((mag_windows & valid_live_windows).sum()),
        'same_sample_accepted_pair_count': len(paired),
        'same_sample_common_frame_sine_min_diagnostic': float(sine.min()) if len(sine) else None,
        'S_calls': int(samples.s_calls.sum()), 'S_applied': len(se),
        'S_R_diagonal_min': se[['R00','R11','R22']].min().tolist() if len(se) else None,
        'S_R_diagonal_max': se[['R00','R11','R22']].max().tolist() if len(se) else None,
        'source_samples_pruned': 0,
        'runtime_Live_implies_Normal_Live': False,
        'continuous_BRMM_admission_certified': False,
        'transported_vector_PE_certified': False,
        'P3_certified_by_this_audit': False,
        'P4_or_P5_promoted': False,
        'remaining_admission_obligations': [
            'continuous no-DC/recurrence and source primitives',
            'transported asynchronous vector PE on the complete word',
            'source-uniform entry, Q/floor/tuner invariant attachment',
            'nonlinear projection and every-prefix chart/domain retention'],
    }, masks


def compact_report(report):
    """Small checked-in summary; the artifact retains all interval statistics."""
    result = {k: v for k, v in report.items() if k not in ('cases', 'overlay_manifest')}
    result['shipping_source_sha256'] = {
        p: v['shipping_sha256'] for p, v in report['overlay_manifest']['files'].items()}
    result['cases'] = []
    for case in report['cases']:
        reduced = {k: v for k, v in case.items() if k != 'violations'}
        reduced['violation_sample_counts'] = {
            reason: {mode: stat['samples'] for mode, stat in modes.items()}
            for reason, modes in case['violations'].items()}
        result['cases'].append(reduced)
    return result


def verify_manifest(manifest, repository):
    if not manifest.get('stripping_recovers_shipping_bytes'):
        raise ValueError('unverified runtime overlay')
    for path, data in manifest['files'].items():
        if digest(repository/path) != data['shipping_sha256']:
            raise ValueError('runtime overlay is detached from shipping source: '+path)


def replay(binary, directory, filename, env):
    command = [str(binary), '--input', filename]
    with (directory/'stdout.txt').open('wb') as out, (directory/'stderr.txt').open('wb') as err:
        result = subprocess.run(command, cwd=directory, env=env, stdout=out, stderr=err)
    stdout = (directory/'stdout.txt').read_text()
    gate_failure = result.returncode == 1 and 'QUALITY_GATE: PASS=0' in stdout
    if result.returncode and not (env.get('W3D_COLLECT_ALL_GATES') == '1' and gate_failure):
        raise ValueError(f'simulator failed ({result.returncode}): ' +
                         (directory/'stderr.txt').read_text()[-4000:])
    outputs = sorted(p for p in directory.iterdir() if p.name.startswith('w3d_') or p.suffix == '.txt')
    if not any(p.name.startswith('w3d_') and p.suffix == '.csv' for p in outputs):
        raise ValueError('ordinary simulator output missing')
    return {'exit_code': result.returncode, 'files': {p.name: digest(p) for p in outputs}}


def run(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    domain = json.loads(DOMAIN.read_text())
    manifest = json.loads(args.manifest.read_text())
    verify_manifest(manifest, ROOT)
    report = {'qualification': 'OU3_BRMM_RUNTIME_NECESSARY_CONDITIONS_V1',
              'source_commit': subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
              'archive_sha256': digest(args.archive), 'domain_sha256': digest(DOMAIN),
              'baseline_binary_sha256': digest(args.baseline), 'traced_binary_sha256': digest(args.traced),
              'overlay_manifest': manifest,
              'audit_producer_sha256': digest(Path(__file__)),
              'observer_header_sha256': digest(Path(__file__).with_name('ou3_brmm_runtime_trace.h')),
              'overlay_producer_sha256': digest(Path(__file__).with_name('ou3_brmm_runtime_overlay.py')),
              'clock': '200 Hz sample index; snapshot follows IMU then asynchronous 25 Hz magnetometer',
              'protocol': 'unchanged simulator, default noise/seeds/magnetometer, all original quality gates',
              'observer_build': 'baseline and trace use identical flags with -ffp-contract=off; native FMA capture is not bitwise transparent',
              'scope': 'sampled necessary-condition audit; no physical admission or P3/P4/P5 promotion',
              'cases': []}
    # No inherited sweep/ablation or scoring overrides may silently alter the
    # reference experiment. Ordinary PATH/compiler/runtime environment survives.
    prefixes = ('OU_', 'OU3_', 'OU_III_', 'SF_', 'W3D_')
    env = {k: v for k,v in os.environ.items() if not k.startswith(prefixes)}
    env['W3D_COLLECT_ALL_GATES'] = '1'
    with zipfile.ZipFile(args.archive) as archive:
        for family in SOURCE.FAMILIES:
            for height in SOURCE.HEIGHTS:
                pattern = rf'wave_data_{family}_H{height:.3f}_L[^_]+_A[^_]+_P[^_]+\.csv'
                names = [n for n in archive.namelist() if re.fullmatch(pattern, Path(n).name)]
                if len(names) != 1:
                    raise ValueError('reference case is missing or ambiguous')
                name = Path(names[0]).name
                case_dir = args.output_dir/Path(name).stem
                case_dir.mkdir(exist_ok=True)
                with tempfile.TemporaryDirectory(prefix='ou3-brmm-') as temp:
                    temp = Path(temp)
                    print('BRMM_RUNTIME_BEGIN', name, flush=True)
                    for label in ('baseline', 'traced'):
                        directory = temp/label
                        directory.mkdir()
                        (directory/name).write_bytes(archive.read(names[0]))
                    baseline = replay(args.baseline, temp/'baseline', name, env)
                    trace_prefix = case_dir/'runtime'
                    traced = replay(args.traced, temp/'traced', name, dict(env, OU3_BRMM_TRACE=str(trace_prefix)))
                    if baseline != traced:
                        for label in ('baseline', 'traced'):
                            for path in (temp/label).glob('*.txt'):
                                shutil.copyfile(path, case_dir/(label+'-'+path.name))
                        (case_dir/'transparency-failure.json').write_text(json.dumps(
                            {'baseline':baseline, 'traced':traced}, indent=2)+'\n')
                        raise ValueError('runtime overlay changed ordinary simulator output: '+name)
                    source_path = temp/'baseline'/name
                    samples_path = Path(str(trace_prefix)+'.samples.csv')
                    events_path = Path(str(trace_prefix)+'.events.csv')
                    result, masks = analyze(pd.read_csv(source_path), pd.read_csv(samples_path), pd.read_csv(events_path), domain)
                    result.update(record=name, csv_sha256=digest(source_path),
                                  ordinary_output_sha256=baseline["files"],
                                  regression_exit_code=baseline["exit_code"],
                                  regression_quality_gates_pass=baseline["exit_code"] == 0,
                                  read_only_runtime_output_bitwise_identical=True)
                    with gzip.open(case_dir/'violating-intervals.csv.gz', 'wt') as stream:
                        stream.write('reason,start_sample,stop_sample_exclusive\n')
                        for reason, mask in masks.items():
                            for start, stop in runs(mask):
                                stream.write(f'{reason},{start},{stop}\n')
                    result['trace_sha256'] = {}
                    for path in (samples_path, events_path):
                        result['trace_sha256'][path.name] = digest(path)
                        with path.open('rb') as src, gzip.open(str(path)+'.gz', 'wb') as dst:
                            shutil.copyfileobj(src, dst)
                        with gzip.open(str(path)+'.gz', 'rb') as check:
                            if hashlib.file_digest(check, 'sha256').hexdigest() != result['trace_sha256'][path.name]:
                                raise ValueError('compressed runtime trace failed roundtrip integrity')
                        path.unlink()
                    for filename in ('stdout.txt', 'stderr.txt'):
                        shutil.copyfile(temp/'baseline'/filename, case_dir/filename)
                report['cases'].append(result)
                (args.output_dir/'runtime-audit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
                (args.output_dir/'runtime-summary.json').write_text(
                    json.dumps(compact_report(report), indent=2, allow_nan=False)+'\n')
                print('BRMM_RUNTIME_RESULT', name, json.dumps({k:result[k] for k in (
                    'samples','S_applied','physical_caps_inside_runtime_Live','passes_observed_sample_checks_only')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--traced', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    for key in ('archive','baseline','traced','manifest','output_dir'):
        setattr(args, key, getattr(args,key).resolve())
    run(args)
