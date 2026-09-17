"""Compare native magnetic service gaps and literal ungauged quotient maps."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from tools.stability.ou3_alt_contraction import carried_storage_rho_diagnostic as D
from tools.stability.ou3_alt_contraction import magnetic_service_formulation as F
from tools.stability.ou3_alt_contraction.informative_service_theorem import InformativeServiceAssumption


def quotient_metric(M, N, Q):
    # inf_n (Q y + N n)^T M (Q y + N n); do not merely delete M rows.
    cross = Q.T @ M @ N
    result = Q.T @ M @ Q - cross @ np.linalg.solve(N.T @ M @ N, cross.T)
    return (result+result.T)/2


def quotient_diagnostic(A, M0=None, M1=None):
    A = np.asarray(A)
    rows = {}
    for label, removed in (('heading_only', (2,)), ('heading_and_axial_bias', (2, 5))):
        keep = [i for i in range(24) if i not in removed]
        N, Q = np.eye(24)[:, removed], np.eye(24)[:, keep]
        leakage = float(np.linalg.norm(Q.T @ A @ N, 2))
        row = {'fibre_to_quotient_leakage': leakage,
               'physical_equivariance_certified': False,
               'axial_bias_is_heading_gauge': False}
        try:
            quotient = F.quotient_map(A, N, N, Q, Q)
        except ValueError:
            row['point_quotient_exists'] = False
        else:
            Aq = quotient['A_quotient']
            motion = [j for j, i in enumerate(keep) if i < 18]
            Z = np.eye(len(keep))[:, motion]
            row.update(point_quotient_exists=True,
                       motion_spectral_radius=float(max(abs(np.linalg.eigvals(Aq[np.ix_(motion, motion)])))),
                       identity_storage_ratio=F.projected_storage_ratio(Aq, np.eye(len(keep)), np.eye(len(keep)), Z))
        if row['point_quotient_exists'] and M0 is not None and M1 is not None:
            row['compatible_quotient_storage_ratio'] = F.projected_storage_ratio(
                Aq, quotient_metric(M0, N, Q), quotient_metric(M1, N, Q), Z)
        rows[label] = row
    return rows


def analyze_history(path, metadata):
    windows = {}
    for samples in (600, 1200):
        assumption = InformativeServiceAssumption(1.0, samples*.005, .01, 1.0)
        rows = D.analyze(path, metadata, window_samples=samples, include_map=True)
        for row in rows:
            matrix = row.pop('joint24_tangent_map')
            M0 = D.covariance_metric(row.pop('covariance_before'))
            M1 = D.covariance_metric(row.pop('covariance_after'))
            if metadata['mag_stride'] == 0:
                row['ungauged_quotients'] = quotient_diagnostic(matrix, M0, M1)
            row['conditional_service_audit'] = assumption.audit_native_window(row, metadata['dt_s'])
        windows[str(samples)] = {'assumption': assumption.declaration(), 'rows': rows}
    return windows


def run(work):
    D.PLAN.require_theorem_task(
        obligation='conditional service and quotient complete-word feasibility',
        evidence_kind='analytic_stationary_source_native_tangent', complete_physical_word=False,
        requested_phase='feasibility_diagnostic')
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    manifest = D.NATIVE.make_overlay(work/'include')
    binaries = {}
    for observed in (True, False):
        binary = work/('observed' if observed else 'plain')
        cmd = [os.environ.get('CXX', 'g++'), '-std=c++20', '-O1', '-ffp-contract=off',
               '-fno-fast-math', '-DEIGEN_DONT_VECTORIZE']
        if observed:
            cmd += ['-I'+str(work/'include')]
        cmd += ['-I'+str(D.NATIVE.eigen_include()), '-I'+str(D.NATIVE.ROOT/'src'),
                str(D.SOURCE), '-o', str(binary)]
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=240)
        binary.chmod(binary.stat().st_mode | 0o111)
        binaries[observed] = binary
    histories = {}
    # Predeclared experiment grid. These are callback schedules, never claimed
    # to guarantee acceptance; actual accepted H/R/event ancestry is audited.
    for mode in ('H', 'A', 'HA'):
        for stride in (0, 8, 40, 200):
            key = f'{mode}:stride{stride}'
            metadata = None
            paths = {}
            for observed in (True, False):
                path = work/(key.replace(':', '-')+('-observed' if observed else '-plain'))
                proc = subprocess.run([str(binaries[observed]), mode, '1200',
                                       str(path)+'.txt', str(path)+'.bin', str(stride)],
                                      check=True, capture_output=True, text=True, timeout=180)
                current = json.loads(proc.stdout)
                if metadata is not None and current != metadata:
                    raise ValueError('instrumentation changed metadata')
                metadata = current
                paths[observed] = path
            if not D.NATIVE.same_files(Path(str(paths[True])+'.bin'), Path(str(paths[False])+'.bin')):
                raise ValueError('passive instrumentation changed states')
            windows = analyze_history(str(paths[True])+'.txt', metadata)
            histories[key] = {'native': metadata, 'windows': windows}
            print(key+' complete', flush=True)
    report = {'qualification': 'OU3_ALT_NATIVE_SERVICE_REGIME_DIAGNOSTIC_V1',
            'histories': histories,
            'source_manifest': manifest,
            'probe_sha256': hashlib.sha256(D.SOURCE.read_bytes()).hexdigest(),
            'passive_instrumentation_bit_identity': True,
            'source': 'stationary zero BIAS0 horizontal field; real startup at 25 Hz, then declared schedule',
            'metric_law': 'diag(inverse(full_carried_P21), I3)',
            'tangent_only': True, 'metric_fitted': False,
            'source_uniform_rho_certified': False, 'uniform_metric_coercivity_certified': False,
            'infinite_service_recurrence_certified': False,
            'storage_search_allowed': False, 'ALT_LIVE_PASS': False,
            'ALT_STARTUP_PASS': False, 'ALT_END_TO_END_PASS': False}

    annotate_schedule_scope(report)
    failures = validate(report)
    if failures:
        raise ValueError('native service diagnostic invalid: '+repr(failures))
    return report


def annotate_schedule_scope(report):
    """Keep exploratory callback outages outside the retained 40 ms profile."""
    for history in report['histories'].values():
        stride = history['native']['mag_stride']
        retained = 0 < stride <= 8
        history['retained_MAG_CALL_SCHEDULE_v1_satisfied'] = retained
        history['schedule_scope'] = ('retained_runtime_profile' if retained else
                                     'outage_stress_outside_retained_callback_profile')
        for group in history['windows'].values():
            for row in group['rows']:
                row['full_declared_profile_point_pass'] = bool(retained and
                    row['conditional_service_audit']['finite_window_membership_point_pass'])


def validate(report):
    failures = []
    for key in ('source_uniform_rho_certified', 'uniform_metric_coercivity_certified',
                'infinite_service_recurrence_certified', 'storage_search_allowed',
                'ALT_LIVE_PASS', 'ALT_STARTUP_PASS', 'ALT_END_TO_END_PASS', 'metric_fitted'):
        if report.get(key) is not False:
            failures.append(key+' incorrectly promoted')
    if report.get('passive_instrumentation_bit_identity') is not True:
        failures.append('native parity absent')
    if report.get('tangent_only') is not True:
        failures.append('tangent scope absent')
    histories = report.get('histories', {})
    if set(histories) != {f'{m}:stride{s}' for m in ('H', 'A', 'HA') for s in (0, 8, 40, 200)}:
        failures.append('predeclared family incomplete')
    for key, history in histories.items():
        if history.get('retained_MAG_CALL_SCHEDULE_v1_satisfied') != (0 < history['native']['mag_stride'] <= 8):
            failures.append(key+' callback profile misclassified')
        if set(history['windows']) != {'600', '1200'}:
            failures.append(key+' horizons incomplete')
        for size, group in history['windows'].items():
            samples = int(size)
            if len(group['rows']) != 1200//samples:
                failures.append(key+' window count incorrect')
            for i, row in enumerate(group['rows']):
                if (row['start_sample'], row['end_sample']) != (i*samples, (i+1)*samples):
                    failures.append(key+' detached window')
                mode = history['native']['mode']
                expected_mode = 'H' if mode == 'H' or (mode == 'HA' and row['end_sample'] <= 600) else 'A'
                if row['mode_at_end'] != expected_mode:
                    failures.append(key+' native mode mismatch')
                expected_release = int(mode == 'HA' and row['start_sample'] < 601 <= row['end_sample'])
                if row['release_edges'] != expected_release:
                    failures.append(key+' release ancestry mismatch')
                if any(row['event_counts'].get(kind) != samples for kind in ('2', '12')):
                    failures.append(key+' incomplete prediction/acceleration word')
                if row['double_vs_extended_terminal_difference'] > 1e-8 or row['energy_accounting_residual'] > 1e-8:
                    failures.append(key+' precision/accounting mismatch')
                if not np.isfinite(row['storage_ratio']['rho_point']):
                    failures.append(key+' nonfinite rho')
                if row['conditional_service_audit']['source_uniform_service_certified'] is not False:
                    failures.append(key+' service promoted')
    # A ratio above one or failed point membership is evidence, not a schema
    # failure: preserve it for the mandatory mathematical critic pass.
    return failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.work)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
