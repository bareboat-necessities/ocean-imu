#!/usr/bin/env python3
"""Finite native history feasibility, not a nonlinear source-uniform certificate.

Quiet physical BIAS0 source, a predeclared gyro-disturbance prefix, followed by
zero sensor disturbances. Host rounding remains part of every measured word.
No state/covariance installation, metric fitting, or eigenvalue flooring.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import mpmath as mp
import numpy as np
from tools.stability.ou3_alt_contraction import shipping_finite_identity as NATIVE
from tools.stability.ou3_alt_contraction.binary32_covariance_coercivity import exact_spd

SOURCE = NATIVE.ROOT / 'tests/ou3_alt_contraction/finite_error_storage_probe.cpp'


def nonlinear_budget(tangent_norm, remainder_gain, input_metric_loss,
                     output_metric_gain, young, supply_gain):
    """Exact conditional comparison coefficients; no native premise is inferred.

    See ou3-alt-finite-error-storage.md for the quantified norm hypotheses.
    Inputs are exact rationals; accepting floats would hide budget rounding.
    """
    values = (tangent_norm, remainder_gain, input_metric_loss,
              output_metric_gain, young, supply_gain)
    if any(not isinstance(x, (int, Fraction)) for x in values):
        raise TypeError('exact rational budget required')
    a, eps, di, do, eta, b = map(Fraction, values)
    if min(a, eps, di, do, b) < 0 or di >= 1 or eta <= 0:
        raise ValueError('invalid metric/remainder budget')
    rho = (1+do)/(1-di)*(1+eta)*(a+eps)**2
    supply = (1+do)*(1+1/eta)*b*b
    return rho, supply


def analyze(path):
    data = np.loadtxt(path)
    if data.shape != (1801, 464) or not np.array_equal(data[:, 0], np.arange(1801)):
        raise ValueError('complete 1800-step native trace required')
    errors = data[:, 2:23]
    covariances = data[:, 23:].reshape(-1, 21, 21)
    if not np.isfinite(data).all() or not np.array_equal(covariances, covariances.transpose(0, 2, 1)):
        raise ValueError('finite exactly symmetric covariance required')
    np.linalg.cholesky(covariances)
    energy = np.einsum('bi,bi->b', errors, np.linalg.solve(covariances, errors[..., None])[..., 0])
    endpoints = []
    with mp.workdps(60):
        for index in (0, 600, 1200, 1800):
            P = covariances[index]
            audit = exact_spd(P.tolist())
            # Native float covariance/error values roundtrip through binary64.
            # mp.mpf(float) imports that binary64 value exactly at this precision.
            p = mp.matrix(P.tolist())
            e = mp.matrix(errors[index].tolist())
            v = (e.T * mp.lu_solve(p, e))[0]
            endpoints.append({'sample': index, 'V_60digits': mp.nstr(v, 60),
                              'exact_covariance_SPD': audit['exact_SPD'],
                              'minimum_exact_pivot': str(audit['minimum_exact_LDL_pivot']),
                              'active_bias': bool(data[index, 1]),
                              'error_joint24': errors[index].tolist() + [0.0]*3})
        words = []
        for j in range(3):
            v0 = mp.mpf(endpoints[j]['V_60digits'])
            v1 = mp.mpf(endpoints[j+1]['V_60digits'])
            if v0 <= 0:
                raise ValueError('nonzero finite initial storage required')
            words.append({'start_sample': 600*j, 'ratio_60digits': mp.nstr(v1/v0, 60),
                          'prefix_max_over_start': float(max(energy[600*j:600*(j+1)+1])/float(v0))})
    return {'endpoints': endpoints, 'words': words,
            'native_trace_sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'sample_boundary_float_cholesky_checks': len(data)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = work / 'finite-error-probe'
    command = ['g++', '-std=c++20', '-O1', '-ffp-contract=off', '-fno-fast-math',
               '-DEIGEN_DONT_VECTORIZE', '-I'+str(NATIVE.eigen_include()),
               '-I'+str(NATIVE.ROOT/'src'), str(SOURCE), '-o', str(binary)]
    subprocess.run(command, check=True)
    binary.chmod(binary.stat().st_mode | 0o111)
    rows = []
    for mode in ('H', 'A', 'HA'):
        for axis in range(3):
            for amplitude in ('0.001', '0.01'):
                path = work / f'{mode}-{axis}-{amplitude}.txt'
                metadata = json.loads(subprocess.check_output(
                    [str(binary), mode, str(axis), amplitude, str(path)], text=True))
                row = {'mode': mode, 'axis': axis, 'requested_pulse_rad_s': amplitude,
                       'native': metadata, **analyze(path)}
                rows.append(row)
                print(mode, axis, amplitude, max(float(w['ratio_60digits']) for w in row['words']), flush=True)
    report = {'qualification': 'OU3_ALT_FINITE_ERROR_STORAGE_DIAGNOSTIC_V1',
              'source': 'stationary identity attitude, zero physical BIAS0, horizontal 32uT field',
              'sensor_prefix': '100 samples of single-axis gyro disturbance; then zero',
              'magnetic_calls': 'every eight global samples; acceptance is not inferred',
              'dt_binary32_s': float(np.float32(.005)), 'word_samples': 600,
              'probe_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'shipping_tree_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD:src'], cwd=NATIVE.ROOT, text=True).strip(),
              'compiler': subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
              'compile_flags': command[1:6],
              'rows': rows,
              'worst_measured_ratio': max(float(w['ratio_60digits']) for r in rows for w in r['words']),
              'finite_nonlinear_source_uniform_rho_proved': False,
              'infinite_SPD_preservation_proved': False,
              'P4_PASS': False, 'P5_MAY_START': False,
              'limitations': ['Finite quiet BIAS0 histories only; not BRMM/BIAS1/2 coverage.',
                              '60-digit evaluation is not an outward-rounded enclosure; Cayley coordinates are serialized binary64 quotients.',
                              'Host execution does not qualify target FCR/compiler/libm.',
                              'Zero sensor supply after pulse does not remove machine roundoff.',
                              'Prefix checks are at IMU boundaries, not every internal operation.',
                              'No basin, universal capture, or informative acceptance theorem.']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')

if __name__ == '__main__':
    main()
