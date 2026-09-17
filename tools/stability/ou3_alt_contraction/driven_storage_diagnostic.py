#!/usr/bin/env python3
"""Non-promoting forced nonlinear words with analytically bounded source supply.

Reported C requirements are necessary lower bounds from finite histories, never
fitted theorem constants. The declared supply uses fixed SI coordinate scales.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import mpmath as mp
import numpy as np
from tools.stability.ou3_alt_contraction import shipping_finite_identity as NATIVE
from tools.stability.ou3_alt_contraction.binary32_covariance_coercivity import exact_spd
import ou3_p4_bias0_family as B0
import ou3_p4_bias1_family as B1
import ou3_p4_bias2_family as B2
import ou3_brmm_physical_wave_source as WAVE
import ou3_brmm_physical_wave_condition as PHYS

SOURCE = NATIVE.ROOT/'tests/ou3_alt_contraction/driven_storage_probe.cpp'
RHO_TARGETS = ('0.95', '0.98', '0.99')


def source_certificate(omega, family):
    """Exact all-time member bounds, independent of any replay or error state."""
    w = F(omega)
    if w not in (F(1,2), F(1)) or family not in (0, 1, 2):
        raise ValueError('predeclared source member required')
    contracts = [m.build() for m in (B0, B1, B2)]
    for module, contract in zip((B0, B1, B2), contracts):
        if module.validate(contract):
            raise ValueError('upstream bias contract invalid')
    contract = contracts[family]
    domain = json.loads((NATIVE.ROOT/'tools/stability/ou3_proof_operating_domain.json').read_text())['complete_brmm_physical_envelope']
    caps = {
        'wave_position_norm_upper_m': F(1,4),
        'wave_velocity_norm_upper_mps': w/4,
        'wave_acceleration_norm_upper_mps2': w*w/4,
        'body_rate_norm_upper_deg_s': F(0),
        'centered_primitive_D_S_upper_m_s': F(1,2)/w,
        # Hs=4*std(z)<=4*|z|<=4*|p| is sufficient here.
        'significant_wave_height_Hs_upper_m': F(1),
    }
    flo,fhi = map(lambda x: F(str(x)), domain['frequency_support_hz'])
    if any(bound > F(str(domain[key])) for key,bound in caps.items()) or not flo <= w*F(7,44) <= w/6 <= fhi:
        raise ValueError('analytic wave member exceeds COMPLETE-BRMM envelope')
    # 3 < pi < 22/7 bounds the exact sinusoid derivative and frequency.
    wb_upper = F(44, 7*600)
    h = F(1,200)  # native binary32 h is smaller
    if family == 0:
        if not contract['gauss_markov_tau_true_s'][0] <= 600 <= contract['gauss_markov_tau_true_s'][1]:
            raise ValueError('BIAS0 persistent root outside family')
        magnitude = F(25,1000)
        rate = F(5,1000)*wb_upper
        driver = h*(rate+magnitude/F(600))
        if F(str(contract['turn_on_offset_component_abs_upper_mps2'])) < F(2,100) or F(str(contract['non_gauss_markov_component_abs_upper_mps2'])) < F(5,1000) or F(str(contract['non_gauss_markov_rate_abs_upper_mps3'])) < rate:
            raise ValueError('BIAS0 component envelope changed')
    elif family == 1:
        magnitude = F(95,1000)
        rate = F(8,100)/1200 + F(15,1000)*wb_upper
        driver = h*F(15,1000)*(wb_upper+F(1,1200))
        if not (contract['tau_true_s'][0] <= 1200 <= contract['tau_true_s'][1] and contract['sinusoid_period_s'][0] <= 600 <= contract['sinusoid_period_s'][1] and F(str(contract['root_component_abs_upper_mps2'])) >= F(8,100) and F(str(contract['sinusoid_component_amplitude_abs_upper_mps2'])) >= F(15,1000)):
            raise ValueError('BIAS1 parameter envelope changed')
    else:
        magnitude = F(6,100)
        rate = F(1,100)*wb_upper
        driver = h*rate
        if not contract['non_relaxing_limit_admitted'] or F(str(contract['variation_rate_abs_upper_mps3'])) < rate:
            raise ValueError('BIAS2 drift envelope changed')
    if magnitude > F(str(contract['true_bias_component_abs_upper_mps2'])) or driver > F(str(contract['driver_increment_component_abs_upper_mps2'])):
        raise ValueError('member exceeds bias contract')
    wave = WAVE.spectral_certificate([WAVE.SpectralBand(w,w,F(1,4))])
    physical = PHYS.qualify_with_certificate(wave)
    # xi=(p,v,S_centered,a,beta), each divided by one in its declared SI unit.
    # ||direction||=1; S uses its one actual Live origin, so |S|<=2A/w.
    bounds = (F(1,4), w/4, F(1,2)/w, w*w/4, magnitude)
    D2 = sum((x*x for x in bounds), F(0))
    return {'family': f'BIAS{family}', 'omega_rad_s': str(w),
            'source_qualification': physical,
            'analytic_wave_bounds_within_COMPLETE_BRMM': True,
            'frequency_hz_enclosure': [str(w*F(7,44)), str(w/6)],
            'bias_magnitude_bound': str(magnitude), 'bias_derivative_bound': str(rate),
            'bias_driver_bound': str(driver),
            'bias_root_tau_s': (600 if family == 0 else 1200 if family == 1 else None),
            'bias2_phi_one': family == 2,
            'physical_supply_coordinates': ['p', 'v', 'S_centered', 'a', 'beta'],
            'physical_supply_squared_bound': str(D2),
            'machine_arithmetic_supply_bound_proved': False}


def necessary_supply(v0, v1, rho, D2):
    """C >= max(0,V1-rho*V0)/D2 for this endpoint pair and declared D."""
    if v0 < 0 or v1 < 0 or not 0 <= rho < 1 or D2 <= 0:
        raise ValueError('nonnegative energy, strict rho and positive supply required')
    return max(mp.mpf(0), v1-rho*v0)/D2


def analyze(path, certificate):
    data = np.loadtxt(path)
    if data.shape != (1801,482) or not np.array_equal(data[:,0], np.arange(1801)) or not np.isfinite(data).all():
        raise ValueError('complete finite source-bound native trace required')
    errors = data[:,2:23]
    beta = data[:,23:26]
    P = data[:,26:467].reshape(-1,21,21)
    physical = data[:,467:482]
    if not np.array_equal(P, P.transpose(0,2,1)):
        raise ValueError('asymmetric covariance')
    np.linalg.cholesky(P)
    D2 = F(certificate['physical_supply_squared_bound'])
    if np.max(np.sum(physical*physical,axis=1)) > float(D2):
        raise ValueError('physical trace exceeds analytic supply bound')
    energies = np.einsum('bi,bi->b', errors, np.linalg.solve(P,errors[...,None])[...,0])+np.sum(beta*beta,axis=1)
    endpoints = []
    with mp.workdps(60):
        for i in (0,600,1200,1800):
            audit = exact_spd(P[i].tolist())
            e = mp.matrix(errors[i].tolist())
            b = mp.matrix(beta[i].tolist())
            dual = mp.lu_solve(mp.matrix(P[i].tolist()),e)
            value = (e.T*dual)[0]+(b.T*b)[0]
            # Signed endpoint allocation: cross terms are retained, not dropped.
            groups = ('attitude','gyro_bias','velocity','position','S','acceleration','accel_bias')
            contributions = {name: mp.nstr(sum(e[k]*dual[k] for k in range(3*j,3*j+3)),60)
                             for j,name in enumerate(groups)}
            contributions['physical_beta'] = mp.nstr((b.T*b)[0],60)
            endpoints.append({'sample': i, 'V_60digits': mp.nstr(value,60),
                              'error_joint24': data[i,2:26].tolist(),
                              'active_bias': bool(data[i,1]),
                              'signed_storage_contributions': contributions,
                              'exact_covariance_SPD': audit['exact_SPD']})
        bound = mp.mpf(D2.numerator)/D2.denominator
        words = []
        for j in range(3):
            v0,v1 = [mp.mpf(endpoints[k]['V_60digits']) for k in (j,j+1)]
            words.append({'start_sample': j*600,
                          'raw_ratio_60digits': mp.nstr(v1/v0,60),
                          'prefix_max_over_start': float(max(energies[600*j:600*(j+1)+1])/float(v0)),
                          'necessary_C_by_rho': {rho: mp.nstr(necessary_supply(v0,v1,mp.mpf(rho),bound),60) for rho in RHO_TARGETS}})
    return {'endpoints': endpoints, 'words': words,
            'trace_sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--binary',type=Path)
    args = parser.parse_args()
    work = args.work_directory.resolve(); work.mkdir(parents=True,exist_ok=True)
    binary = args.binary.resolve() if args.binary else work/'driven-probe'
    if not args.binary:
        subprocess.run(['g++','-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                        '-DEIGEN_DONT_VECTORIZE','-I'+str(NATIVE.eigen_include()),
                        '-I'+str(NATIVE.ROOT/'src'),str(SOURCE),'-o',str(binary)],check=True)
    binary.chmod(binary.stat().st_mode | 0o111)
    rows = []
    for mode in ('H','A','HA'):
        for family in range(3):
            for omega in ('0.5','1'):
                cert = source_certificate(omega,family)
                path = work/f'{mode}-BIAS{family}-w{omega}.txt'
                meta = json.loads(subprocess.check_output([str(binary),mode,str(family),omega,str(path)],text=True))
                row = {'mode':mode,'family':family,'omega':omega,'native':meta,
                       'source_certificate':cert,**analyze(path,cert)}
                rows.append(row)
                print(mode,family,omega,max(float(w['raw_ratio_60digits']) for w in row['words']),flush=True)
    report = {'qualification':'OU3_ALT_DRIVEN_FINITE_STORAGE_DIAGNOSTIC_V1',
              'probe_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
              'prebuilt_binary_supplied':bool(args.binary),
              'shipping_tree_sha':subprocess.check_output(['git','rev-parse','HEAD:src'],cwd=NATIVE.ROOT,text=True).strip(),
              'rows':rows,'rho_targets':RHO_TARGETS,
              'necessary_C_is_not_a_sufficient_bound':True,
              'supply_is_predeclared_physical_coordinate_envelope_not_error_fit':True,
              'canonical_word_supply_attachment_proved':False,
              'full_source_uniform_nonlinear_rho_proved':False,
              'P4_PASS':False,'P5_MAY_START':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')

if __name__ == '__main__':
    main()
