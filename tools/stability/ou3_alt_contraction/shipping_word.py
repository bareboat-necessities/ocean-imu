#!/usr/bin/env python3
"""Actual paired shipping-word attachment experiment; never a proof gate.

Temporary read-only hooks expose actual P/N/S/K/residuals; no deployed source
is edited. A continuous analytic source drives the full private startup and
frontend. Selected finite perturbations are falsification probes, not a cover
of COMPLETE-BRMM, a hard entry set, or a metric training set.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/"tools/stability"))
from tools.stability.ou3_alt_contraction.shipping_graph import binding_report
HEADER = Path('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h')
FUSION = Path('src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h')
# Binary payload, always row-major; Eigen's in-memory order is not serialized.
FIELDS = [('q', 4), ('x', 21), ('P', 441), ('N', 63), ('S', 9),
          ('K', 63), ('r', 3), ('R', 9), ('source', 27), ('frontend', 14),
          ('measured',3), ('mag_reference',3), ('parameters',6), ('hard_iron',3)]
DTYPE = np.dtype([('kind','<u4'),('sample','<u4'),('active','<u4')] +
                 [(name,'<f4',(size,)) for name,size in FIELDS])
KIND = {1:'prediction_entry',2:'prediction_exit',3:'floor_exit',
        10:'accelerometer',11:'magnetometer',12:'S_zero',
        20:'accelerometer_exit',21:'magnetometer_exit',22:'S_zero_exit',
        30:'bias_mode_entry',31:'bias_mode_exit',40:'sample_exit',
        50:'accelerometer_mean',51:'magnetometer_mean',52:'S_zero_mean',
        60:'accelerometer_Joseph',61:'magnetometer_Joseph',62:'S_zero_Joseph',
        70:'quaternion_covariance_reset',80:'bias_projection',
        110:'accelerometer_attempt',111:'magnetometer_attempt',112:'S_zero_attempt',
        120:'accelerometer_return',121:'magnetometer_return',122:'S_zero_return'}


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'shipping hook anchor changed: expected once: {old!r}')
    return text.replace(old, new, 1)


def instrument_header(text):
    """Insert pure observation calls at verified lexical anchors, fail closed."""
    if 'static constexpr T tempC_ref = T(35.0);' not in text:
        raise ValueError('probe temperature no longer equals shipping reference')
    # Method-scoped replacement keeps unrelated pseudo-measurements untouched.
    for method, kind, noise, measured in [
            ('measurement_update_acc_only',10,'Racc','f_meas'),
            ('measurement_update_mag_only',11,'Rmag','mag_meas'),
            ('applyIntegralZeroPseudoMeas',12,'R_S','Vector3::Zero()')]:
        marker = f'::'+method+'('
        start = text.index(marker)
        stop = text.find('\ntemplate<', start+len(marker))
        if stop == -1: stop = len(text)
        body = text[start:stop]
        brace = body.index('{')
        body = (body[:brace+1] + f'\n    alt::point({kind+100}, *this);\n'
                + f'    alt::ExitGuard alt_guard({kind+110}, *this);\n' + body[brace+1:])
        anchor = 'xext.noalias() += K * r;'
        body = once(body, anchor,
            f'alt::solve({kind}, *this, r, PCt, S_mat, K, {noise}, {measured});\n    '
            +anchor+f'\n    alt::point({kind+40}, *this);')
        anchor = 'joseph_update3_(K, S_mat, PCt);'
        body = once(body, anchor, anchor+f'\n    alt::point({kind+50}, *this);')
        anchor = 'applyQuaternionCorrectionFromErrorState();'
        body = once(body, anchor, anchor+f'\n    alt::point({kind+10}, *this);')
        text = text[:start]+body+text[stop:]
    anchor = '    apply_pending_aw_covariance_inflation_();\n    symmetrize_Pext_();   // Symmetry hygiene'
    text = once(text, anchor,
        '    alt::point(2, *this);\n'+anchor+'\n    alt::point(3, *this);')
    start = text.index('::time_update(')
    brace = text.index('{', start)
    text = text[:brace+1]+'\n    alt::point(1, *this);'+text[brace+1:]
    anchor = '    void set_acc_bias_updates_enabled(bool en) {\n'
    text = once(text, anchor, anchor+'        alt::point(30, *this);\n')
    anchor = '        acc_bias_updates_enabled_ = en;\n'
    text = once(text, anchor, anchor+'        alt::point(31, *this);\n')
    anchor = '    project_acc_bias_();\n}'
    text = once(text, anchor, '    alt::point(70, *this);\n' + anchor[:-2]
                + '\n    alt::point(80, *this);\n}')
    return text


def instrument_fusion(text):
    anchor='            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);'
    return once(text,anchor,'            alt::set_hard_iron(mag_hard_iron_body_uT_);\n'+anchor)


def gain_pair_identity(a, b):
    """Actual row-solve graph with explicit defects; no inverse and no r0=0.

    Eigen LDLT's default Lower view solves the symmetric matrix reconstructed
    from S's lower triangle. The raw (possibly asymmetric) S remains separate
    for shipping Joseph. E_i = Ssolve_i K_i^T - N_i^T is a measured numerical
    defect, not a rigorous all-input error bound.
    """
    def unpack(r):
        raw = r['S'].reshape(3,3).astype(float)
        s = np.tril(raw)+np.tril(raw,-1).T
        return s, r['N'].reshape(21,3).astype(float), r['K'].reshape(21,3).astype(float), r['r'].astype(float)
    s0,n0,k0,r0 = unpack(a); s1,n1,k1,r1 = unpack(b)
    ds,dn,dk,dr = s1-s0,n1-n0,k1-k0,r1-r0
    e0,e1 = s0@k0.T-n0.T, s1@k1.T-n1.T
    equality = s1@dk.T+ds@k0.T-dn.T-(e1-e0)
    full = k1@r1-k0@r0
    frozen = k0@dr
    gain_term = dk@r1  # exact alternate decomposition full=K0*dr+dK*r1
    return dict(equality_max_abs=float(np.max(np.abs(equality))),
                correction_identity_max_abs=float(np.max(np.abs(full-frozen-gain_term))),
                nominal_residual_norm=float(np.linalg.norm(r0)),
                dN_norm=float(np.linalg.norm(dn)), dS_norm=float(np.linalg.norm(ds)),
                endogenous_correction_norm=float(np.linalg.norm(gain_term)),
                frozen_correction_norm=float(np.linalg.norm(frozen)),
                total_correction_norm=float(np.linalg.norm(full)),
                row_solve_defect_max_abs=float(max(np.max(np.abs(e0)),np.max(np.abs(e1)))),
                raw_S_asymmetry_max_abs=float(max(np.max(np.abs(a['S'].reshape(3,3)-a['S'].reshape(3,3).T)),np.max(np.abs(b['S'].reshape(3,3)-b['S'].reshape(3,3).T)))))


def high_precision_pair(a,b,dps):
    """Recompute the *finite row-solve identity*, not an omitted nonlinear word."""
    import mpmath as mp
    with mp.workdps(dps):
        def mat(values,rows,cols):
            vals=[mp.mpf(float(v)) for v in values]
            return mp.matrix([[vals[i*cols+j] for j in range(cols)] for i in range(rows)])
        def unpack(r):
            s=mat(r['S'],3,3)
            for i in range(3):
                for j in range(i+1,3):s[i,j]=s[j,i]
            return s,mat(r['N'],21,3),mat(r['K'],21,3),mat(r['r'],3,1)
        s0,n0,k0,r0=unpack(a);s1,n1,k1,r1=unpack(b)
        e0=s0*k0.T-n0.T;e1=s1*k1.T-n1.T
        eq=s1*(k1-k0).T+(s1-s0)*k0.T-(n1-n0).T-(e1-e0)
        gap=k1*r1-k0*r0-k0*(r1-r0)-(k1-k0)*r1
        return {'digits':dps,'row_graph_residual':str(max(abs(v) for v in eq)),
                'correction_graph_residual':str(max(abs(v) for v in gap)),
                'scope':'finite recorded row-solve graph, NOT high-precision nonlinear re-execution'}


def read_trace(path):
    if path.stat().st_size % DTYPE.itemsize:
        raise ValueError(f'truncated trace {path}')
    rows=np.fromfile(path,dtype=DTYPE)
    if not len(rows): raise ValueError(f'empty trace {path}')
    if not np.isin(rows['kind'],list(KIND)).all():raise ValueError('unknown trace event')
    for field,_ in FIELDS:
        if not np.isfinite(rows[field]).all():
            raise ValueError(f'nonfinite {field} in {path.name}')
    return rows


def reference_error(row):
    """Diagnostic true-minus-estimate joint24 at a sample boundary.

    The reference orientation uses the recorded physical Euler angles and the
    principal SO(3) log. Their evaluation/roundoff and chart domain are not
    certified. Temporary uninjected error-state slots are not endpoint states.
    """
    from scipy.spatial.transform import Rotation
    if int(row['kind'])!=40:raise ValueError('reference error requires a sample boundary')
    source=row['source'].astype(float);x=row['x'].astype(float)
    q=row['q'].astype(float)
    estimate=Rotation.from_quat(np.r_[q[1:],q[0]])
    truth=Rotation.from_euler('xyz',source[24:27]).inv()
    theta=(truth*estimate.inv()).as_rotvec()
    return np.r_[theta,-x[3:6],source[3:6]-x[6:9],source[:3]-x[9:12],
                 source[9:12]-x[12:15],source[6:9]-x[15:18],
                 source[18:21]-x[18:21],source[18:21]]


def state_difference(a,b):
    return reference_error(b)-reference_error(a)


def analyze_pair(base,other):
    # Mismatched guard paths are an explicit attachment failure, never zipped
    # away or repaired by dropping the extra event.
    same=np.array_equal(base[['kind','sample','active']],other[['kind','sample','active']])
    result={'event_paths_equal':same,'records_baseline':len(base),'records_perturbed':len(other)}
    if not same:
        result['classification']='HYBRID_PAIR_PATH_MISMATCH_REQUIRES_BRANCH_GRAPH'
        return result
    solve=np.flatnonzero(np.isin(base['kind'],[10,11,12]))
    if not len(solve): raise ValueError('word has no actual measurement solves')
    stats=[gain_pair_identity(base[i],other[i]) for i in solve]
    worst=int(np.argmax([s['endogenous_correction_norm'] for s in stats])); idx=int(solve[worst])
    start,end=state_difference(base[0],other[0]),state_difference(base[-1],other[-1])
    result.update(qualification='ACTUAL_FINITE_PAIRED_WORD_ATTACHMENT_DIAGNOSTIC',
        solve_count=len(solve), counts={KIND[k]:int(np.count_nonzero(base['kind']==k)) for k in KIND},
        initial_joint24_increment=start.tolist(), final_joint24_increment=end.tolist(),
        initial_covariance_increment_norm=float(np.linalg.norm(base[0]['P'].astype(float)-other[0]['P'].astype(float))),
        max_frontend_increment=float(np.max(np.abs(base['frontend'].astype(float)-other['frontend'].astype(float)))),
        max_mag_reference_increment=float(np.max(np.abs(base['mag_reference'].astype(float)-other['mag_reference'].astype(float)))),
        max_hard_iron_increment=float(np.max(np.abs(base['hard_iron'].astype(float)-other['hard_iron'].astype(float)))),
        largest_endogenous_term={**stats[worst],'sample':int(base[idx]['sample']),'event':KIND[int(base[idx]['kind'])]},
        maximum_row_identity_residual=max(s['equality_max_abs'] for s in stats),
        maximum_correction_identity_residual=max(s['correction_identity_max_abs'] for s in stats),
        maximum_measured_row_solve_defect=max(s['row_solve_defect_max_abs'] for s in stats),
        checks_80_120=[high_precision_pair(base[idx],other[idx],d) for d in (80,120)],
        source_uniform_cover=False, hard_entry_qualified=False, common_metric_searched=False,
        rho=None, contraction_margin=None, maximum_generalized_direction=None,
        prefix_chart_retention_certified=False, finite_precision_enclosure=False)
    # A raw Euclidean quantity with mixed units is deliberately not called rho.
    return result


def source_contract():
    """Sufficient analytic kinematic bounds, not the complete theorem's PE cover."""
    import math
    w=math.pi/4
    p=math.sqrt(0.3**2+0.2**2+1.0**2)
    return {'type':'continuous harmonic kinematic subfamily; analytic bounds, no sampled maxima',
        'position_norm_upper_m':p,'velocity_norm_upper_mps':w*p,
        'acceleration_norm_upper_mps2':w*w*p,'centered_primitive_any_u_t_upper_ms':2*p/w,
        'body_rate_norm_upper_radps':w*(.12+1.25*.09+.75*.08),
        'translation_frequency_hz':.125,
        'S_origin':'first actual Live sample, retained across all word boundaries',
        'sensor_equation':'acc=R_bw^T*(a-g*ez)+beta; gyro=exact Euler body rate; mag=R_bw^T*B',
        'lever_arm':0,'temperature_C':35,
        'BRMM_full_admission':False,'PE_admission':False,
        'kinematic_envelope_sufficient':p<=8.1 and w*p<=5.5 and w*w*p<=8.8 and 2*p/w<=1100,
        'probes_are_not_source_cover':True}


def bias_contracts():
    """Read-only definitions plus analytic membership of the chosen examples."""
    import importlib, math
    result={}
    for i in range(3):
        module=importlib.import_module(f'tools.stability.ou3_p4_bias{i}_family')
        definition=module.build();failures=module.validate(definition)
        if failures:raise ValueError(f'BIAS{i} definition failed: {failures}')
        if i==0:
            ok=(definition['gauss_markov_component_abs_upper_mps2']>=.005
                and definition['gauss_markov_tau_true_s'][0]<=900<=definition['gauss_markov_tau_true_s'][1]
                and definition['turn_on_offset_component_abs_upper_mps2']>=.025
                and definition['non_gauss_markov_component_abs_upper_mps2']>=.002
                and definition['non_gauss_markov_rate_abs_upper_mps3']>=.002*2*math.pi/900)
            channels='GM=.005*exp(-t/900), offset<=.025, nonGM sine .002/900s; thermal=strain=0'
        elif i==1:
            ok=(definition['root_component_abs_upper_mps2']>=.05
                and definition['tau_true_s'][0]<=1200<=definition['tau_true_s'][1]
                and definition['sinusoid_component_amplitude_abs_upper_mps2']>=.01
                and definition['sinusoid_period_s'][0]<=600<=definition['sinusoid_period_s'][1])
            channels='root<=.05, tau=1200s, sine .01/600s, deterministic mismatch=0'
        else:
            ok=(definition['true_bias_component_abs_upper_mps2']>=.045
                and definition['variation_rate_abs_upper_mps3']>=.005*2*math.pi/900)
            channels='phi_true=1; offset<=.04 plus .005/900s sine; magnitude<=.045'
        if not ok:raise ValueError(f'BIAS{i} analytic example no longer fits definition')
        result[f'BIAS{i}']={'definition_validator_pass':True,'analytic_example_parameter_membership':ok,
            'channels':channels,'declared_driver_component_upper':definition['driver_increment_component_abs_upper_mps2'],
            'hardware_admission':False,'fresh_Live_admission':False,'binary32_source_roundoff_enclosed':False}
    return result


def run(output,work,eigen,verify_uninstrumented=False):
    work.mkdir(parents=True,exist_ok=True)
    overlay=work/'include'; dst=overlay/'kalman_ou_iii'; dst.mkdir(parents=True,exist_ok=True)
    (dst/HEADER.name).write_text(instrument_header((ROOT/HEADER).read_text()))
    (dst/FUSION.name).write_text(instrument_fusion((ROOT/FUSION).read_text()))
    binary=work/'shipping-word'
    command=[os.environ.get('CXX','g++'),'-O1','-DEIGEN_UNROLLING_LIMIT=0','-std=c++20','-DEIGEN_NON_ARDUINO',
        '-I'+str(overlay),'-I'+str(ROOT/'src'),'-isystem',str(eigen),
        str(ROOT/'tests/ou3_alt_contraction/shipping_word.cpp'),'-o',str(binary)]
    fingerprints={str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
        for p in (HEADER,FUSION,Path('src/kalman_ou_common/KalmanOUCoreMath.h'),
                  Path('tests/ou3_alt_contraction/shipping_word.cpp'))}
    producer={'source_sha256':fingerprints,'compiler_command':command,
              'compiler_version':subprocess.check_output([command[0],'--version'],text=True).splitlines()[0],
              'temporary_Kalman_sha256':hashlib.sha256((dst/HEADER.name).read_bytes()).hexdigest(),
              'temporary_Fusion_sha256':hashlib.sha256((dst/FUSION.name).read_bytes()).hexdigest()}
    (work/'producer.json').write_text(json.dumps(producer,indent=2)+'\n')
    subprocess.run(command,check=True)
    subprocess.run([str(binary),str(work)],check=True)
    if verify_uninstrumented:
        plain=work/'uninstrumented';plain.mkdir(exist_ok=True)
        plain_command=[v for v in command if v!='-I'+str(overlay)]
        plain_command[-1]=str(plain/'shipping-word')
        subprocess.run(plain_command,check=True)
        subprocess.run([plain_command[-1],str(plain)],check=True)
        producer['transparency']=verify_transparency(work,plain)
        (work/'producer.json').write_text(json.dumps(producer,indent=2)+'\n')
    return analyze(output,work,command)


def verify_transparency(hooked,plain):
    manifest=json.loads((hooked/'manifest.json').read_text())
    if manifest!=json.loads((plain/'manifest.json').read_text()):
        raise ValueError('observation hooks changed actual runtime boundaries')
    samples=0
    for word in manifest['words']:
        for direction in ('base','tilt','held_bias','covariance'):
            name=word['name']+'-'+direction+'.bin'
            a=read_trace(hooked/name);b=read_trace(plain/name)
            a=a[a['kind']==40];b=b[b['kind']==40]
            if a.tobytes()!=b.tobytes():raise ValueError('hooks changed sample state: '+name)
            samples+=len(a)
    return {'sample_records_compared':samples,'all_sample_records_bit_identical':True,
            'scope':'these executions and this compiler, NOT an all-input instrumentation theorem'}


def analyze(output,work,command=None):
    manifest=json.loads((work/'manifest.json').read_text())
    report={'schema':1,'qualification':'ACTUAL_SHIPPING_FINITE_INCREMENT_ATTACHMENT_EXPERIMENT',
        'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
        'P4_promoted':False,'P5_promoted':False,'source':source_contract(),
        'source_sha256':{str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (HEADER,FUSION)},
        'execution':manifest,'pairs':{},'operand_bindings':{},'reference_words':{},'compiler_command':command,
        'analysis_sha256':{str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
            (Path('tools/stability/ou3_alt_contraction/shipping_word.py'),
             Path('tools/stability/ou3_alt_contraction/shipping_graph.py'))},
        'producer':json.loads((work/'producer.json').read_text()) if (work/'producer.json').exists() else None,
        'bias_definitions':bias_contracts(),
        'unclosed':['source-uniform joint24 word and reachable covariance/frontend pair relation',
                    'coercive compatible storage with justified bias/source/roundoff supply',
                    'full nonlinear high-precision complete-word feasibility and every-prefix charts',
                    'fresh-Live hard membership, all hybrid guards and deployment roundoff enclosure']}
    for word in manifest['words']:
        stem=word['name'];base=read_trace(work/(stem+'-base.bin'))
        family=int(stem[4])
        report['reference_words'][stem]={'initial_joint24_error':reference_error(base[0]).tolist(),
            'final_joint24_error':reference_error(base[-1]).tolist(),
            'qualification':'diagnostic physical reference chart, no hard membership or coercivity claim'}
        report['operand_bindings'][stem+'/base']=binding_report(base,family,manifest['dt'])
        for direction in ('tilt','held_bias','covariance'):
            other=read_trace(work/(stem+'-'+direction+'.bin'))
            report['pairs'][stem+'/'+direction]=analyze_pair(base,other)
            report['operand_bindings'][stem+'/'+direction]=binding_report(other,family,manifest['dt'])
    if not all(v['measurement_binding_regression_pass'] for v in report['operand_bindings'].values()):
        report['execution_failure']='measurement operand binding regression failed'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'output':str(output),'pairs':len(report['pairs']),
                      'ALT_LIVE_PASS':False,'source_uniform_cover':False}))
    if 'execution_failure' in report:raise RuntimeError(report['execution_failure'])
    return report


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--work',type=Path,required=True)
    ap.add_argument('--verify-uninstrumented',action='store_true',
        help='also compile without hooks and demand bit-identical sample-boundary traces')
    ap.add_argument('--eigen',type=Path,default=Path('/usr/include/eigen3'))
    a=ap.parse_args()
    run(a.output.resolve(),a.work.resolve(),a.eigen.resolve(),a.verify_uninstrumented)
