"""Literal moving-word signed balance; no source-uniform certificate.

The existing actual-S functional is rewritten by the forced data adjoint.
Only a temporary shipping header is observed. All 18 Euclidean means, actual
rank-three gains, innovations, rotations, references, OU predictions and
post-injection bias projections are retained. The attitude recurrence is NOT
replaced by this frozen-coefficient identity.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct
import subprocess

from .ag_readout_source_diagnostic import HEADER, REPO
from .signed_temporal import atoms, jet

RECORD = struct.Struct('<II117f')
ROOT_STEP = 45000
END_STEP = 57800  # 64 seconds after the carried 225-second root


def instrument(source):
    def once(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError('signed balance observer anchor changed: '+old[:60])
        source = source.replace(old, new)

    once('    Vector12& x_lin_next = x_lin_next_scratch_;',
         '    const auto balance_before = xext.eval();\n'
         '    Vector12& x_lin_next = x_lin_next_scratch_;')
    once('    apply_pending_aw_covariance_inflation_();', '''    if (recording) {
        const T balance_phi = acc_bias_updates_enabled_
            ? std::exp(-Ts / std::max(T(1e-3), tau_bacc_)) : T(1);
        balance_prediction(F_LL, balance_phi, balance_before, xext);
    }
    apply_pending_aw_covariance_inflation_();''')
    for start, stop, kind, measured in (
        ('::measurement_update_acc_only(', '::measurement_update_mag_only(', 1, 'acc_meas'),
        ('::measurement_update_mag_only(', '::accelerometer_measurement_func(', 2, 'mag_meas'),
        ('::applyIntegralZeroPseudoMeas()', '::measurement_update_position_pseudo(', 3, 'Vector3::Zero()'),
    ):
        a, b = source.index(start), source.index(stop)
        chunk = source[a:b]
        anchor = '    xext.noalias() += K * r;'
        if chunk.count(anchor) != 1:
            raise ValueError('signed balance correction anchor changed')
        chunk = chunk.replace(anchor, '    const auto balance_before = xext.eval();\n'+anchor+
                              f'\n    balance_correction({kind}, K, S_mat, r, '
                              f'R_wb(), v2ref, {measured}, balance_before, xext);')
        source = source[:a]+chunk+source[b:]
    once('::applyQuaternionCorrectionFromErrorState()\n{',
         '::applyQuaternionCorrectionFromErrorState()\n{\n'
         '    const auto balance_before_projection = xext.eval();')
    once('    project_acc_bias_();', '''    project_acc_bias_();
    balance_projection(dtheta, balance_before_projection, xext);''')
    return source


def driver_source():
    source = (REPO/'tools/stability/ag_readout_source.cpp').read_text()
    for old, new in [('k<=45064', f'k<={END_STEP}'), ('applied!=8', 'applied!=1600')]:
        if source.count(old) != 1:
            raise ValueError('signed balance driver anchor changed')
        source = source.replace(old, new)
    tap = r'''
#include <cstdint>
#include <fstream>
static std::ofstream balance_trace;
static std::uint32_t balance_step;
template<class X> static void balance_copy(float* a, const X& x) {
    for(int i=0; i<18; ++i) a[i]=x(3+i);
}
static void balance_write(std::uint32_t kind, const float* a) {
    if (!recording) return;
    balance_trace.write(reinterpret_cast<const char*>(&kind), sizeof(kind));
    balance_trace.write(reinterpret_cast<const char*>(&balance_step), sizeof(balance_step));
    balance_trace.write(reinterpret_cast<const char*>(a), 117*sizeof(float));
}
template<class A, class X> static void balance_prediction(const A& f, float phi,
                                                         const X& before, const X& after) {
    if (!recording) return;
    float a[117]={}; balance_copy(a,before); balance_copy(a+18,after);
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) a[36+4*i+j]=f(3*i,3*j);
    a[52]=phi; balance_write(0,a);
}
template<class K, class S, class V, class R, class B, class M, class X>
static void balance_correction(int kind, const K& k, const S& s, const V& r,
                               const R& rot, const B& field, const M& measured,
                               const X& before, const X& after) {
    if (!recording) return;
    float a[117]={}; balance_copy(a,before); balance_copy(a+18,after);
    for(int i=0;i<18;++i) for(int j=0;j<3;++j) a[36+3*i+j]=k(3+i,j);
    for(int i=0;i<3;++i) {
        a[90+i]=r(i); a[102+i]=field(i); a[105+i]=measured(i);
        for(int j=0;j<3;++j) { a[93+3*i+j]=rot(i,j); a[108+3*i+j]=s(i,j); }
    }
    balance_write(kind,a);
}
template<class V, class X>
static void balance_projection(const V& d, const X& before, const X& after) {
    if (!recording) return;
    float a[117]={}; balance_copy(a,before); balance_copy(a+18,after);
    for(int i=0;i<3;++i) a[90+i]=d(i);
    balance_write(4,a);
}
'''
    source = source.replace('#define private public', tap+'\n#define private public')
    source = source.replace('    const bool wave =',
                            '    balance_trace.open(argv[2],std::ios::binary);\n    const bool wave =')
    source = source.replace('        const double t =', '        balance_step=k;\n        const double t =')
    source = source.replace('    const auto& m=filter.raw().mekf();',
                            '    balance_trace.flush();\n    if (!balance_trace) return 4;\n'
                            '    const auto& m=filter.raw().mekf();')
    source = source.replace('<< ",\\\"root_covariance\\\":" << root',
                            '<< ",\\\"committed_field\\\":" << matrix_json(m.v2ref) '
                            '<< ",\\\"root_covariance\\\":" << root')
    return source


def export(directory, eigen):
    directory.mkdir(parents=True, exist_ok=True)
    include = directory/'include'/HEADER.parent.relative_to('src')
    include.mkdir(parents=True, exist_ok=True)
    header = instrument((REPO/HEADER).read_text())
    (include/HEADER.name).write_text(header)
    wrapper = 'SeaStateFusionFilter_OU_III.h'
    (include/wrapper).write_text((REPO/HEADER.parent/wrapper).read_text())
    driver = directory/'observer.cpp'
    driver.write_text(driver_source())
    records = []
    trace = directory/'balance.bin'
    for name, inc in [('observer', ['-I'+str(directory/'include')]), ('control', [])]:
        binary = directory/name
        subprocess.run(['g++', '-O2', '-std=c++20', *inc, '-I'+str(REPO/'src'),
                        '-isystem', str(eigen), str(driver), '-o', str(binary)], check=True)
        records.append(json.loads(subprocess.check_output(
            [str(binary), 'wave', str(trace if name == 'observer' else directory/'control.bin')], text=True)))
    observed, control = records
    if observed != control:
        raise ArithmeticError('signed balance observer changed terminal output')
    with trace.open('rb') as stream:
        trace_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    meta = {'native': observed, 'literal_terminal_parity': True,
            'shipping_header_sha256': hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest(),
            'driver_sha256': hashlib.sha256(driver_source().encode()).hexdigest(),
            'instrumented_header_sha256': hashlib.sha256(header.encode()).hexdigest(),
            'trace_sha256': trace_hash}
    (directory/'export.json').write_text(json.dumps(meta, indent=2, sort_keys=True)+'\n')
    return meta


def records(directory):
    result = []
    state = None
    with (directory/'balance.bin').open('rb') as stream:
        while data := stream.read(RECORD.size):
            if len(data) != RECORD.size:
                raise ArithmeticError('truncated signed balance trace')
            kind, step, *a = RECORD.unpack(data)
            if state is not None and state != a[:18]:
                raise ArithmeticError('unrecorded Euclidean mean change')
            state = a[18:36]
            result.append((kind, step, a))
    terminal = json.loads((directory/'export.json').read_text())['native']['terminal_state'][3:]
    if state != [r[0] for r in terminal]:
        raise ArithmeticError('signed balance terminal mean mismatch')
    return result


def analyze(directory):
    import mpmath as mp
    meta = json.loads((directory/'export.json').read_text())
    trace = records(directory)
    sites = [(step, i) for i, (kind, step, _) in enumerate(trace) if kind == 3]
    selected = [sites[round(j*(len(sites)-1)/3)] for j in range(4)]
    knots = [F(step-ROOT_STEP, 200) for step, _ in selected]
    atom_by_index = {i: v for (_, i), v in zip(selected, atoms(knots))}
    site_by_step = dict(selected)
    # Include operations at the first/last time with literal left/right jets.
    start = next(i for i, (_, step, _) in enumerate(trace) if step == selected[0][0])
    stop = next((i for i, (_, step, _) in enumerate(trace) if step > selected[-1][0]), len(trace))
    with mp.workdps(80):
        def number(x):
            return mp.mpf(x.numerator)/x.denominator if isinstance(x, F) else mp.mpf(x)

        def matrix(a, rows, cols):
            return mp.matrix([[number(a[cols*i+j]) for j in range(cols)] for i in range(rows)])

        zero = lambda: mp.zeros(3, 1)
        z, direct = mp.zeros(3, 18), zero()
        terms = {key: zero() for key in ('acceleration', 'joint_rotation', 'measurement',
                                         'gravity_representation', 'innovation_arithmetic',
                                         'mean_arithmetic', 'projection')}
        weights = {'acc': [], 'mag': []}
        g = number(struct.unpack('<f', struct.pack('<f', 9.80665))[0])
        physical_g = mp.mpf('9.80665')
        ez, field = mp.matrix([0, 0, 1]), mp.matrix([60, 0, 30])
        counts = [0]*5
        for index in reversed(range(start, stop)):
            kind, step, a = trace[index]
            counts[kind] += 1
            before, after = matrix(a[:18], 18, 1), matrix(a[18:36], 18, 1)
            if kind == 0:
                f = mp.eye(18)
                for i in range(4):
                    for j in range(4):
                        for axis in range(3):
                            f[3+3*i+axis, 3+3*j+axis] = number(a[36+4*i+j])
                for axis in range(3):
                    f[15+axis, 15+axis] = number(a[52])
                terms['mean_arithmetic'] += z*(after-f*before)
                z = z*f
                continue
            if kind == 4:
                terms['projection'] += z*(after-before)
                continue
            t = F(step-ROOT_STEP, 200)
            side = 'left' if step in site_by_step and index < site_by_step[step] else 'right'
            k = matrix(a[36:90], 18, 3)
            r, rot = matrix(a[90:93], 3, 1), matrix(a[93:102], 3, 3)
            nominal_field, measured = matrix(a[102:105], 3, 1), matrix(a[105:108], 3, 1)
            w = (-number(jet(knots, t))*k[3:6, :]
                 +number(jet(knots, t, 1))*k[6:9, :]
                 -number(jet(knots, t, 2, side))*k[9:12, :])
            if index in atom_by_index:
                w += number(atom_by_index[index])*mp.eye(3)
            direct += w*r
            ell = w+z*k
            h = mp.zeros(3, 18)
            time = number(F(step, 200))
            roll = mp.mpf('.02')*mp.sin(time/2)
            q = mp.matrix([[1, 0, 0], [0, mp.cos(roll), mp.sin(roll)],
                           [0, -mp.sin(roll), mp.cos(roll)]])
            accel = mp.matrix([0, 0, -mp.mpf('.144')*mp.sin(mp.mpf('.6')*time)])
            if kind == 1:
                h[:, 12:15], h[:, 15:18] = rot, mp.eye(3)
                y = measured+g*rot*ez
                terms['acceleration'] += ell*q*accel
                terms['joint_rotation'] += g*ell*(rot-q)*ez
                terms['measurement'] += ell*(measured-q*(accel-physical_g*ez))
                terms['gravity_representation'] += (g-physical_g)*ell*q*ez
                weights['acc'].append((time, ell.copy(), ell*q))
            elif kind == 2:
                y = measured-rot*nominal_field
                terms['joint_rotation'] += ell*(q*field-rot*nominal_field)
                terms['measurement'] += ell*(measured-q*field)
                weights['mag'].append((time, ell.copy(), ell*q))
            elif kind == 3:
                h[:, 9:12] = mp.eye(3)
                y = zero()
            else:
                raise ValueError('unknown signed balance operation')
            terms['innovation_arithmetic'] += ell*(r-y+h*before)
            terms['mean_arithmetic'] += z*(after-before-k*r)
            z -= ell*h
        root = z*matrix(trace[start][2][:18], 18, 1)
        reconstructed = root+sum(terms.values(), zero())
        error = mp.norm(direct-reconstructed)
        if error > mp.mpf('1e-65'):
            raise ArithmeticError('forced data adjoint identity failed: '+str(error))
        # Supply bounds use signed total/tail sums and physical velocity before
        # norms. They are CONDITIONAL on this exported multiplier sequence.
        norm = lambda a: mp.norm(a)  # Frobenius upper bound on operator norm
        acc = list(reversed(weights['acc']))
        h = mp.mpf('.005')
        c = [wq/h for _, _, wq in acc]
        variation = sum((norm(b-a) for a, b in zip(c, c[1:])), mp.mpf(0))
        velocity_supply = mp.mpf('5.5')*(norm(c[0])+norm(c[-1])+variation)
        quadrature_supply = mp.mpf(50)*sum((norm(x)*h*h for x in c), mp.mpf(0))
        # One motivated refinement of the sampling supply: first combine the
        # signed weights over actual S-to-S intervals, then use the physical
        # velocity integral. The jerk error compares each sample with that
        # SAME interval's continuous mean, not independent acceleration noise.
        boundaries = [number(F(step, 200)) for step, i in sites if start <= i < stop]
        boundaries.append(acc[-1][0]+h)
        grouped, grouped_error, cursor = [], mp.mpf(0), 0
        for left, right in zip(boundaries, boundaries[1:]):
            total = mp.zeros(3)
            while cursor < len(acc) and acc[cursor][0] < right:
                time, _, wq = acc[cursor]
                if time < left:
                    raise ArithmeticError('sample outside its physical S interval')
                total += wq
                grouped_error += 100*norm(wq)*((time-left)**2+(right-time)**2)/(2*(right-left))
                cursor += 1
            grouped.append(total/(right-left))
        if cursor != len(acc):
            raise ArithmeticError('uncovered physical acceleration samples')
        grouped_velocity = mp.mpf('5.5')*(norm(grouped[0])+norm(grouped[-1])
            +sum((norm(b-a) for a, b in zip(grouped, grouped[1:])), mp.mpf(0)))
        tail = mp.zeros(3)
        bias_rate_action = mp.mpf(0)
        for _, w, _ in reversed(acc):
            bias_rate_action += norm(tail)*h
            tail += w
        bias_supply = mp.mpf('0.22516660498395405')*norm(tail)+mp.mpf('.001')*bias_rate_action
        acc_residual_supply = mp.mpf('.3')*sum((norm(w) for _, w, _ in acc), mp.mpf(0))
        mag_residual_supply = mp.mpf(7)*sum((norm(w) for _, w, _ in weights['mag']), mp.mpf(0))
        # This is the explicitly tested failed relaxation, not an independent
        # nominal coefficient certificate: discard only the joint rotation
        # cancellation AFTER obtaining the literal same-history multipliers.
        rotation_supply = 2*g*sum((norm(w) for _, w, _ in acc), mp.mpf(0))
        # Reference is checked constant on every actual correction, never
        # identified with the physical field or assigned its horizontal floor.
        b = mp.matrix([number(x[0]) for x in meta['native']['committed_field']])
        if any(a[102:105] != [x[0] for x in meta['native']['committed_field']]
               for kind, _, a in trace[start:stop] if kind in (1, 2, 3)):
            raise ArithmeticError('reference changed on signed balance word')
        rotation_supply += (norm(field)+norm(b))*sum((norm(w) for _, w, _ in weights['mag']), mp.mpf(0))
        threshold = g*mp.sqrt(b[0]**2+b[1]**2)/norm(b)
        fmt = lambda x: mp.nstr(x, 32)
        return {**{key: value for key, value in meta.items() if key != 'native'},
                'qualification': 'OU3_CARRIED_FORCED_DATA_ADJOINT_DIAGNOSTIC_V1',
                'profile': 'roll .02 sin(t/2), vertical displacement .4 sin(.6t); construction to 289 s',
                'balance_time_interval_s': [str(F(trace[start][1]-1, 200)), str(F(trace[stop-1][1], 200))],
                'knots_s': list(map(str, knots)), 'operation_counts_prediction_acc_mag_S_projection': counts,
                'decimal_digits': 80, 'identity_residual_norm': fmt(error),
                'direct_signed_functional': [fmt(x) for x in direct],
                'root_term': [fmt(x) for x in root],
                'signed_terms': {key: [fmt(x) for x in value] for key, value in terms.items()},
                'root_multiplier_block_frobenius_norms': [fmt(norm(z[:, i:i+3])) for i in range(0, 18, 3)],
                'root_multiplier_blocks': ['BG', 'V', 'P', 'S', 'AW', 'BA'],
                'physical_gravity': fmt(physical_g), 'literal_float_gravity': fmt(g),
                'conditional_supplies': {
                    'physical_velocity': fmt(velocity_supply), 'jerk_quadrature': fmt(quadrature_supply),
                    'S_interval_velocity': fmt(grouped_velocity), 'S_interval_jerk': fmt(grouped_error),
                    'physical_accel_bias_Abel': fmt(bias_supply), 'accel_fast_residual': fmt(acc_residual_supply),
                    'gravity_representation': fmt(abs(g-physical_g)*sum((norm(w) for _, w, _ in acc), mp.mpf(0))),
                    'magnetic_residual': fmt(mag_residual_supply), 'rotation_triangle': fmt(rotation_supply)},
                'recorded_nominal_projected_gravity_threshold': fmt(threshold),
                'S_interval_acceleration_margin_before_other_supplies': fmt(threshold-grouped_velocity-grouped_error),
                'rotation_triangle_margin_before_other_supplies': fmt(threshold-rotation_supply),
                'source_uniform_verified': False, 'all_time_magnetic_service_certified': False,
                'symbolic_attitude_excitation_parameters_qualified': False,
                'Delta_col_certified': False, 'Delta_gyr_certified': False, 'theorem_closed': False}


def verify_diagnostic(report):
    """Scope/fingerprint/rounded-balance checks, NOT a trajectory enclosure.

    Native reproduction is a separate workflow step. No finite decimal
    residual, recorded gain or reference is a source-uniform certificate.
    """
    if report['qualification'] != 'OU3_CARRIED_FORCED_DATA_ADJOINT_DIAGNOSTIC_V1':
        raise ValueError('invalid forced data diagnostic qualification')
    for flag in ('source_uniform_verified', 'all_time_magnetic_service_certified',
                 'symbolic_attitude_excitation_parameters_qualified',
                 'Delta_col_certified', 'Delta_gyr_certified', 'theorem_closed'):
        if report[flag] is not False:
            raise ValueError('forced data diagnostic cannot promote '+flag)
    if report['literal_terminal_parity'] is not True or report['decimal_digits'] < 80:
        raise ValueError('literal parity and high precision required')
    source = (REPO/HEADER).read_text()
    for key, value in (('shipping_header_sha256', source), ('driver_sha256', driver_source()),
                       ('instrumented_header_sha256', instrument(source))):
        if report[key] != hashlib.sha256(value.encode()).hexdigest():
            raise ValueError('forced data diagnostic fingerprint changed: '+key)
    if not 0 <= F(report['identity_residual_norm']) < F('1e-65'):
        raise ValueError('forced data diagnostic identity residual')
    for axis in range(3):
        total = F(report['root_term'][axis])+sum(F(v[axis]) for v in report['signed_terms'].values())
        if abs(total-F(report['direct_signed_functional'][axis])) > F('1e-30'):
            raise ValueError('rounded signed balance inconsistent')
    threshold = F(report['recorded_nominal_projected_gravity_threshold'])
    bounds = report['conditional_supplies']
    for key, charge in (
        ('rotation_triangle_margin_before_other_supplies', F(bounds['rotation_triangle'])),
        ('S_interval_acceleration_margin_before_other_supplies',
         F(bounds['S_interval_velocity'])+F(bounds['S_interval_jerk'])),
    ):
        if abs(F(report[key])-(threshold-charge)) > F('1e-29'):
            raise ValueError('rounded failed supply margin inconsistent')
    return True


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--export', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.export:
        export(args.directory, args.eigen)
    report = analyze(args.directory)
    verify_diagnostic(report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
