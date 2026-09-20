"""Carried six-mean action: a construction-linked exclusion attempt.

The controlling use is to exclude gyro aliases and nominal force/field
collinearity before constructing a uniform historical AG reader. This is
not the reader noise action B_W. No source-uniform theorem is promoted.

Only a temporary header is instrumented. Binary float operands, actual
accepted gains and innovations are exported from the complete construction.
The same recurrence runs at 80 digits, then with outward interval arithmetic.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct
import subprocess
from dataclasses import dataclass

from .construction_history_diagnostic import REPO, driver_source

HEADER = Path('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h')
RECORD = struct.Struct('<I43f')


@dataclass(frozen=True)
class RationalInterval:
    """Exact endpoint arithmetic for the small committed factor certificate."""
    lo: F
    hi: F

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('reversed interval')

    @staticmethod
    def cast(value):
        return value if isinstance(value, RationalInterval) else RationalInterval(F(value), F(value))

    def __add__(self, other):
        b = self.cast(other)
        return RationalInterval(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return RationalInterval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        b = self.cast(other)
        products = [x*y for x in (self.lo, self.hi) for y in (b.lo, b.hi)]
        return RationalInterval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.cast(other)
        if b.lo <= 0 <= b.hi:
            raise ArithmeticError('interval denominator contains zero')
        return self*RationalInterval(1/b.hi, 1/b.lo)

    def __gt__(self, other):
        return self.lo > self.cast(other).hi


def verify_summary(report, field):
    """Recheck the decisive inequalities with exact rational interval LDL."""
    def read(pair):
        return RationalInterval(*(F(x) for x in pair))
    z, one = RationalInterval.cast(0), RationalInterval.cast(1)
    d = [[read(x) for x in row] for row in report['joint_action']]
    if len(d) != 6 or any(len(row) != 6 for row in d):
        raise ValueError('six-coordinate full matrix required')
    if any(d[i][j] != d[j][i] for i in range(6) for j in range(6)):
        raise ValueError('symmetric joint action required')
    ldlt(d, z, one)
    e = read(report['innovation_energy'])
    rounding = read(report['gyro_roundoff_prefix_l1'])
    # This one-rad/s ceiling is proved from the exported full matrix action;
    # it is not imposed on the source admissible domain.
    radius = F(1)-rounding
    if not radius > 0:
        raise ArithmeticError('positive gyro barrier radius required')
    barrier = [[(radius*radius if i == j else z)-e*d[i][j] for j in range(3)] for i in range(3)]
    _, pivots = ldlt(barrier, z, one)
    b = [F(x) for x in field]
    p = max(range(3), key=lambda i: abs(b[i]))
    if b[p] == 0:
        raise ValueError('nonzero committed field required')
    c = []
    for i in range(3):
        if i != p:
            row = [F(0)]*3
            row[i], row[p] = b[p], -b[i]
            c.append(row)
    eta = [read(x) for x in report['mean_roundoff'][3:]]
    gravity = [F(0), F(0), F(struct.unpack('<f', struct.pack('<f', 9.80665))[0])]
    gram = [[sum((c[i][u]*d[3+u][3+v]*c[j][v] for u in range(3) for v in range(3)), z) for j in range(2)] for i in range(2)]
    target = [sum((c[i][j]*(gravity[j]-eta[j]) for j in range(3)), z) for i in range(2)]
    cost = inverse_action(gram, target, z, one)
    margin = cost-e
    if margin.hi >= 0:
        raise ArithmeticError('strict failure of the collinearity exclusion not certified')
    # Small outward decimal bounds keep the evidence reviewable; the exact
    # input endpoints and all rational cross terms above remain in the audit.
    def down(x):
        return F(x*10**6//1, 10**6)
    return {'qualification': 'OU3_RECORDED_CONSTRUCTION_MEAN_EXCLUSION_V1',
            'joint_action_SPD_verified': True,
            'gyro_barrier_pivot_lower': [str(down(x.lo)) for x in pivots],
            'collinearity_exclusion_margin_lower': str(down(margin.lo)),
            'collinearity_exclusion_margin_upper': str(-down(-margin.hi)),
            'collinear_energy_upper': str(-down(-cost.hi)),
            'finite_recorded_gyro_increment_below_three_verified': True,
            'finite_recorded_gyro_mean_norm_below_one_verified': True,
            'finite_recorded_angular_increment_upper': str(F('.006')*(F('0.6108652381980153')+F('.04')+1)),
            'nominal_collinearity_excluded_by_this_action': False,
            'source_uniform_verified': False, 'theorem_closed': False}


def recorded_force_separation(directory):
    """Exact finite-prefix force test; no free innovations or alternate means.

The default wrapper holds the committed reference after refinement. Check
every actual pre-accelerometer state in the same 400--600 s recorded tail.
"""
    native = json.loads((directory/'export.json').read_text())['native']
    if native['refined_step'] >= 80001:
        raise ValueError('committed reference not established before this tail')
    b = [F(x[0]) for x in native['committed_field']]
    g = F(struct.unpack('<f', struct.pack('<f', 9.80665))[0])
    b2 = sum(x*x for x in b)
    predictions = accepted = checked = 0
    margin, worst_force = None, None
    with (directory/'mean-action.bin').open('rb') as stream:
        while data := stream.read(RECORD.size):
            kind, *a = RECORD.unpack(data)
            if kind == 0:
                predictions += 1
            elif a[42] == 1:
                accepted += 1
                if native['live_step']+predictions <= 80000:
                    continue
                f = [F(x) for x in a[33:36]]
                f[2] -= g
                f2 = sum(x*x for x in f)
                dot = sum(f[i]*b[i] for i in range(3))
                # sin(angle)^2 > (2/5)^2, with exact binary operands.
                candidate = 21*f2*b2-25*dot*dot
                if margin is None or candidate < margin:
                    margin, worst_force = candidate, f
                checked += 1
    if accepted != predictions or checked != 40000 or not margin > 0:
        raise ArithmeticError('recorded sustained force separation not certified')
    return {'qualification': 'OU3_RECORDED_FORCE_SEPARATION_V1',
            'pre_accelerometer_prefixes': checked,
            'all_predictions_have_an_applied_accelerometer_correction': True,
            'sine_lower': '2/5', 'polynomial_margin_lower': str(margin),
            'worst_recorded_world_force': [str(x) for x in worst_force],
            'recorded_tail_verified': True, 'source_uniform_verified': False,
            'theorem_closed': False}


def instrument(source):
    """Tap only the closed (b_g,a_w) mean subsystem; retain all gain columns."""
    anchor = '    Vector12& x_lin_next = x_lin_next_scratch_;'
    if source.count(anchor) != 1:
        raise ValueError('prediction observer anchor changed')
    source = source.replace(anchor, '    const auto action_before = xext.eval();\n'+anchor)
    anchor = '    xext.template segment<3>(OFF_AW) = x_lin_next.template segment<3>(9);'
    if source.count(anchor) != 1:
        raise ValueError('prediction completion anchor changed')
    source = source.replace(anchor, anchor+'\n    mean_prediction(F_LL(9,9), action_before, xext);')
    for start, stop in (
        ('::measurement_update_acc_only(', '::measurement_update_mag_only('),
        ('::measurement_update_mag_only(', '::accelerometer_measurement_func('),
        ('::applyIntegralZeroPseudoMeas()', '::measurement_update_position_pseudo('),
    ):
        a, b = source.index(start), source.index(stop)
        chunk = source[a:b]
        anchor = '    xext.noalias() += K * r;'
        if chunk.count(anchor) != 1:
            raise ValueError('accepted correction observer anchor changed')
        is_acc = 1 if start == '::measurement_update_acc_only(' else 0
        chunk = chunk.replace(anchor, '    const auto action_before = xext.eval();\n'+anchor+
                              f'\n    mean_correction(K, S_mat, r, action_before, xext, {is_acc});')
        source = source[:a]+chunk+source[b:]
    return source


def observer_source():
    s = driver_source()
    tap = r'''
#include <cstdint>
#include <fstream>
static std::ofstream mean_trace;
template<class X> static void mean_copy(float* dst, const X& x) {
    for(int i=0; i<3; ++i) { dst[i]=x(3+i); dst[3+i]=x(15+i); }
}
static void mean_write(std::uint32_t kind, const float* a) {
    mean_trace.write(reinterpret_cast<const char*>(&kind), sizeof(kind));
    mean_trace.write(reinterpret_cast<const char*>(a), 43*sizeof(float));
}
template<class X> static void mean_prediction(float phi, const X& before, const X& after) {
    float a[43]={}; a[0]=phi; mean_copy(a+1,before); mean_copy(a+7,after); mean_write(0,a);
}
template<class K, class S, class R, class X>
static void mean_correction(const K& k, const S& s, const R& r, const X& before, const X& after, int is_acc) {
    float a[43]={};
    for(int i=0;i<6;++i) for(int j=0;j<3;++j) a[3*i+j]=k(i<3?3+i:12+i,j);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) a[18+3*i+j]=s(i,j);
    for(int i=0;i<3;++i) a[27+i]=r(i);
    mean_copy(a+30,before); mean_copy(a+36,after); a[42]=static_cast<float>(is_acc); mean_write(1,a);
}
'''
    s = s.replace('#define private public', tap+'\n#define private public')
    s = s.replace('    const bool wave =', '    mean_trace.open(argv[2],std::ios::binary);\n    const bool wave =')
    # Record the actually committed reference after refinement as well.
    s = s.replace('<< ",\\\"root_covariance\\\":" << root',
                  '<< ",\\\"committed_field\\\":" << matrix_json(m.v2ref) << ",\\\"root_covariance\\\":" << root')
    return s


def export(directory, eigen):
    """Run the full literal history and check untapped-control parity."""
    directory.mkdir(parents=True, exist_ok=True)
    include = directory/'include'/HEADER.parent.relative_to('src')
    include.mkdir(parents=True, exist_ok=True)
    header = instrument((REPO/HEADER).read_text())
    (include/HEADER.name).write_text(header)
    wrapper = 'SeaStateFusionFilter_OU_III.h'
    (include/wrapper).write_text((REPO/HEADER.parent/wrapper).read_text())
    source = observer_source()
    driver = directory/'observer.cpp'
    driver.write_text(source)
    binary = directory/'observer'
    subprocess.run(['g++', '-O2', '-std=c++20', '-I'+str(directory/'include'),
                    '-I'+str(REPO/'src'), '-isystem', str(eigen), str(driver), '-o', str(binary)], check=True)
    trace = directory/'mean-action.bin'
    native = json.loads(subprocess.check_output([str(binary), 'wave', str(trace)], text=True))
    control_driver, control_binary = directory/'control.cpp', directory/'control'
    control_driver.write_text(driver_source())
    subprocess.run(['g++', '-O2', '-std=c++20', '-I'+str(REPO/'src'),
                    '-isystem', str(eigen), str(control_driver), '-o', str(control_binary)], check=True)
    control = json.loads(subprocess.check_output([str(control_binary), 'wave'], text=True))
    if {k: native[k] for k in control} != control:
        raise ArithmeticError('observer differs from the fresh untapped construction control')
    with trace.open('rb') as stream:
        trace_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    meta = {'native': native, 'observer_control_parity': True,
            'observer_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'instrumented_header_sha256': hashlib.sha256(header.encode()).hexdigest(),
            'trace_sha256': trace_hash}
    (directory/'export.json').write_text(json.dumps(meta, indent=2, sort_keys=True)+'\n')
    return meta


def ldlt(a, zero=F(0), one=F(1)):
    """Small exact/interval factorization; rejects an unproved positive pivot."""
    n = len(a)
    l = [[one if i == j else zero for j in range(n)] for i in range(n)]
    d = []
    for i in range(n):
        p = a[i][i]-sum((l[i][k]*l[i][k]*d[k] for k in range(i)), zero)
        if not p > 0:
            raise ArithmeticError('positive factor pivot not certified')
        d.append(p)
        for j in range(i+1, n):
            l[j][i] = (a[j][i]-sum((l[j][k]*l[i][k]*d[k] for k in range(i)), zero))/p
    return l, d


def inverse_action(a, x, zero=F(0), one=F(1)):
    l, d = ldlt(a, zero, one)
    y = []
    for i in range(len(x)):
        y.append(x[i]-sum((l[i][j]*y[j] for j in range(i)), zero))
    return sum((y[i]*y[i]/d[i] for i in range(len(x))), zero)


def correction(d, eta, energy, k, s, r, before, after, zero=F(0), one=F(1)):
    """Accumulate the full six-coordinate action with three factor columns."""
    l, pivots = ldlt(s, zero, one)
    z = [[sum((k[i][v]*l[v][j] for v in range(j, 3)), zero) for j in range(3)] for i in range(6)]
    for i in range(6):
        for j in range(i+1):
            d[i][j] += sum((z[i][v]*pivots[v]*z[j][v] for v in range(3)), zero)
            d[j][i] = d[i][j]
        eta[i] += after[i]-before[i]-sum((k[i][v]*r[v] for v in range(3)), zero)
    return energy+inverse_action(s, r, zero, one)


def analyze(directory, intervals=False):
    import mpmath as mp
    ctx = mp.iv if intervals else mp.mp
    ctx.dps = 40 if intervals else 80
    number = ctx.mpf
    zero, one = number(0), number(1)
    d = [[zero for _ in range(6)] for _ in range(6)]
    eta, energy = [zero]*6, zero
    state = [0.0]*6
    counts = [0, 0]
    eta_bg_prefix = zero
    with (directory/'mean-action.bin').open('rb') as stream:
        while data := stream.read(RECORD.size):
            if len(data) != RECORD.size:
                raise ValueError('truncated construction trace')
            kind, *a = RECORD.unpack(data)
            if kind not in (0, 1):
                raise ValueError('unknown operation')
            counts[kind] += 1
            before, after = (a[1:7], a[7:13]) if kind == 0 else (a[30:36], a[36:42])
            if before != state:
                raise ArithmeticError('unrecorded mean change across a hard event')
            state = after
            before, after = list(map(number, before)), list(map(number, after))
            if kind == 0:
                f = [one]*3+[number(a[0])]*3
                for i in range(6):
                    eta[i] = f[i]*eta[i]+after[i]-f[i]*before[i]
                    for j in range(i+1):
                        d[i][j] *= f[i]*f[j]
                        d[j][i] = d[i][j]
            else:
                k = [[number(a[3*i+j]) for j in range(3)] for i in range(6)]
                s = [[(number(a[18+3*i+j])+number(a[18+3*j+i]))/2 for j in range(3)] for i in range(3)]
                energy = correction(d, eta, energy, k, s, list(map(number, a[27:30])), before, after, zero, one)
            # BG has identity predictor: its roundoff is a sum. Bound every
            # prefix of that defect without replacing the action matrix.
            bound = sum((abs(x) for x in eta[:3]), zero)
            if intervals:
                eta_bg_prefix = number(max(eta_bg_prefix.b, bound.b))
            else:
                eta_bg_prefix = max(eta_bg_prefix, bound)
    native = json.loads((directory/'export.json').read_text())['native']
    target_state = [x[0] for x in native['terminal_state'][3:6]+native['terminal_state'][15:18]]
    if state != target_state:
        raise ArithmeticError('terminal mean mismatch')
    # h_max * (Omega + B_g + N_g + |b_hat_g|) < 3 < pi.
    radius = number(500)-number('0.6108652381980153')-number('.02')-number('.02')-eta_bg_prefix
    barrier = [[(radius*radius if i == j else zero)-energy*d[i][j] for j in range(3)] for i in range(3)]
    _, bg_pivots = ldlt(barrier, zero, one)
    # The actual committed B, not the physical field or a substituted heading.
    field = [number(row[0]) for row in native['committed_field']]
    pivot = max(range(3), key=lambda i: abs(native['committed_field'][i][0]))
    c = []
    for i in range(3):
        if i == pivot:
            continue
        row = [zero]*3
        row[i], row[pivot] = field[pivot], -field[i]
        c.append(row)
    gram = [[sum((c[i][u]*d[3+u][3+v]*c[j][v] for u in range(3) for v in range(3)), zero) for j in range(2)] for i in range(2)]
    grav = [zero, zero, number(float(struct.unpack('<f', struct.pack('<f', 9.80665))[0]))]
    target = [sum((c[i][j]*(grav[j]-eta[3+j]) for j in range(3)), zero) for i in range(2)]
    collinear_cost = inverse_action(gram, target, zero, one)
    def encode(x):
        if intervals:
            return [str(F(*mp.libmp.to_rational(x._mpi_[i]))) for i in (0, 1)]
        return mp.nstr(x, 35)
    return {'arithmetic': 'outward intervals, 40 decimal digits' if intervals else '80 decimal digits, non-promoting',
            'predictions': counts[0], 'accepted_rank_three_corrections': counts[1],
            'innovation_energy': encode(energy), 'joint_action': [[encode(x) for x in row] for row in d],
            'mean_roundoff': [encode(x) for x in eta], 'gyro_roundoff_prefix_l1': encode(eta_bg_prefix),
            'gyro_barrier_LDL_pivots': [encode(x) for x in bg_pivots],
            'collinear_minimum_energy': encode(collinear_cost),
            'collinearity_exclusion_margin': encode(collinear_cost-energy),
            'finite_recorded_gyro_barrier_verified': intervals,
            'collinearity_excluded': bool(collinear_cost > energy),
            'source_uniform_verified': False, 'all_time_magnetic_service_verified': False,
            'theorem_closed': False}


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    p.add_argument('--export', action='store_true')
    p.add_argument('--intervals', action='store_true')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.export:
        export(args.directory, args.eigen)
    result = analyze(args.directory, args.intervals)
    if args.intervals:
        field = [row[0] for row in json.loads((args.directory/'export.json').read_text())['native']['committed_field']]
        result = {'enclosure': result, 'exact_summary': verify_summary(result, field),
                  'recorded_force_separation': recorded_force_separation(args.directory)}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
