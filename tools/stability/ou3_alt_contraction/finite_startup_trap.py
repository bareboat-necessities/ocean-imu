"""A source-controlled antialigned Mahony tail, in the named scalar profile.

This checks invariant inequalities conditional on rounding majorants still to audit.
``verify_prefix`` separately computes that state from the actual first packet;
neither the quaternion nor the integral is installed in the native wrapper.
The controller generates sensor residuals, not estimator writes. Its physical
circle/yaw is the existing docking source. See ou3-alt-startup-trap.md for the
proposed rounding majorants and the distinction between this profile and target firmware.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import ou3_fast_inv_sqrt_interval as INV
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_startup_timeout_alignment_obstruction as DOCK
from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V

QUALIFICATION = 'OU3_ALT_SCALAR_BINARY32_ANTIALIGNED_TAIL_V1'
ROOT_STEP = 30002
U = F(1, 2**24)
ETA = F(1, 2**150)
RADIUS = F(1, 10000)
GAIN = F(16)
MAGNITUDE = F('10.2264552')
DT = SEED.rn(F('.005'))
KP, KI = SEED.rn(F('.2')), SEED.rn(F('.02'))
CFG = V.Config(KP, KI, SEED.rn(F('9.80665')), F(20))


def hex_value(value):
    # Python's hexadecimal conversion is exact for these binary32 inputs.
    return F.from_float(float.fromhex(value))


ROOT_Q = tuple(map(hex_value, ('-0x1.22c4cp-4', '-0x1.312dd6p-1',
                              '-0x1.9370ccp-1', '0x1.fc8a84p-4')))
ROOT_I = tuple(map(hex_value, ('0x1.0fae42p-3', '0x1.d04b2ap-5', '-0x1.cb569ep-4')))


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def down(q):
    """Exact direction of a nonzero stored quaternion; no approximate sqrt."""
    w, x, y, z = q
    s = dot(q, q)
    if s <= 0:
        raise ValueError('nonzero quaternion required')
    return (2*(x*z-w*y)/s, 2*(y*z+w*x)/s, (w*w-x*x-y*y+z*z)/s)


REFERENCE = down(ROOT_Q)


def half_error(q, acc):
    """The same literal binary32 feedback prefix as Mahony_AHRS::update."""
    a = M.Arithmetic()
    add, sub, mul = a.add, a.sub, a.mul
    w, x, y, z = map(M.exact_bits, q)
    ax, ay, az = (M.exact_bits(-v) for v in acc)
    inverse = a.invsqrt(add(add(mul(ax, ax), mul(ay, ay)), mul(az, az)))
    ax, ay, az = (mul(v, inverse) for v in (ax, ay, az))
    vx = sub(mul(x, z), mul(w, y))
    vy = add(mul(w, x), mul(y, z))
    vz = mul(M.exact_bits(F(1, 2)), add(sub(sub(mul(w, w), mul(x, x)), mul(y, y)), mul(z, z)))
    return tuple(map(M.value, (sub(mul(ay, vz), mul(az, vy)),
                              sub(mul(az, vx), mul(ax, vz)),
                              sub(mul(ax, vy), mul(ay, vx)))))


def tail_packet(state):
    """One explicitly specified source, recursively driven by its own shadow.

    The exact construction is followed by the same binary32 API conversion.
    It is independent of any native estimator state installation.
    """
    h = down(state.q)
    acc = tuple(SEED.rn(MAGNITUDE*x) for x in h)
    error = half_error(state.q, acc)
    total = tuple(i+g for i, g in zip(ROOT_I, (F(0), F(0), F('.6'))))
    mu = dot(total, h)
    correction = cross(REFERENCE, h)
    gyro = tuple(SEED.rn(mu*h[j] + GAIN*correction[j] - ROOT_I[j] - KP*error[j])
                 for j in range(3))
    return gyro, acc


@lru_cache(maxsize=None)
def normalization_shell(lo, hi):
    """Exhaust every binary32 norm word in this interval, in 128 bit cells.

    The input is a rounded sum of at most four squares. Correlated division by
    that same sum gives the output norm bounds, including component rounding.
    """
    start = M.round_bits(F(lo)) - 1
    end = M.round_bits(F(hi)) + 1
    count = 128
    width = (end-start+1+count-1)//count
    lower, upper, previous = F(2), F(0), start-1
    for left in range(start, end+1, width):
        right = min(end, left+width-1)
        assert left == previous+1
        x, y = INV._cell_inverse_sqrt(left, right)
        assert y.lo > 0
        low = (1-U)**2/(1+U)**4 * F(x.lo)*F(y.lo)**2 - 1000*ETA
        high = (1+U)**2/(1-U)**4 * F(x.hi)*F(y.hi)**2 + 1000*ETA
        lower, upper, previous = min(lower, low), max(upper, high), right
    assert previous == end
    return lower, upper


def certificate():
    """Rational candidate inequalities; the rounding majorants remain obligations."""
    DOCK.docking_source_bounds()
    root = Path(__file__).resolve().parents[3]
    for name, digest in SENSOR.SOURCES.items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
            raise ValueError('startup source changed; re-audit the invariant: '+name)
    assert dot(REFERENCE, REFERENCE) == 1
    assert F('.99') < dot(ROOT_Q, ROOT_Q) < F('1.01')
    qlo, qhi = normalization_shell(F('.98'), F('1.02'))
    alo, ahi = normalization_shell(F(104), F(106))
    assert F('.99') < qlo <= qhi < F('1.01')
    assert F('.99') < alo <= ahi < F('1.01')
    # Feedback cross products are parallel in exact arithmetic. The actual
    # rounded acceleration/half-gravity/cross graph has norm error <16u.
    feedback = 16*U
    integral_increment = (1+U)**2*KI*DT*feedback + 10*ETA
    half_cells = []
    for value in ROOT_I:
        bits = M.exact_bits(abs(value))
        half_cells.append(min(abs(value)-M.value(bits-1), M.value(bits+1)-abs(value))/2)
    assert integral_increment < min(half_cells)

    # Exact normalized Euler step with w=mu*h+k(h0 cross h):
    # d_next^2/d^2 = 1-(1+c)(dt*k-dt^2*k^2*c/2)/(1+dt^2*|w|^2/4).
    # |mu|<1, |w|<1, c>=1-r^2/2. The deliberately padded .93 contraction
    # absorbs the small axial Euler rotation. 64u covers the complete rounded
    # gyro/add/Euler/normalization direction defect, listed as unaudited obligations in the document.
    total = tuple(i+g for i, g in zip(ROOT_I, (F(0), F(0), F('.6'))))
    assert dot(total, total) + GAIN**2*RADIUS**2 < 1
    cmin = 1-RADIUS**2/2
    ratio2 = 1-(1+cmin)*(DT*GAIN-DT**2*GAIN**2/2)/(1+DT**2/4)
    assert ratio2 < F('.93')**2
    rounding = 64*U
    successor_radius = F('.93')*RADIUS + rounding
    assert successor_radius < RADIUS

    # Physical angular rate remains (0,0,.6). The integral and centre are the
    # attained ones, not freely selected error coordinates.
    mu0 = dot(total, REFERENCE)
    base = tuple(mu0*h-t for h, t in zip(REFERENCE, total))
    assert dot(base, base) < F('.000004')**2
    gyro_residual = F('.000004') + (2+GAIN)*RADIUS + KP*feedback + 2*U
    assert gyro_residual < F('.002') < F('.02')
    # Centre-to-physical packet error is checked against the analytic circle
    # by verify_prefix; the all-time increment about that centre is below .0011.
    accel_increment = MAGNITUDE*RADIUS + U*MAGNITUDE + 2*ETA
    assert accel_increment < F('.0011')
    accel_residual = F('.0071') + accel_increment
    assert accel_residual < DOCK.DOCKING_ACCEL_RESIDUAL
    return {
        'qualification': QUALIFICATION, 'arithmetic_profile': M.PROFILE,
        'root_step': ROOT_STEP, 'radius': RADIUS,
        'exact_Euler_direction_ratio_squared_upper': ratio2,
        'rounded_direction_defect_upper': rounding,
        'successor_radius_upper': successor_radius,
        'invariant_radius_margin': RADIUS-successor_radius,
        'feedback_norm_upper': feedback, 'integral_increment_upper': integral_increment,
        'smallest_integral_rounding_half_cell': min(half_cells),
        'quaternion_norm_squared_enclosure': (qlo, qhi),
        'normalized_accel_norm_squared_enclosure': (alo, ahi),
        'gyro_residual_norm_upper': gyro_residual,
        'accel_residual_norm_upper': accel_residual,
        'conditional_observer_tail_invariant_closed': False,
        'conditional_invariant_inequalities_verified': True,
        'rounding_majorants_audit_closed': False,
        'root_reachability_must_be_verified_separately': True,
        'indefinite_wrapper_guard_composition_closed': False,
        'eventual_finite_startup_refuted': False,
        'target_firmware_correspondence_closed': False,
        'source_uniform_rho_certified': False, 'storage_search_allowed': False,
        'ALT_STARTUP_PASS': False, 'ALT_LIVE_PASS': False, 'ALT_END_TO_END_PASS': False,
    }


def verify_prefix(path, *, tail_steps=0, tail_output=None):
    """Exact reset-to-root calculation plus independent analytic packet audit."""
    certificate()
    lines = Path(path).read_text().splitlines(keepends=True)[:ROOT_STEP]
    if len(lines) != ROOT_STEP:
        raise ValueError('complete reset-to-docking packet prefix required')
    state = V.State()
    operations = 0
    for ordinal, line in enumerate(lines, 1):
        words = line.split()
        if len(words) != 7 or int(words[0]) != ordinal:
            raise ValueError('nonconsecutive source prefix')
        values = tuple(map(hex_value, words[1:]))
        event = SEED.step(state, CFG, dt=DT, gyro=values[:3], acc=values[3:])
        state = event.vertical.state
        operations += len(event.operations)
    if state.q != ROOT_Q or state.integral != ROOT_I:
        raise ValueError('source prefix did not reach the certified docking state')
    # Extend with the new rational source controller. Audit requires at least
    # 30602 packets; this also checks its centre against the actual circle.
    generated = []
    for ordinal in range(ROOT_STEP+1, ROOT_STEP+max(600, tail_steps)+1):
        gyro, acc = tail_packet(state)
        event = M.step_initialized(state, CFG, dt=DT, gyro=gyro, acc=acc)
        state = event.vertical.state
        if state.integral != ROOT_I:
            raise ArithmeticError('integral left the proved rounding deadzone')
        delta = tuple(x-y for x, y in zip(down(state.q), REFERENCE))
        if dot(delta, delta) > RADIUS**2:
            raise ArithmeticError('tail left invariant direction tube')
        generated.append(str(ordinal)+' '+ ' '.join(float(x).hex() for x in (*gyro, *acc))+'\n')
    audit = DOCK.audit_native_packets(iter(lines+generated), docking=True)
    # Exact Taylor enclosures for the fixed body's physical acceleration at
    # the docking time; after docking the physical body packet is constant.
    yaw, _ = DOCK.docking_yaw_and_rate(DOCK.DOCKING_TIME)
    sine, cosine = DOCK._sin_cos_interval(DOCK.DOCKING_OMEGA*DOCK.DOCKING_TIME-yaw)
    physical = ((DOCK.DOCKING_A*cosine[0], DOCK.DOCKING_A*cosine[1]),
                (DOCK.DOCKING_A*sine[0], DOCK.DOCKING_A*sine[1]),
                (-DOCK.SENSOR.G, -DOCK.SENSOR.G))
    centre_error2 = sum(max(abs(MAGNITUDE*h-lo), abs(MAGNITUDE*h-hi))**2
                        for h, (lo, hi) in zip(REFERENCE, physical))
    if centre_error2 >= F('.0071')**2:
        raise ArithmeticError('tail centre detached from the same physical circle')
    if tail_output:
        Path(tail_output).write_text(''.join(lines+generated))
    return {
        'prefix_sha256': hashlib.sha256(''.join(lines).encode()).hexdigest(),
        'root_quaternion': ROOT_Q, 'root_integral': ROOT_I,
        'exact_prefix_operations': operations,
        'reset_to_root_exact_binary32_replay_closed': True,
        'physical_centre_error_norm_squared_upper': centre_error2,
        'finite_packet_source_audit': audit,
        'generated_tail_samples': len(generated),
        'finite_tail_final_quaternion': state.q,
        'initial_state_installed_in_filter': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--tail-output', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = {'tail': certificate(), 'entry': verify_prefix(args.prefix, tail_output=args.tail_output)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, default=str, indent=2, sort_keys=True)+'\n')
    print('Exact scalar-profile docking and conditional tail inequalities verified.')


if __name__ == '__main__':
    main()
