"""Carried-source audit of world-frame rows, injections, service information
and the AW covariance block.

Non-promoting finite diagnostic. It derives a temporary observer from the
existing read-only historical-reader observer and a driver from the existing
driver, compiles an untapped control and requires literal terminal parity.
The shipping sources and committed observers are never edited. 80-digit
values audit algebra on exported float operands; they are not enclosures.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import tempfile

from .ag_readout_source_diagnostic import HEADER, REPO, instrument
from .lin_path_certificate import small_x_source_defect
from .world_frame import nominal_attitude_column_floor

DRIVER = REPO/'tools/stability/ag_readout_source.cpp'
PROFILES = ('0', 'wave', 'collinear', 'collinear-service')
# The driver calls update(.005f, ...): gyro blocks integrate the float32 step.
STEP = struct.unpack('<f', struct.pack('<f', .005))[0]


def observer_source(source):
    """Existing observer plus quaternion and pre-correction state taps."""
    source = instrument(source)
    taps = [
        ('readout_prediction(F_AA, F_LL, Q_AA, Q_LL, trace_phi, (Q_bacc_*trace_qscale).eval());',
         'readout_prediction(F_AA, F_LL, Q_AA, Q_LL, trace_phi, (Q_bacc_*trace_qscale).eval());'
         ' readout_quaternion(qref.coeffs());', 1),
        ('readout_reset(dtheta); }\n',
         'readout_reset(dtheta); readout_quaternion(qref.coeffs()); }\n', 1),
        ('        readout_correction("',
         '        readout_state(qref.coeffs(),xext,Pext.template topLeftCorner<3,3>().eval(),r,S_mat);\n'
         '        readout_correction("', 3),
    ]
    for old, new, count in taps:
        if source.count(old) != count:
            raise ValueError('world-frame observer anchor changed: '+old[:48])
        source = source.replace(old, new)
    return source


MAIN = r'''
template<class A> static void readout_quaternion(const A& q) {
    if (recording) events.push_back("{\"kind\":\"quaternion\",\"q\":"+matrix_json(q)+'}');
}
template<class A, class B, class C, class D, class E>
static void readout_state(const A& q, const B& x, const C& ptt, const D& r, const E& s) {
    if (recording) events.push_back("{\"kind\":\"state\",\"q\":"+matrix_json(q)+",\"x\":"
        +matrix_json(x)+",\"P_theta\":"+matrix_json(ptt)+",\"r\":"+matrix_json(r)
        +",\"S\":"+matrix_json(s)+'}');
}
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

int main(int argc, char** argv) {
    const std::string profile = argc > 1 ? argv[1] : "0";
    const bool wave = profile == "wave", collinear = profile == "collinear";
    // Same collinear motion with 25-Hz magnetic corrections throughout.
    const bool service = profile == "collinear-service", moving = collinear || service;
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.12f);
    cfg.sigma_g.setConstant(.00135f);
    cfg.sigma_m.setConstant(.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    // Collinear witness: a(t)=c cos(2 pi t), c=g(e_z-(e_z.b)b), b=(7,0,24)/25.
    const double cx = -168.0/625.0*9.80665, cz = 49.0/625.0*9.80665;
    const int start = 45001, stop = moving ? 45800 : 45064;
    int live=-1, refined=-1, active=-1, applied=0;
    std::string root;
    for (int k=1; k<=stop; ++k) {
        if (k==start) {
            if (active<0 || k-active<3400) return 2;
            const auto& m=filter.raw().mekf();
            root="\"root_covariance\":"+matrix_json(m.Pext)+",\"root_state\":"+matrix_json(m.xext)
                +",\"root_quaternion\":"+matrix_json(m.qref.coeffs())+",\"reference\":"+matrix_json(m.v2ref)
                +",\"gravity\":"+std::to_string(static_cast<double>(m.gravity_magnitude_));
            recording=true;
        }
        const double t = static_cast<double>(k)*.005;
        float roll=0.0f, rate=0.0f;
        Eigen::Vector3f force(0.0f,0.0f,-g_std), field(75.0f,0.0f,0.0f);
        if (wave) {
            roll = static_cast<float>(.02*std::sin(.5*t));
            rate = static_cast<float>(.01*std::cos(.5*t));
            force.z() = static_cast<float>(-.144*std::sin(.6*t))-g_std;
            field = Eigen::Vector3f(60.0f,0.0f,30.0f);
        } else if (moving) {
            roll = static_cast<float>(.01*std::sin(.5*t));
            rate = static_cast<float>(.005*std::cos(.5*t));
            const double c = std::cos(2.0*3.14159265358979323846*t);
            force = Eigen::Vector3f(static_cast<float>(cx*c),0.0f,static_cast<float>(cz*c)-g_std);
            field = Eigen::Vector3f(21.0f,0.0f,72.0f);
        }
        const Eigen::Matrix3f rwb=Eigen::AngleAxisf(-roll,Eigen::Vector3f::UnitX()).toRotationMatrix();
        filter.update(.005f,Eigen::Vector3f(rate,0,0),rwb*force);
        const bool due = collinear && recording ? k%200==0 : k%8==0;
        if (due) {
            filter.updateMag((rwb*field).eval());
            if (recording && filter.raw().mekf().lastMagDiag().accepted) ++applied;
        }
        if (live<0 && filter.isLive()) live=k;
        if (refined<0 && filter.hasRefinedMagReference()) refined=k;
        if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;
    }
    const auto& m=filter.raw().mekf();
    if (!m.Pext.allFinite() || !m.xext.allFinite() || applied!=(collinear ? 4 : service ? 100 : 8)) return 3;
    std::cout << "{\"live_step\":" << live << ",\"refined_step\":" << refined
              << ",\"active_step\":" << active << ",\"applied_magnetic_updates\":" << applied
              << ',' << root
              << ",\"terminal_covariance\":" << matrix_json(m.Pext)
              << ",\"terminal_state\":" << matrix_json(m.xext)
              << ",\"terminal_quaternion\":" << matrix_json(m.qref.coeffs())
              << ",\"events\":[";
    for (std::size_t i=0; i<events.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << events[i];
    }
    std::cout << "]}\n";
}
'''


def driver_source():
    """Existing driver helpers with extra taps and the three input profiles."""
    source = DRIVER.read_text()
    anchor = '#define private public'
    if source.count(anchor) != 1:
        raise ValueError('readout driver anchor changed')
    return source[:source.index(anchor)]+MAIN.lstrip('\n')


def _truth(profile, time):
    """World acceleration of the driver input (NED, heading of its true field)."""
    if profile == 'wave':
        return (0.0, 0.0, -0.144*math.sin(0.6*time))
    if profile in ('collinear', 'collinear-service'):
        c = 9.80665*math.cos(2*math.pi*time)
        return (-168/625*c, 0.0, 49/625*c)
    return (0.0, 0.0, 0.0)


def analyze(trace, profile, dps=80):
    """World rows, same-cell floors, injection budgets, aggregate rank and AW block."""
    import mpmath as mp
    with mp.workdps(dps):
        num = lambda x: mp.mpf(float(x))
        mat = lambda a: mp.matrix([[num(x) for x in row] for row in a])
        col = lambda a: mp.matrix([num(row[0]) for row in a])
        eye = mp.eye(3)

        def rotation(q):
            x, y, z, w = (num(v[0]) for v in q)
            n = w*w+x*x+y*y+z*z
            return mp.matrix([[(w*w+x*x-y*y-z*z)/n, 2*(x*y-w*z)/n, 2*(x*z+w*y)/n],
                              [2*(x*y+w*z)/n, (w*w-x*x+y*y-z*z)/n, 2*(y*z-w*x)/n],
                              [2*(x*z-w*y)/n, 2*(y*z+w*x)/n, (w*w-x*x-y*y+z*z)/n]])

        def skew(v):
            return mp.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

        def unskew(m):
            return mp.matrix([m[2, 1], m[0, 2], m[1, 0]])

        def cross(u, v):
            return mp.matrix([u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]])

        def hstack(a, b):
            return mp.matrix([[a[r, c] for c in range(a.cols)]+[b[r, c] for c in range(b.cols)]
                              for r in range(a.rows)])

        def vstack(a, b):
            return mp.matrix([[a[r, c] for c in range(a.cols)] for r in range(a.rows)]
                             + [[b[r, c] for c in range(b.cols)] for r in range(b.rows)])

        def smin(m):
            return min(mp.svd(m, compute_uv=False))

        def lam_minus(ff, bb, kappa):
            return (ff+bb)/2-mp.sqrt((ff-bb)**2/4+ff*bb*kappa**2)

        def exp_rotation(v):
            theta = mp.norm(v)
            if theta == 0:
                return mp.eye(3)
            k = skew(v/theta)
            return mp.eye(3)+mp.sin(theta)*k+(1-mp.cos(theta))*k*k

        def rotation_angle(m):
            return mp.atan2(mp.norm(unskew((m-m.T)/2)), (m[0, 0]+m[1, 1]+m[2, 2]-1)/2)

        def psi(theta):
            rho = 1/mp.sqrt(1+theta**2/4)
            return theta-mp.atan(theta/2)+mp.atan((1-rho)/(2*mp.sqrt(rho)))

        g_c = num(trace['gravity'])
        reference = col(trace['reference'])
        # The estimator world heading is that of its committed reference; the
        # driver truth is rotated about e_z into that frame, never re-tilted.
        heading = mp.atan2(reference[1], reference[0])
        spin = mp.matrix([[mp.cos(heading), -mp.sin(heading), 0],
                          [mp.sin(heading), mp.cos(heading), 0], [0, 0, 1]])
        r0 = rotation(trace['root_quaternion'])
        rot = r0
        body_a, body_b = mp.eye(3), mp.zeros(3)
        world_a, world_b = mp.eye(3), mp.zeros(3)
        pred_defect = reset_defect = row_defect = force_defect = field_defect = mp.mpf(0)
        injection_float_angle = mp.mpf(0)
        max_rate = mp.mpf(0)
        steps, rows, anchor, groups, injections, ops, aw_errors = 0, [], None, [], [], [], []
        signed_world_injection, injection_norm_sum = mp.zeros(3, 1), mp.mpf(0)

        def service_root(rotation_now):
            # Normalized heading (1 rad) and axial gyro bias (.02 rad/s) about
            # the world vertical, in body error coordinates at the window root.
            axis = rotation_now*mp.matrix([0, 0, 1])
            e = mp.zeros(21, 2)
            for j in range(3):
                e[j, 0], e[3+j, 1] = axis[j], axis[j]/50
            return e
        service = service_root(r0)
        service_information = []

        # AW marginal: 3x3 blocks only, P_aw at offset 15 and a_w at LIN 9.
        def block(matrix, a, b):
            return mp.matrix([[num(x) for x in row[a:b]] for row in matrix[a:b]])

        def spectrum(m):
            return mp.eigsy((m+m.T)/2, eigvals_only=True)

        def positive_part(m):
            values, vectors = mp.eigsy((m+m.T)/2)
            return vectors*mp.diag([max(0, x) for x in values])*vectors.T
        eps_q = small_x_source_defect()[0]
        one_eps = 1+mp.mpf(eps_q.numerator)/eps_q.denominator
        threshold = nominal_attitude_column_floor(16, 0, Fraction(1, 5))['aw_error_threshold']
        aw_now = block(trace['root_covariance'], 15, 18)
        aw_last = aw_now
        aw_seed = max(spectrum(aw_now))
        aw_step = aw_delta = None
        awc = {'recon': mp.mpf(0), 'iso': mp.mpf(0), 'floor': mp.inf, 'sigma2': mp.mpf(0),
               'ratio': mp.mpf(0), 'ceiling': mp.mpf(0), 'peak': mp.mpf(0), 'trough': mp.inf,
               'storage': mp.mpf(0), 'syncs': 0}
        pending = None
        events = trace['events']
        i = 0
        while i < len(events):
            e = events[i]
            kind = e['kind']
            if kind in ('prediction', 'reset'):
                if i+1 >= len(events) or events[i+1]['kind'] != 'quaternion':
                    raise ValueError('quaternion tap missing after '+kind)
                nxt = rotation(events[i+1]['q'])
                if kind == 'prediction':
                    lin = e['F_LIN']
                    if any(lin[9+a][j] != 0 for a in range(3) for j in range(12) if j != 9+a):
                        raise ValueError('AW prediction row is not phi e_a')
                    phi = mp.diag([num(lin[9+a][9+a]) for a in range(3)])
                    q_aa = block(e['Q_LIN'], 9, 12)
                    sigma2_q = max(q_aa[a, a]/(1-phi[a, a]**2) for a in range(3))
                    aw_step, aw_delta = (phi, q_aa, sigma2_q), None
                    f = mat(e['F_AG'])
                    rs, bs = f[:3, :3], f[:3, 3:]
                    w = nxt.T*rs*rot
                    pred_defect = max(pred_defect, mp.mnorm(w-eye, 1))
                    body_a, body_b = rs*body_a, rs*body_b+bs
                    world_a, world_b = w*world_a, w*world_b+nxt.T*bs
                    full = mp.eye(21)
                    full[:6, :6], full[6:18, 6:18] = f, mat(e['F_LIN'])
                    for j in range(18, 21):
                        full[j, j] = num(e['phi_BA'])
                    service = full*service
                    # atan2 keeps float32 angles whose cosine rounds to one.
                    axial = mp.norm(unskew((rs-rs.T)/2))
                    angle = mp.atan2(axial, (rs[0, 0]+rs[1, 1]+rs[2, 2]-1)/2)
                    max_rate = max(max_rate, angle/num(STEP))
                    ops.append(('predict', angle, None))
                    steps, anchor = steps+1, None
                else:
                    d = col(e['d'])
                    gmat = eye+skew(d)/2
                    n = nxt.T*gmat*rot
                    x = rot.T*d
                    signed_world_injection += x
                    injection_norm_sum += mp.norm(d)
                    expected = eye+((x.T*x)[0]*eye-x*x.T)/4
                    reset_defect = max(reset_defect, mp.mnorm(n.T*n-expected, 1))
                    # psi(theta) bounds the turning of the exact injection
                    # Exp([d]x); the exported float mean injection differs by
                    # the rotation angle delta, charged by the triangle inequality.
                    delta = rotation_angle(nxt*rot.T*exp_rotation(d).T)
                    injection_float_angle = max(injection_float_angle, delta)
                    body_a, body_b = gmat*body_a, gmat*body_b
                    world_a, world_b = n*world_a, n*world_b
                    service[:3, :] = gmat*service[:3, :]
                    if pending is not None and pending[4] == 'mag':
                        # A new one-correction window starts after this reset.
                        service = service_root(nxt)
                    bound = None
                    if pending is not None:
                        k_theta, s, r, p_theta, sensor = pending
                        nis = (r.T*mp.inverse(s)*r)[0]
                        ksk = k_theta*s*k_theta.T
                        gain_action = max(mp.eigsy((ksk+ksk.T)/2, eigvals_only=True))
                        prior = max(mp.eigsy((p_theta+p_theta.T)/2, eigvals_only=True))
                        bound = mp.sqrt(nis*prior)
                        # Zero innovation gives zero injection exactly; 0/0 is
                        # recorded as ratio zero only when the injection is zero.
                        ratio = lambda num_, den: (num_/den if den > 0 else
                                                   (mp.mpf(0) if num_ == 0 else mp.inf))
                        injections.append({'sensor': sensor, 'norm': mp.norm(d),
                                           'consistency': mp.norm(k_theta*r-d),
                                           'gain_ratio': ratio(mp.norm(d), mp.sqrt(nis*gain_action)),
                                           'prior_ratio': ratio(mp.norm(d), bound),
                                           'prior_bound': bound, 'NIS': nis})
                        pending = None
                    ops.append(('reset', mp.norm(d), bound))
                    if anchor is not None:
                        anchor['angle'] += psi(mp.norm(d))+delta
                        anchor['local'] = gmat*anchor['local']
                        anchor['world_local'] = n*anchor['world_local']
                rot = nxt
                i += 2
                continue
            if kind == 'state':
                correction = events[i+1]
                if correction['kind'] != 'correction':
                    raise ValueError('state tap must precede its correction')
                if mp.mnorm(rotation(e['q'])-rot, 1) > mp.mpf(10)**-60:
                    raise ValueError('pre-correction quaternion differs from carried rotation')
                sensor = correction['sensor']
                gain, s_act = mat(correction['K']), mat(e['S'])
                trough = max(spectrum(aw_now))
                pending = (gain[:3, :], s_act, col(e['r']), mat(e['P_theta']), sensor)
                h_full = mat(correction['H'])
                if sensor == 'mag':
                    g_rows = mp.inverse(mp.cholesky((s_act+s_act.T)/2))*h_full*service
                    info = g_rows.T*g_rows
                    service_information.append(min(mp.eigsy((info+info.T)/2, eigvals_only=True)))
                service = (mp.eye(21)-gain*h_full)*service
                # K=PC'S^-1 on the AW rows: the Joseph AW block drops by K_a S K_a'.
                k_aw = gain[15:18, :]
                aw_now = aw_now-k_aw*s_act*k_aw.T
                if sensor in ('acc', 'mag'):
                    h = mat(correction['H'])[:, :6]
                    if mp.mnorm(h[:, 3:], 1) != 0:
                        raise ValueError('zero-lever-arm gyro columns required')
                    ha = h[:, :3]
                    body_rows = ha*hstack(body_a, body_b)
                    world_vec = rot.T*unskew(-ha)
                    world_rows = -rot*skew(world_vec)*hstack(world_a*r0.T, world_b)
                    row_defect = max(row_defect, mp.mnorm(body_rows-world_rows, 1))
                    rows.append(body_rows)
                    if sensor == 'acc':
                        xs = [num(v[0]) for v in e['x']]
                        aw = mp.matrix(xs[15:18])
                        force = aw-g_c*mp.matrix([0, 0, 1])
                        force_defect = max(force_defect, mp.norm(force-world_vec))
                        # Driver time is the double product k*.005.
                        time = (45000+steps)*.005
                        truth = spin*mp.matrix([num(v) for v in _truth(profile, time)])
                        aw_errors.append(mp.norm(aw-truth))
                        awc['trough'] = min(awc['trough'], trough)
                        awc['storage'] = max(awc['storage'], aw_errors[-1]**2/trough)
                        anchor = {'C_acc': ha, 'local': mp.eye(3), 'world_local': mp.eye(3),
                                  'force': world_vec, 'angle': mp.mpf(0), 'op': len(ops),
                                  'transport': (body_a.copy(), body_b.copy())}
                    else:
                        field_defect = max(field_defect, mp.norm(world_vec-reference))
                        if anchor is not None:
                            # World vectors are read from the literal rows, so the
                            # comparison uses exactly the exported operands.
                            c = vstack(anchor['C_acc'], ha*anchor['local'])
                            f, b, n = anchor['force'], world_vec, anchor['world_local']
                            ff, bb = (f.T*f)[0], (b.T*b)[0]
                            sine = min(1, mp.norm(cross(f, b))/mp.sqrt(ff*bb))
                            k = mp.inverse(n)*b
                            turned = mp.acos(min(1, abs((k.T*b)[0])/(mp.norm(k)*mp.norm(b))))
                            kappa = mp.cos(max(0, mp.asin(sine)-anchor['angle']))
                            groups.append({'actual': smin(c), 'C': c, 'sine': sine, 'op': anchor['op'],
                                           'floor': mp.sqrt(max(0, lam_minus(ff, bb, kappa))),
                                           'reset_angle_slack': anchor['angle']-turned,
                                           'transport': anchor['transport']})
                i += 1
                continue
            if kind == 'sync':
                aw_delta = mat(e['Q'])
            elif kind == 'sync_completion':
                before, after = block(e['before'], 15, 18), block(e['after'], 15, 18)
                phi, q_aa, sigma2_q = aw_step
                # Corrections since the last completion were charged K_a S K_a'.
                predicted = phi*aw_now*phi.T+q_aa
                awc['recon'] = max(awc['recon'], mp.mnorm(predicted-before, 1)/mp.mnorm(before, 1))
                target = mp.mpf(0)
                if aw_delta is not None:
                    awc['syncs'] += 1
                    if mp.mnorm(aw_delta, 1) > 0:
                        # Isotropic target: after=sigma^2 I+(before-sigma^2 I)_+.
                        target = min(spectrum(after))
                        iso = target*eye+positive_part(before-target*eye)
                        awc['iso'] = max(awc['iso'], mp.mnorm(after-iso, 1)/mp.mnorm(after, 1))
                    awc['floor'] = min(awc['floor'], min(spectrum(after))/sigma2_q)
                peak = max(spectrum(after))
                bound = max(phi[a, a]**2 for a in range(3))*max(spectrum(aw_last))+max(spectrum(q_aa))
                awc['ratio'] = max(awc['ratio'], peak/max(bound, target))
                awc['sigma2'] = max(awc['sigma2'], sigma2_q, target)
                awc['ceiling'] = max(awc['ceiling'], peak/max(aw_seed, one_eps*awc['sigma2']))
                awc['peak'] = max(awc['peak'], peak)
                aw_now = aw_last = after
            elif kind not in ('correction', 'quaternion'):
                raise ValueError('unexpected event '+kind)
            if kind in ('correction', 'sync', 'sync_completion', 'quaternion'):
                i += 1
                continue
        if len(groups) < 2 or not injections:
            raise ValueError('two same-cell groups and applied injections required')
        first, last = groups[0], groups[-1]
        a0, b0 = first['transport']
        a1, b1 = last['transport']
        between_a = a1*mp.inverse(a0)
        between_b = b1-between_a*b0
        actual_rows = vstack(hstack(first['C'], mp.zeros(6, 3)),
                             hstack(last['C']*between_a, last['C']*between_b))
        c_floor = min(first['floor'], last['floor'])
        attitude_norm = max(mp.svd(between_a, compute_uv=False))
        # Existing body-frame inverse recurrence, with each reset charged by its
        # NIS/prior-covariance bound instead of its realized norm.
        beta = defect = elapsed = mp.mpf(0)
        between = ops[first['op']:last['op']]
        injection_actual = sum((v for kind, v, b in between if kind == 'reset'), mp.mpf(0))
        injection_bound = sum((b for kind, v, b in between if kind == 'reset' and b is not None), mp.mpf(0))
        for kind, value, bound in between:
            if kind == 'predict':
                h = num(STEP)
                defect += h*beta+h*(value/2+value**2/3)
                elapsed += h
                beta = min(2, beta+value+value**2/2)
            else:
                beta = min(2, beta+min(1, (bound if bound is not None else value)/2))
        gyro_floor = max(0, elapsed-defect)
        s_world = c_floor/(1+(attitude_norm+1)/gyro_floor) if gyro_floor > 0 else mp.mpf(0)
        gram = mp.zeros(6)
        for row in rows:
            gram += row.T*row
        aggregate = mp.sqrt(max(0, min(mp.eigsy((gram+gram.T)/2, eigvals_only=True))))
        attitude = mp.sqrt(max(0, min(mp.eigsy((gram[:3, :3]+gram[:3, :3].T)/2, eigvals_only=True))))
        fmt = lambda v: mp.nstr(v, 24)
        return {
            'applied_rows': 3*len(rows), 'same_cell_groups': len(groups),
            'prediction_world_discrepancy_max': fmt(pred_defect),
            'reset_gram_identity_defect_max': fmt(reset_defect),
            'reset_injection_float_angle_max': fmt(injection_float_angle),
            'row_factorization_defect_max': fmt(row_defect),
            'acc_force_state_vs_row_defect_max': fmt(force_defect),
            'mag_row_vs_reference_defect_max': fmt(field_defect),
            'maximum_nominal_rotation_rate_rad_s': fmt(max_rate),
            'same_cell_actual_singular_min': fmt(min(gr['actual'] for gr in groups)),
            'same_cell_actual_singular_max': fmt(max(gr['actual'] for gr in groups)),
            'same_cell_world_floor_min': fmt(min(gr['floor'] for gr in groups)),
            'same_cell_floor_to_actual_ratio_max': fmt(max(gr['floor']/gr['actual'] for gr in groups)),
            'same_cell_force_field_sine_min': fmt(min(gr['sine'] for gr in groups)),
            'same_cell_reset_angle_bound_slack_min': fmt(min(gr['reset_angle_slack'] for gr in groups)),
            'injections': len(injections),
            'injection_norm_max_rad': fmt(max(inj['norm'] for inj in injections)),
            'injection_gain_consistency_max': fmt(max(inj['consistency'] for inj in injections)),
            'injection_to_gain_action_bound_ratio_max': fmt(max(inj['gain_ratio'] for inj in injections)),
            'injection_to_prior_covariance_bound_ratio_max': fmt(max(inj['prior_ratio'] for inj in injections)),
            'injection_prior_covariance_bound_max_rad': fmt(max(inj['prior_bound'] for inj in injections)),
            'NIS_max': fmt(max(inj['NIS'] for inj in injections)),
            'two_group_world_c_floor': fmt(c_floor),
            'interanchor_attitude_norm': fmt(attitude_norm),
            'interanchor_gyro_actual_singular': fmt(smin(between_b)),
            'interanchor_gyro_NIS_budget_floor': fmt(gyro_floor),
            'interanchor_injection_sum_actual_rad': fmt(injection_actual),
            'interanchor_injection_sum_NIS_prior_bound_rad': fmt(injection_bound),
            'interanchor_seconds': fmt(elapsed),
            'two_group_budget_from_world_and_NIS': fmt(s_world),
            'two_group_actual_singular': fmt(smin(actual_rows)),
            'aggregate_six_column_singular': fmt(aggregate),
            'aggregate_attitude_column_singular': fmt(attitude),
            'window_world_attitude_transport_distance': fmt(max(mp.svd(world_a-eye, compute_uv=False))),
            'window_signed_world_injection_norm_rad': fmt(mp.norm(signed_world_injection)),
            'window_injection_norm_sum_rad': fmt(injection_norm_sum),
            'aw_tracking_error_max_mps2': fmt(max(aw_errors)),
            'one_correction_window_service_lambda_min_max': fmt(max(service_information)),
            'aw_syncs': awc['syncs'],
            'aw_prediction_reconstruction_defect_max': fmt(awc['recon']),
            'aw_sync_isotropy_defect_max': fmt(awc['iso']),
            'aw_post_sync_floor_ratio_min': fmt(awc['floor']),
            'aw_stationary_variance_max': fmt(awc['sigma2']),
            'aw_lemma_step_ratio_max': fmt(awc['ratio']),
            'aw_ceiling_ratio_max': fmt(awc['ceiling']),
            'aw_covariance_lambda_max_max': fmt(awc['peak']),
            'aw_covariance_lambda_max_min_at_acc': fmt(awc['trough']),
            'aw_storage_lower_bound_max': fmt(awc['storage']),
            # Uniform storage route: sup V * sup lambda_max(P_aw) < threshold^2.
            'aw_uniform_storage_route_ratio': fmt(awc['storage']*awc['peak']/(mp.mpf(threshold.numerator)/threshold.denominator)**2),
            'source_uniform_verified': False, 'real_trajectory_enclosed': False,
        }


METRICS = (
    'prediction_world_discrepancy_max', 'reset_gram_identity_defect_max', 'reset_injection_float_angle_max',
    'row_factorization_defect_max', 'acc_force_state_vs_row_defect_max',
    'mag_row_vs_reference_defect_max', 'maximum_nominal_rotation_rate_rad_s',
    'same_cell_actual_singular_min', 'same_cell_actual_singular_max', 'same_cell_world_floor_min',
    'same_cell_floor_to_actual_ratio_max', 'same_cell_force_field_sine_min',
    'same_cell_reset_angle_bound_slack_min', 'injection_norm_max_rad',
    'injection_gain_consistency_max', 'injection_to_gain_action_bound_ratio_max',
    'injection_to_prior_covariance_bound_ratio_max', 'injection_prior_covariance_bound_max_rad',
    'NIS_max', 'two_group_world_c_floor', 'interanchor_attitude_norm',
    'interanchor_gyro_actual_singular', 'interanchor_gyro_NIS_budget_floor',
    'interanchor_injection_sum_actual_rad', 'interanchor_injection_sum_NIS_prior_bound_rad',
    'interanchor_seconds', 'two_group_budget_from_world_and_NIS', 'two_group_actual_singular',
    'aggregate_six_column_singular', 'aggregate_attitude_column_singular',
    'window_world_attitude_transport_distance', 'window_signed_world_injection_norm_rad',
    'window_injection_norm_sum_rad', 'aw_tracking_error_max_mps2',
    'one_correction_window_service_lambda_min_max', 'aw_prediction_reconstruction_defect_max',
    'aw_sync_isotropy_defect_max', 'aw_post_sync_floor_ratio_min', 'aw_stationary_variance_max',
    'aw_lemma_step_ratio_max', 'aw_ceiling_ratio_max', 'aw_covariance_lambda_max_max',
    'aw_covariance_lambda_max_min_at_acc', 'aw_storage_lower_bound_max', 'aw_uniform_storage_route_ratio')
COUNTS = ('applied_rows', 'same_cell_groups', 'injections', 'live_step', 'refined_step', 'active_step',
          'aw_syncs')
CASE_KEYS = frozenset(METRICS+COUNTS+(
    'input_profile', 'literal_terminal_parity', 'reference', 'exported_trace_sha256',
    'source_uniform_verified', 'real_trajectory_enclosed'))
SUMMARY_KEYS = frozenset((
    'qualification', 'shipping_header_sha256', 'observer_header_sha256', 'driver_sha256',
    'decimal_digits', 'cases', 'source_uniform_verified', 'all_time_magnetic_service_certified',
    'theorem_closed'))
# Corollary A: nominal attitude columns stay positive iff epsilon_a is below this.
AW_TRACKING_THRESHOLD = 1.12383


def _require(condition, message):
    if not condition:
        raise ValueError('world-frame diagnostic: '+message)


def verify_diagnostic(summary):
    """Verify provenance, scope and every reported metric of the committed record.

    The unretained trace is bound by exact reproduction (``--expect``); this
    check validates structure, the relations any correct analysis satisfies,
    identity tolerances and the finite conclusions the documentation draws.
    """
    source = (REPO/HEADER).read_text()
    expected = {'shipping_header_sha256': hashlib.sha256(source.encode()).hexdigest(),
                'observer_header_sha256': hashlib.sha256(observer_source(source).encode()).hexdigest(),
                'driver_sha256': hashlib.sha256(driver_source().encode()).hexdigest()}
    if any(summary.get(k) != v for k, v in expected.items()):
        raise ValueError('world-frame source/observer/driver provenance changed')
    if summary.get('qualification') != 'OU3_CARRIED_WORLD_FRAME_DIAGNOSTIC_V1':
        raise ValueError('world-frame diagnostic qualification changed')
    for flag in ('source_uniform_verified', 'all_time_magnetic_service_certified', 'theorem_closed'):
        if summary.get(flag) is not False:
            raise ValueError('finite world-frame diagnostic cannot promote '+flag)
    _require(set(summary) == SUMMARY_KEYS and summary['decimal_digits'] == 80, 'summary fields changed')
    cases = summary.get('cases', [])
    if [c.get('input_profile') for c in cases] != list(PROFILES):
        raise ValueError('all three carried source profiles required')
    for case in cases:
        if (case.get('literal_terminal_parity') is not True or
                case.get('source_uniform_verified') is not False or
                case.get('real_trajectory_enclosed') is not False):
            raise ValueError('world-frame source audit scope/parity failed')
        _require(set(case) == CASE_KEYS, 'case fields changed')
        for key in METRICS:
            _require(isinstance(case[key], str) and math.isfinite(float(case[key])), 'non-finite '+key)
        v = {key: float(case[key]) for key in METRICS}
        n = {key: case[key] for key in COUNTS}
        _require(all(type(x) is int and x > 0 for x in n.values()), 'counts must be positive integers')
        _require(n['applied_rows'] % 3 == 0 and n['same_cell_groups'] >= 2, 'row/group counts')
        _require(n['live_step'] <= n['refined_step'] <= n['active_step'], 'regime step order')
        _require(len(case['reference']) == 3 and all(math.isfinite(x) for x in case['reference']),
                 'committed reference')
        _require(len(case['exported_trace_sha256']) == 64
                 and set(case['exported_trace_sha256']) <= set('0123456789abcdef'), 'trace hash')
        nonnegative = [k for k in METRICS if k != 'one_correction_window_service_lambda_min_max']
        _require(all(v[k] >= 0 for k in nonnegative), 'negative norm, ratio or slack')
        _require(v['one_correction_window_service_lambda_min_max'] >= -1e-60, 'indefinite service Gram')
        # Exact identities hold to 80-digit rounding; float operands to float rounding.
        _require(v['row_factorization_defect_max'] <= 1e-60, 'row factorization defect')
        _require(v['reset_gram_identity_defect_max'] <= 1e-60, 'reset Gram identity defect')
        _require(v['prediction_world_discrepancy_max'] <= 1e-6
                 and v['reset_injection_float_angle_max'] <= 1e-6, 'float rotation discrepancy')
        _require(v['injection_gain_consistency_max'] <= 1e-9, 'injection differs from K_theta r')
        _require(v['acc_force_state_vs_row_defect_max'] <= 1e-4
                 and v['mag_row_vs_reference_defect_max'] <= 1e-4, 'state/reference row defect')
        # Relations any correct analysis of any trace satisfies.
        for key in ('same_cell_floor_to_actual_ratio_max', 'injection_to_gain_action_bound_ratio_max',
                    'injection_to_prior_covariance_bound_ratio_max'):
            if not v[key] <= 1:
                raise ValueError('conditional budget exceeded its actual value: '+key)
        _require(v['same_cell_world_floor_min'] <= v['same_cell_actual_singular_min']
                 <= v['same_cell_actual_singular_max'], 'same-cell floor/actual order')
        _require(v['same_cell_world_floor_min'] <= v['two_group_world_c_floor'], 'two-group c floor')
        _require(v['same_cell_force_field_sine_min'] <= 1, 'sine above one')
        _require(v['two_group_budget_from_world_and_NIS'] <= v['two_group_actual_singular'],
                 'two-group budget above actual')
        _require(v['interanchor_gyro_NIS_budget_floor'] <= v['interanchor_gyro_actual_singular'],
                 'gyro budget above actual')
        _require(v['interanchor_injection_sum_actual_rad'] <= v['interanchor_injection_sum_NIS_prior_bound_rad']
                 and v['interanchor_injection_sum_actual_rad'] <= v['window_injection_norm_sum_rad'],
                 'injection sums')
        _require(v['injection_norm_max_rad'] <= v['injection_prior_covariance_bound_max_rad'],
                 'injection above its bound')
        _require(v['window_signed_world_injection_norm_rad'] <= v['window_injection_norm_sum_rad'],
                 'signed injection sum above norm sum')
        # Cauchy interlacing: a principal block has the larger least eigenvalue.
        _require(v['aggregate_six_column_singular'] <= v['aggregate_attitude_column_singular'],
                 'aggregate attitude block below six-column value')
        _require(v['interanchor_attitude_norm'] >= 1-1e-3 and v['interanchor_seconds'] > 0,
                 'inter-anchor transport')
        _require(v['aw_tracking_error_max_mps2'] < AW_TRACKING_THRESHOLD, 'AW tracking threshold')
        # AW ceiling lemma on the literal blocks: float-level reconstruction,
        # isotropic syncs flooring P_aw at sigma^2, and step/ceiling ratios <=1.
        _require(v['aw_prediction_reconstruction_defect_max'] <= 1e-5
                 and v['aw_sync_isotropy_defect_max'] <= 1e-5, 'AW block reconstruction')
        _require(v['aw_lemma_step_ratio_max'] <= 1+1e-5 and v['aw_ceiling_ratio_max'] <= 1+1e-5,
                 'AW ceiling lemma ratio')
        _require(v['aw_post_sync_floor_ratio_min'] >= 1-1e-4 and n['aw_syncs'] >= 3, 'AW sync floor')
        _require(.0025*(1-1e-3) <= v['aw_stationary_variance_max'] <= 16*(1+1e-5)
                 and v['aw_covariance_lambda_max_min_at_acc'] <= v['aw_covariance_lambda_max_max'] <= 16.001,
                 'AW variance range')
        route = v['aw_storage_lower_bound_max']*v['aw_covariance_lambda_max_max']/AW_TRACKING_THRESHOLD**2
        _require(abs(route-v['aw_uniform_storage_route_ratio']) <= 1e-12*max(1, route), 'storage route ratio')
        profile = case['input_profile']
        if profile == '0':
            _require(v['aw_storage_lower_bound_max'] == v['aw_uniform_storage_route_ratio'] == 0,
                     'quiet AW storage')
            _require(v['NIS_max'] == v['injection_norm_max_rad'] == v['aw_tracking_error_max_mps2'] == 0
                     and v['same_cell_floor_to_actual_ratio_max'] == 1
                     and v['interanchor_gyro_NIS_budget_floor'] == v['interanchor_gyro_actual_singular'],
                     'quiet profile conclusions')
        elif profile == 'wave':
            _require(v['same_cell_actual_singular_min'] > 1 and v['aw_uniform_storage_route_ratio'] < 1
                     and 1 < v['two_group_budget_from_world_and_NIS'], 'wave profile conclusions')
        elif profile == 'collinear-service':
            # Same motion under 25-Hz magnetic corrections: the uniform storage
            # route to Corollary A fails although the actual AW error passes.
            _require(v['aw_uniform_storage_route_ratio'] > 1, 'collinear storage route conclusions')
        else:
            # Degenerate same-cell groups, well-conditioned aggregate rows, a
            # failed one-correction service (mu_M=1) and a failed NIS budget.
            _require(v['same_cell_actual_singular_max'] < 1 < 10 < v['aggregate_six_column_singular'],
                     'collinear same-cell/aggregate conclusions')
            _require(v['one_correction_window_service_lambda_min_max'] < 1e-3, 'collinear service')
            _require(v['interanchor_gyro_NIS_budget_floor'] == v['two_group_budget_from_world_and_NIS'] == 0
                     and v['interanchor_gyro_actual_singular'] > 1
                     and v['interanchor_injection_sum_NIS_prior_bound_rad']
                     > 100*v['interanchor_injection_sum_actual_rad'], 'collinear NIS budget conclusions')
            _require(v['aw_uniform_storage_route_ratio'] > 1, 'collinear storage route conclusions')
            half = v['window_signed_world_injection_norm_rad']/2
            _require(10*v['window_signed_world_injection_norm_rad'] < v['window_injection_norm_sum_rad']
                     and abs(v['window_world_attitude_transport_distance']-half) <= 1e-2*half,
                     'collinear signed world injection conclusions')
    return True


def reproduction_differences(fresh, committed, path=''):
    """Paths where a fresh run differs from the committed record.

    Everything must be equal, including exported trace hashes, except that
    80-digit metric strings may differ below 1e-20 relative (mpmath rounding).
    """
    if isinstance(fresh, dict) and isinstance(committed, dict):
        if set(fresh) != set(committed):
            return [path or '/']
        return [d for k in sorted(fresh) for d in reproduction_differences(fresh[k], committed[k], path+'/'+k)]
    if isinstance(fresh, list) and isinstance(committed, list):
        if len(fresh) != len(committed):
            return [path]
        return [d for i, (a, b) in enumerate(zip(fresh, committed))
                for d in reproduction_differences(a, b, f'{path}/{i}')]
    if path.rsplit('/', 1)[-1] in METRICS and isinstance(fresh, str) and isinstance(committed, str):
        with localcontext() as context:
            context.prec = 60
            a, b = Decimal(fresh), Decimal(committed)
            return [] if abs(a-b) <= Decimal('1e-20')*max(1, abs(a), abs(b)) else [path]
    return [] if type(fresh) is type(committed) and fresh == committed else [path]


def run(eigen):
    source = (REPO/HEADER).read_text()
    observed_source, driver = observer_source(source), driver_source()
    cases = []
    with tempfile.TemporaryDirectory(prefix='ou3-world-source-') as directory:
        tmp = Path(directory)
        include = tmp/'kalman_ou_iii'
        include.mkdir()
        (include/HEADER.name).write_text(observed_source)
        (tmp/'driver.cpp').write_text(driver)
        binaries = []
        for name, inc in [('observed', ['-I'+str(tmp)]), ('control', [])]:
            binary = tmp/name
            subprocess.run(['g++', '-O2', '-std=c++20', *inc, '-I'+str(REPO/'src'),
                            '-isystem', str(eigen), str(tmp/'driver.cpp'), '-o', str(binary)], check=True)
            binaries.append(binary)
        for profile in PROFILES:
            observed, control = [json.loads(subprocess.check_output([str(p), profile], text=True))
                                 for p in binaries]
            if any(control[k] != observed[k] for k in control if k != 'events'):
                raise ArithmeticError('observer changed literal source output')
            cases.append({'input_profile': profile, 'literal_terminal_parity': True,
                          'live_step': observed['live_step'], 'refined_step': observed['refined_step'],
                          'active_step': observed['active_step'],
                          'reference': [row[0] for row in observed['reference']],
                          'exported_trace_sha256': hashlib.sha256(
                              json.dumps(observed, sort_keys=True).encode()).hexdigest(),
                          **analyze(observed, profile)})
    return {'qualification': 'OU3_CARRIED_WORLD_FRAME_DIAGNOSTIC_V1',
            'shipping_header_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'observer_header_sha256': hashlib.sha256(observed_source.encode()).hexdigest(),
            'driver_sha256': hashlib.sha256(driver.encode()).hexdigest(),
            'decimal_digits': 80, 'cases': cases, 'source_uniform_verified': False,
            'all_time_magnetic_service_certified': False, 'theorem_closed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expect', type=Path,
                        help='committed record the fresh run must reproduce')
    args = parser.parse_args()
    result = run(args.eigen)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    verify_diagnostic(result)
    if args.expect is not None:
        differences = reproduction_differences(result, json.loads(args.expect.read_text()))
        if differences:
            raise SystemExit('world-frame record does not reproduce at: '+', '.join(differences))
