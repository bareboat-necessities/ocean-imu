"""Signed world-injection transport and the rotating-frame attitude floor.

Proof: docs/ou3-world-frame-rows.md, sections 10--11. Exact rational checks on
supplied words (rotation angles at 80 digits). These lemmas enter
V_next<=rho V+c_d|d|^2 through the aggregate six-column floor s (injection
modulus of the historical rows, hence B_*, J_AG, P_max, rho_0) and, for the
single-injection bound, the |d|^2 term of the reset supply in eta. No
source-uniform transport bound, floor, B_* or rho_0 is promoted.
"""
from fractions import Fraction as F

from .aw_tracking import G, FRACTION, corollary_a_star
from .matrix_certificates import add, identity, is_psd, matmul, transpose
from .world_frame import _cross, _dot, decimal_lower, quaternion_rotation, skew

GYRO_RESIDUAL = F('0.02')
TILT_RETAINED = F('0.10471975511965977')


def rotation_angle(r, dps=80):
    """Angle of a rational rotation matrix (bi-invariant distance to I)."""
    import mpmath as mp
    with mp.workdps(dps):
        tr = sum(mp.mpf(r[i][i].numerator)/r[i][i].denominator for i in range(3))
        return mp.acos(max(-1, min(1, (tr-1)/2)))


def injection_rotation_identity(word, m0):
    """Lemma I*: exact ordered-product identity for the injection rotations.

    word: chronological ('inject', X) with X=Exp(-x) acting on the world
    attitude error M=R_hat' R by M<-X M, and ('rate', P) the exact prediction
    mismatch M<-P M (angle <=h|omega_tilde| up to the branch defect). With J the
    ordered product of the injections alone, M_N M_0'=J K, K the ordered product
    of the rate rotations each conjugated by the preceding injection product.
    By bi-invariance: angle(J)<=angle(M_N)+angle(M_0)+sum angle(P_j).
    No injection norm is summed.
    """
    m, j, k = [row[:] for row in m0], identity(3), identity(3)
    rates = []
    for kind, t in word:
        if matmul(transpose(t), t) != identity(3):
            raise ValueError('rotation matrices required')
        m = matmul(t, m)
        if kind == 'inject':
            j = matmul(t, j)
        elif kind == 'rate':
            k = matmul(matmul(transpose(j), matmul(t, j)), k)
            rates.append(t)
        else:
            raise ValueError('unknown operation')
    assert matmul(m, transpose(m0)) == matmul(j, k)
    lhs = rotation_angle(j)
    rhs = rotation_angle(m)+rotation_angle(m0)+sum(rotation_angle(p) for p in rates)
    assert lhs <= rhs
    return {'injection_angle': lhs, 'bound': rhs}


def third_order_reset_remainder(x, dps=80):
    """N=Exp(-[x])(I+[x]/2)=I-[x]/2+R, R=sum_{n>=3}(1-n/2)(-X)^n/n!.

    |X^n|<=|x|^n for the skew X gives |R|<=|x|^3/6 when |x|<=1: the quadratic
    term cancels, so the literal world reset factor is a half-angle rotation up
    to a second-order stretch and third-order rotation.
    """
    import mpmath as mp
    with mp.workdps(dps):
        v = mp.matrix([mp.mpf(F(c).numerator)/F(c).denominator for c in x])
        th = mp.norm(v)
        X = mp.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        exp_neg = mp.eye(3)-mp.sin(th)/th*X+(1-mp.cos(th))/th**2*X*X
        rem = exp_neg*(mp.eye(3)+X/2)-mp.eye(3)+X/2
        norm = max(mp.svd(rem, compute_uv=False))
        assert th <= 1 and norm <= th**3/6
        return norm, th**3/6


def rotating_frame_floor(weights, forces, rotations, field_unit, anchor):
    """Corollary A**: exact attitude Gram with rotation factors Q_k in the rows.

    Rows [u_k]x Q_k (weights alpha) and one magnetic row [b]x Q_j. Since
    |u x Q theta|=|Q'u x theta|, convexity gives the floor
    |Q_j v_bar x b|^2/(|v_bar|^2+1), v_bar=sum alpha Q_k'u_k, and
    Q_j v_bar=u_bar+sum alpha(Q_j Q_k'-I)u_k: the loss is a SIGNED weighted mean
    of relative rotations, not max|Q_k-Q_j|.
    """
    gram = matmul(transpose(matmul(skew(field_unit), rotations[anchor])),
                  matmul(skew(field_unit), rotations[anchor]))
    vbar = [F(0)]*3
    for w, u, q in zip(weights, forces, rotations):
        row = matmul(skew(u), q)
        gram = add(gram, matmul(transpose(row), row), F(w))
        qu = [r[0] for r in matmul(transpose(q), [[c] for c in u])]
        vbar = [a+F(w)*c for a, c in zip(vbar, qu)]
    qv = [r[0] for r in matmul(rotations[anchor], [[c] for c in vbar])]
    cross = _cross(qv, field_unit)
    floor = _dot(cross, cross)/(_dot(vbar, vbar)+1)
    assert is_psd(add(gram, identity(3), -floor))
    return floor


def rotating_frame_example():
    """Opposite relative rotations cancel in the signed mean; one-sided do not."""
    b = [F(7, 25), F(0), F(24, 25)]
    q_plus = quaternion_rotation((40, 0, 1, 0))     # about e_y (the b, e_z plane)
    q_minus = transpose(q_plus)
    forces = [[F(0), F(0), F(-1)]]*4
    weights = [F(1, 4)]*4
    signed = rotating_frame_floor(weights, forces, [q_minus, identity(3), identity(3), q_plus], b, 1)
    one_sided = min(rotating_frame_floor(weights, forces, [q, identity(3), identity(3), q], b, 1)
                    for q in (q_plus, q_minus))
    none = corollary_a_star(0, 0, F(7, 25))['normalized_gram_floor']
    assert one_sided < signed <= none
    return {'relative_rotation_rad': '2 atan(1/40)', 'signed_floor': str(decimal_lower(signed)),
            'one_sided_floor': str(decimal_lower(one_sided)), 'no_rotation_floor': str(none)}


def relaxed_reset_example_exclusion():
    """DEAD_END 17 needs single injections |d|=4, 4, 8/3 rad.

    One reset with no intervening prediction has angle(Exp(-x))<=alpha_-+alpha_+
    (Lemma I* with an empty rate product), so each such injection needs a world
    attitude error of at least min(4,2 pi-4)/2 rad before or after it.
    """
    import mpmath as mp
    with mp.workdps(30):
        need = min(mp.mpf(4), 2*mp.pi-4)/2
    return {'required_attitude_error_rad_lower': mp.nstr(need, 12),
            'retained_tilt_domain_rad': str(TILT_RETAINED)}


def feasibility():
    """Worst admissible ratios for two injection formulations on 16-s windows.

    Perturbative: |A~_k-A~_j|<=delta needs delta^2(u_rms^2+1)<gamma*; to first
    order A~ rotates by half the injection rotation, which Lemma I* bounds by
    2 alpha+int|omega_tilde| over the half window (8 s).
    Rotating frame (Corollary A**): the loss is |sum alpha Delta_k e_z| plus
    max|Delta| mean|a_hat|/g against sigma_w; an adversarial residual that
    flips sign at the anchor makes every partial transport one-signed, so the
    signed mean over the half window is T/4 times the rate.
    """
    rows = []
    gamma = corollary_a_star(0, 0)['normalized_gram_floor']
    for label, alpha, rate in (('retained tilt, bias exact', TILT_RETAINED, GYRO_RESIDUAL),
                               ('vanishing radius, bias exact', F(0), GYRO_RESIDUAL),
                               ('vanishing radius, bias error B_g', F(0), 2*GYRO_RESIDUAL),
                               ('vanishing radius, invariant only', F(0), F('0.54'))):
        for u_rms2, mean_aw in ((F(2), F('0.5')), (F('4.6'), F(3))):
            delta = (2*alpha+8*rate)/2
            pert = delta*delta*(u_rms2+1)/gamma
            rot = ((2*alpha+4*rate)/2+delta*mean_aw/G)/FRACTION
            rows.append({'case': label, 'u_rms_squared': str(u_rms2),
                         'mean_abs_aw_mps2': str(mean_aw),
                         'perturbative_ratio_squared': str(decimal_lower(pert)),
                         'rotating_frame_ratio': str(decimal_lower(rot)),
                         'perturbative_feasible': pert < 1, 'rotating_frame_feasible': rot < 1})
    return rows



def literal_g0_transport_feasibility(*, injection_free_floor,
                                     signed_rotation_charge,
                                     cubic_reset_charge=F(0)):
    """Kill test for transporting an injection-free G0 singular floor.

    Weyl gives sigma_min(O_literal)>=sigma_min(O_0)-||Delta O||.  This helper
    records the only margin a perturbation-style literal extension would have.
    It does not manufacture ||Delta O|| from injection norm sums: callers must
    supply a signed/half-angle source bound.  A nonpositive margin kills this
    formulation and requires a direct rotating-frame Gram proof.
    """
    floor,charge,cubic=map(F,(injection_free_floor,signed_rotation_charge,cubic_reset_charge))
    if min(floor,charge,cubic)<0:
        raise ValueError('nonnegative transport quantities required')
    margin=floor-charge-cubic
    return {'injection_free_singular_floor':str(floor),
            'signed_rotation_charge':str(charge),
            'cubic_reset_charge':str(cubic),
            'literal_singular_floor_margin':str(margin),
            'perturbative_literal_G0_feasible':margin>0,
            'requires_direct_rotating_frame_if_failed':margin<=0,
            'norm_summed_injection_budget_allowed':False}


def certificate():
    m0 = quaternion_rotation((30, 1, 2, -1))
    x1, x2 = quaternion_rotation((25, -1, 0, 2)), quaternion_rotation((50, 1, -1, 0))
    r1, r2, r3 = (quaternion_rotation(q) for q in ((200, 0, 1, 1), (150, 1, 0, -1), (300, 1, 1, 0)))
    lemma = injection_rotation_identity(
        [('rate', r1), ('inject', x1), ('rate', r2), ('inject', x2), ('rate', r3)], m0)
    single = injection_rotation_identity([('inject', x1)], m0)
    rem, bound = third_order_reset_remainder([F(3, 10), F(-1, 5), F(1, 2)])
    small, small_bound = third_order_reset_remainder([F(1, 100), F(0), F(-1, 200)])
    import mpmath as mp
    fmt = lambda v: mp.nstr(v, 20)
    table = feasibility()
    return {
        'qualification': 'OU3_SIGNED_INJECTION_TRANSPORT_V1',
        'scope': 'real-arithmetic regular A21 world frame; supplied exact words, angles at 80 digits',
        'lemma_I_star': {
            'ordered_product_identity_exact': True,
            'supplied_injection_angle': fmt(lemma['injection_angle']),
            'supplied_attitude_and_rate_bound': fmt(lemma['bound']),
            'single_supplied_injection_angle': fmt(single['injection_angle']),
            'single_injection_bound': fmt(single['bound']),
            'norm_sum_of_injections_used': False,
        },
        'reset_factor_third_order': {
            'quadratic_term_cancels': True,
            'supplied_remainder_norm': fmt(rem), 'supplied_cubic_bound': fmt(bound),
            'small_remainder_norm': fmt(small), 'small_cubic_bound': fmt(small_bound),
        },
        'half_angle_product_remainder_source_bound': False,
        'rotating_frame_attitude_floor': rotating_frame_example(),
        'relaxed_reset_cancellation_exclusion': relaxed_reset_example_exclusion(),
        'injection_formulation_feasibility': table,
        'perturbative_charge_feasible_somewhere': any(r['perturbative_feasible'] for r in table),
        'rotating_frame_feasible_at_vanishing_radius_bias_exact': all(
            r['rotating_frame_feasible'] for r in table if r['case'] == 'vanishing radius, bias exact'),
        'source_uniform_signed_injection_bound': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
