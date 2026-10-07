"""Literal CoG measurement feedback in a moving nominal-world frame.

All 21 states and the same P/K/Joseph operands remain. This is an exact change
of frame, including its derivative, not an invariant-EKF replacement. Post-
injection frame transport, reference/noise changes and branch jumps remain.
"""
from __future__ import annotations
from fractions import Fraction as F

from .information_shear_word import zeros, trace
from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, transpose
from .planar_linked_riccati_mean import product
from .world_frame import skew


BODY_BLOCKS = (0, 3, 18)  # attitude, BG, BA; v,p,S,AW stay world coordinates


def frame(rotation):
    """T maps nominal-world tangent components to shipping body components."""
    if len(rotation) != 3 or any(len(row) != 3 for row in rotation):
        raise ValueError('3 by 3 rotation required')
    r = rotation
    determinant = sum((r[0][j]*(r[1][(j+1)%3]*r[2][(j+2)%3]
                               -r[1][(j+2)%3]*r[2][(j+1)%3]) for j in range(3)), F(0))
    if product(transpose(r), r) != identity(3) or determinant != 1:
        raise ValueError('exact SO(3) rotation required; arithmetic defects are separate')
    t = identity(21)
    for offset in BODY_BLOCKS:
        for i in range(3):
            t[offset+i][offset:offset+3] = r[i][:]
    return t


def connection(omega):
    omega = skew(omega)
    out = zeros(21, 21)
    for offset in BODY_BLOCKS:
        for i in range(3):
            out[offset+i][offset:offset+3] = omega[i][:]
    return omega, out


def world_rows(force_world, reference):
    """force_world is actual nominal aw-g, not physical acceleration."""
    acc, mag = zeros(3, 21), zeros(3, 21)
    sf, sb = skew(force_world), skew(reference)
    for i in range(3):
        acc[i][:3], mag[i][:3] = [-x for x in sf[i]], [-x for x in sb[i]]
        acc[i][15+i] = acc[i][18+i] = F(1)
    return acc, mag


def row_differentials(daw, dreference):
    acc, mag = zeros(3, 21), zeros(3, 21)
    for rows, variation in ((acc, daw), (mag, dreference)):
        for i, row in enumerate(skew(variation)):
            rows[i][:3] = [-x for x in row]
    return acc, mag


def pullback_operands(rotation, P, H, Rnoise, residual):
    t, rt = frame(rotation), transpose(rotation)
    return {'T': t, 'P': product(transpose(t), P, t),
            'H': product(rt, H, t), 'R': product(rt, Rnoise, rotation),
            'r': product(rt, residual)}


def pullback_differentials(rotation, omega, P, H, Rnoise, residual, dP, dH, dR, dr):
    """omega is nominal R' dR, NOT the physical comparison-angle variation."""
    b = pullback_operands(rotation, P, H, Rnoise, residual)
    w, o = connection(omega)
    t, rt = b['T'], transpose(rotation)
    return {
        'dP': add(add(product(transpose(t), dP, t), product(o, b['P']), -1), product(b['P'], o)),
        'dH': add(add(product(rt, dH, t), product(w, b['H']), -1), product(b['H'], o)),
        'dR': add(add(product(rt, dR, rotation), product(w, b['R']), -1), product(b['R'], w)),
        'dr': add(product(rt, dr), product(w, b['r']), -1),
    }


def correction(P, H, Rnoise, residual, held_ba=False):
    """Literal PCt/Joseph for default CoG acc/mag/S branches.

    H is the full innovation row. Held BA omits BA gain columns AND rows;
    this is not replaced by an unmasked Riccati covariance identity.
    """
    ldlt(P)
    ldlt(Rnoise)
    Hgain = [row[:] for row in H]
    if held_ba:
        for row in Hgain:
            row[18:21] = [F(0)]*3
    S = add(product(H, P, transpose(H)), Rnoise)
    B = product(P, transpose(Hgain))
    if held_ba:
        B[18:21] = zeros(3, len(H))
    K = product(B, inverse(S))
    C = add(add(add(P, product(K, transpose(B)), -1), product(B, transpose(K)), -1),
            product(K, S, transpose(K)))
    return {'S': S, 'B': B, 'K': K, 'C': C, 'increment': product(K, residual)}


def covariance_row_charge(P, H, Rnoise, dH):
    """Exact linked Fisher charge, with no independent covariance/gain maxima.

    Covered unmasked optimal correction only. Same-history noise/reference
    ports add separately; held BA uses its reached active reduction.
    """
    ldlt(P)
    ldlt(Rnoise)
    S = add(product(H, P, transpose(H)), Rnoise)
    K = product(P, transpose(H), inverse(S))
    C = add(P, product(K, H, P), -1)
    Jc = inverse(C)
    W = add(inverse(Rnoise), inverse(S), -1)
    assert product(transpose(K), Jc, K) == W
    q = add(product(K, dH, C), product(C, transpose(dH), transpose(K)))
    charge = trace(product(Jc, q, Jc, q))
    upper = 4*trace(product(C, transpose(dH), W, dH))
    assert 0 <= charge <= upper
    return {'row_port': [[-x for x in row] for row in q],
            'Fisher_charge': charge, 'linked_upper': upper,
            'innovation_information_increment': W}


def magnetic_connection_charge(P, H, Rnoise, residual, e, omega):
    """Fixed world reference, isotropic fixed R: exact remaining shear port.

    dr=-H de in shipping coordinates. Frame dependence contributes
    -K([omega] r+H Omega e); omega is the NOMINAL frame variation.
    """
    scale = Rnoise[0][0]
    if Rnoise != [[scale*F(i == j) for j in range(3)] for i in range(3)]:
        raise ValueError('isotropic noise needed for zero rotational noise port')
    w, o = connection(omega)
    S = add(product(H, P, transpose(H)), Rnoise)
    K = product(P, transpose(H), inverse(S))
    C = add(P, product(K, H, P), -1)
    v = add(product(w, residual), product(H, o, e))
    f = [[-x for x in row] for row in product(K, v)]
    W = add(inverse(Rnoise), inverse(S), -1)
    actual = product(transpose(f), inverse(C), f)[0][0]
    linked = product(transpose(v), W, v)[0][0]
    assert actual == linked
    return {'port': f, 'measurement_port': v, 'energy': actual,
            'linked_energy': linked, 'information_increment': W}


def transported_factor(before_rotation, after_rotation, factor):
    """Actual prediction/reset factor in endpoint frames, no isometry assumption."""
    return product(transpose(frame(after_rotation)), factor, frame(before_rotation))


def aw_connection(daw):
    """Nilpotent AW/attitude connection in the nominal-world frame."""
    out = zeros(21, 21)
    for i, row in enumerate(skew(daw)):
        out[15+i][:3] = [-x for x in row]
    return out


def aw_shear(aw):
    """L=I+N(aw), L^-1=I-N(aw), N(a)N(b)=0 for all a,b.

    H_acc L^-1=[+[g]x,0,0,0,0,I,I]. No nominal AW cap is
    needed for invertibility. Endpoint/frame derivatives are NOT discarded.
    """
    return add(identity(21), aw_connection(aw))


def aw_shear_differentials(aw, daw, P, dP, e, de):
    """Complete moving-frame differential; not a fixed congruence tangent."""
    L, Gamma = aw_shear(aw), aw_connection(daw)
    p = product(L, P, transpose(L))
    ep = product(L, e)
    dp = add(add(product(L, dP, transpose(L)), product(Gamma, p)), product(p, transpose(Gamma)))
    dep = add(product(L, de), product(Gamma, ep))
    eta = add(dep, product(dp, inverse(p), ep), -1)
    original = add(de, product(dP, inverse(P), e), -1)
    linked = add(product(L, original), product(p, transpose(Gamma), inverse(p), ep), -1)
    assert eta == linked
    return {'L': L, 'connection': Gamma, 'P': p, 'dP': dp, 'e': ep,
            'de': dep, 'eta': eta}


def row_connection_coboundary(P, H, Rnoise, Gamma):
    """dH=H Gamma: exact covariance row port is a metric coboundary.

    Valid optimal branch only. Reached held BA uses its active effective-noise
    reduction. dR and post-correction connection changes remain separate.
    """
    ldlt(P)
    ldlt(Rnoise)
    S = add(product(H, P, transpose(H)), Rnoise)
    K = product(P, transpose(H), inverse(S))
    A = add(identity(len(P)), product(K, H), -1)
    C = product(A, P)
    dH = product(H, Gamma)
    U = [[-x for x in row] for row in add(product(K, dH, C), product(C, transpose(dH), transpose(K)))]
    root = add(product(Gamma, P), product(P, transpose(Gamma)))
    terminal = add(product(Gamma, C), product(C, transpose(Gamma)))
    coboundary = add(product(A, root, transpose(A)), terminal, -1)
    assert U == coboundary
    return {'row_port': U, 'root_connection': root, 'terminal_connection': terminal}


def fixed_row_joint_balance(P, H, Rnoise, eta, dP, mismatch, weight=F(1)):
    """Exact pre-injection joint loss after ALL frame terms enter mismatch.

    eta+=A eta+K mismatch, dP+=A dP A'. Noise/reference variations
    and the next moving-frame jump are not covered by this substep alone.
    The signed work remains endogenous; mismatch is not independent noise.
    """
    ldlt(P)
    ldlt(Rnoise)
    if weight <= 0:
        raise ValueError('positive covariance weight required')
    S = add(product(H, P, transpose(H)), Rnoise)
    K = product(P, transpose(H), inverse(S))
    A = add(identity(len(P)), product(K, H), -1)
    C, dp = product(A, P), product(A, dP, transpose(A))
    ep = add(product(A, eta), product(K, mismatch))
    J, Jc = inverse(P), inverse(C)
    initial = product(transpose(eta), J, eta)[0][0]+weight*trace(product(J, dP, J, dP))
    terminal = product(transpose(ep), Jc, ep)[0][0]+weight*trace(product(Jc, dp, Jc, dp))
    lossP = trace(product(J, dP, J, dP))-trace(product(Jc, dp, Jc, dp))
    corrected = add(product(H, eta), mismatch, -1)
    action = product(transpose(corrected), inverse(S), corrected)[0][0]+weight*lossP
    supply = product(transpose(mismatch), inverse(Rnoise), mismatch)[0][0]
    assert terminal-initial == -action+supply
    return {'initial': initial, 'terminal': terminal, 'corrected_action': action,
            'linked_mismatch_supply': supply, 'covariance_loss': lossP}


def planar_acc_mismatch_charge(P, e, daw, omega, divided_force, source, Rnoise):
    """Same-history physical substitution, not an independent operand box.

    On the fixed-input central planar comparison, e_theta=t e_y and
    divided_force=(R_y(-t)-I)(a_phys-g)/t, continuously extended at t=0.
    source retains delivered FAST/model forcing; it is not an invented noise
    channel. The helper checks algebra, not the physical qualification of inputs.
    """
    ldlt(P)
    ldlt(Rnoise)
    if e[0][0] or e[2][0] or daw[1]:
        raise ValueError('central planar pitch/AW comparison required')
    Jy = skew([0, 1, 0])
    B = zeros(3, 21)
    for i in range(3):
        B[i][1] = daw[i]-omega*divided_force[i]
        B[i][15+i] = omega
    m0 = product(Jy, B, e)
    q = [[-omega*x for x in row] for row in product(Jy, [[x] for x in source])]
    m = add(m0, q)
    precision = inverse(Rnoise)
    E = product(transpose(e), inverse(P), e)[0][0]
    linked_form = trace(product(precision, Jy, B, P, transpose(B), transpose(Jy)))
    actual = product(transpose(m0), precision, m0)[0][0]
    cross = 2*product(transpose(m0), precision, q)[0][0]
    square = product(transpose(q), precision, q)[0][0]
    total = product(transpose(m), precision, m)[0][0]
    assert 0 <= actual <= E*linked_form and total == actual+cross+square
    return {'physical_curvature_row': B, 'mismatch': m, 'comparison_energy': E,
            'linked_covariance_form': linked_form, 'curvature_energy': actual,
            'curvature_upper': E*linked_form, 'source_cross': cross,
            'source_square': square, 'full_mismatch_energy': total}


def moving_frame_word_ports(B, P, C, e, ep, d, U, v, L, Lnext, connection0, connection1):
    """Exact port conjugacy, including the generated-covariance gain score.

    connection_i=L_i^-1 dL_i. All operands/ports are linked to ONE history.
    This helper changes coordinates; it does not certify a shipping domain.
    No invertibility of B or optimal-correction hypothesis is needed.
    """
    Li = inverse(L)
    Z0 = add(product(connection0, P), product(P, transpose(connection0)))
    Z1 = add(product(connection1, C), product(C, transpose(connection1)))
    Uinside = add(add(U, Z1), product(B, Z0, transpose(B)), -1)
    vinside = add(add(v, product(connection1, ep)), product(B, connection0, e), -1)
    vinside = add(vinside, product(B, Z0, d), -1)
    return {'base': product(Lnext, B, Li),
            'gain_score': product(transpose(Li), d),
            'covariance_port': product(Lnext, Uinside, transpose(Lnext)),
            'mean_port': product(Lnext, vinside)}


def aw_frame_boundary_work(P, e, eta, dP, daw, weight=F(1)):
    """Exact three-dimensional endpoint quadratic, NOT an absorption claim.

    W_tilde-W=2 h' daw+daw' Q daw. Internal connection work telescopes
    only together with transformed action/ports, never on its own.
    Q uses the actual conditional AW precision J_aw,aw, not P_aw,aw^-1.
    """
    ldlt(P)
    if len(P) != 21 or weight <= 0 or dP != transpose(dP):
        raise ValueError('complete symmetric tangent and positive weight required')
    J, Gamma = inverse(P), aw_connection(daw)
    Je = product(J, e)
    Cg = add(product(Gamma, P), product(P, transpose(Gamma)))
    shifted_eta = add(eta, product(P, transpose(Gamma), Je), -1)
    shifted_dP = add(dP, Cg)
    old = product(transpose(eta), J, eta)[0][0]+weight*trace(product(J, dP, J, dP))
    new = product(transpose(shifted_eta), J, shifted_eta)[0][0]
    new += weight*trace(product(J, shifted_dP, J, shifted_dP))
    generators = [aw_connection([F(i == j) for i in range(3)]) for j in range(3)]
    h, Q = zeros(3, 1), zeros(3, 3)
    for j, Gj in enumerate(generators):
        h[j][0] = -product(transpose(Je), Gj, eta)[0][0]+2*weight*trace(product(dP, J, Gj))
        for k, Gk in enumerate(generators):
            Q[j][k] = product(transpose(Je), Gj, P, transpose(Gk), Je)[0][0]
            Q[j][k] += 2*weight*trace(product(J, Gj, P, transpose(Gk)))
    ldlt(Q)
    a = [[x] for x in daw]
    work = 2*product(transpose(h), a)[0][0]+product(transpose(a), Q, a)[0][0]
    minimizer = [[-x for x in row] for row in product(inverse(Q), h)]
    lower = -product(transpose(h), inverse(Q), h)[0][0]
    completed = product(transpose(add(a, minimizer, -1)), Q, add(a, minimizer, -1))[0][0]+lower
    assert new-old == work == completed and work >= lower >= -old
    return {'linear_coefficient': h, 'quadratic_coefficient': Q,
            'signed_boundary_work': work, 'completed_square_lower': lower,
            'unconstrained_minimizer': minimizer,
            'initial_storage': old, 'transformed_storage': new}


def aw_conditional_sync(P, e, increment, daw):
    """Actual PSD AW increment improves conditional precision, not a new floor.

    The caller must supply the shipping same-history increment. This algebra
    does not qualify a face/target, branch derivative or inherited upper bound.
    The connection-square decrease is bookkeeping, not an additional W loss.
    """
    ldlt(P)
    if len(P) != 21 or not is_psd(increment) or len(increment) != 3:
        raise ValueError('full covariance and actual 3 by 3 PSD AW increment required')
    a, o = list(range(15, 18)), list(range(15))+list(range(18, 21))
    def block(rows, cols):
        return [[P[i][j] for j in cols] for i in rows]
    regression = product(block(a, o), inverse(block(o, o)))
    conditional = add(block(a, a), product(regression, block(o, a)), -1)
    Jaa, Jnext = inverse(conditional), inverse(add(conditional, increment))
    precision_loss = add(Jaa, Jnext, -1)
    assert is_psd(precision_loss)
    residual = add([e[i] for i in a], product(regression, [e[i] for i in o]), -1)
    V = product(transpose(residual), Jaa, residual)[0][0]
    Vnext = product(transpose(residual), Jnext, residual)[0][0]
    sa = skew(daw)
    connection_loss = trace(product(precision_loss, sa, block(range(3), range(3)), transpose(sa)))
    assert Vnext <= V and connection_loss >= 0
    return {'conditional_covariance': conditional, 'conditional_precision': Jaa,
            'next_conditional_precision': Jnext, 'conditional_comparison_energy': V,
            'next_conditional_comparison_energy': Vnext,
            'Fisher_connection_square_decrease': connection_loss}


def conditional_aw_storage(P, e, eta, dP, weight=F(1)):
    """Exact conditional/complementary joint storage, all cross blocks retained.

    A Gaussian block identity for the actual covariance, not a replacement
    observer or a new physical compatibility quotient.
    """
    ldlt(P)
    if len(P) != 21 or dP != transpose(dP) or weight <= 0:
        raise ValueError('full symmetric covariance tangent and positive weight required')
    a, o = list(range(15, 18)), list(range(15))+list(range(18, 21))
    def block(M, rows, cols):
        return [[M[i][j] for j in cols] for i in rows]
    B, dB = block(P, o, o), block(dP, o, o)
    Binv = inverse(B)
    T = product(block(P, a, o), Binv)
    C = add(block(P, a, a), product(T, block(P, o, a)), -1)
    dT = product(add(block(dP, a, o), product(T, dB), -1), Binv)
    dC = add(add(add(block(dP, a, a), product(block(dP, a, o), transpose(T)), -1),
                 product(T, block(dP, o, a)), -1), product(T, dB, transpose(T)))
    J, Cinv = inverse(P), inverse(C)
    j, s = product(J, e), product(J, eta)
    sa, so = [s[i] for i in a], [s[i] for i in o]
    ell = add(so, product(transpose(T), sa))
    conditional = product(transpose(sa), C, sa)[0][0]+weight*trace(product(Cinv, dC, Cinv, dC))
    complementary = product(transpose(ell), B, ell)[0][0]+weight*trace(product(Binv, dB, Binv, dB))
    complementary += 2*weight*trace(product(Cinv, dT, B, transpose(dT)))
    total = product(transpose(eta), J, eta)[0][0]+weight*trace(product(J, dP, J, dP))
    assert total == conditional+complementary and min(conditional, complementary) >= 0
    return {'regression': T, 'conditional_covariance': C, 'd_regression': dT,
            'd_conditional_covariance': dC, 'conditional_information': [j[i] for i in a],
            'conditional_score': sa, 'conditional_storage': conditional,
            'complementary_storage': complementary, 'total_storage': total}


def aw_precision_process_balance(P, Fmap, process, phi, daw):
    """Same-operation conditional precision / OU-S-chain action balance.

    Nominal-root default one-way tuner scope gives daw+=phi*daw. Source
    variations of phi and both physical/chart discrepancies remain ports.
    The signed integrated-chain and AG terms are not replaced by maxima.
    """
    ldlt(P)
    if len(P) != 21 or not is_psd(process):
        raise ValueError('complete covariance and actual PSD process required')
    a, o = list(range(15, 18)), list(range(15))+list(range(18, 21))
    for i in range(3):
        if Fmap[15+i] != [phi*F(j == 15+i) for j in range(21)]:
            raise ValueError('literal shared scalar OU AW row required')
    C = add(product(Fmap, P, transpose(Fmap)), process)
    J, Jnext = inverse(P), inverse(C)
    loss = add(J, product(transpose(Fmap), Jnext, Fmap), -1)
    assert is_psd(loss)
    def block(M, rows, cols):
        return [[M[i][j] for j in cols] for i in rows]
    Jaa, Jaa_next, Qa = block(J, a, a), block(Jnext, a, a), block(loss, a, a)
    B = block(Fmap, o, a)
    cross = product(block(Jnext, a, o), B)
    chain = add([[phi*x for x in row] for row in add(cross, transpose(cross))],
                product(transpose(B), block(Jnext, o, o), B))
    assert [[phi*phi*x for x in row] for row in Jaa_next] == add(add(Jaa, Qa, -1), chain, -1)
    S = skew(daw)
    old_theta, new_theta = block(P, range(3), range(3)), block(C, range(3), range(3))
    old = trace(product(Jaa, S, old_theta, transpose(S)))
    new = phi*phi*trace(product(Jaa_next, S, new_theta, transpose(S)))
    positive_loss = trace(product(Qa, S, new_theta, transpose(S)))
    chain_work = -trace(product(chain, S, new_theta, transpose(S)))
    ag_work = trace(product(Jaa, S, add(new_theta, old_theta, -1), transpose(S)))
    assert positive_loss >= 0 and new-old == -positive_loss+chain_work+ag_work
    return {'conditional_precision': Jaa, 'next_conditional_precision': Jaa_next,
            'AW_process_action': Qa, 'integrated_chain_precision_term': chain,
            'connection_before': old, 'connection_after': new,
            'positive_connection_loss': positive_loss,
            'signed_chain_work': chain_work, 'signed_AG_work': ag_work}


def qualified_aw_precision_ceiling():
    """Corollary of the existing scaled LIN path certificate, not new sampling."""
    from .lin_matrix_certificate import action_matrix, SCALES
    value = action_matrix()[3][3]/F(SCALES[3])**2
    upper = F(15942618)
    assert 0 < value < upper
    return {'exact_upper': str(value), 'integer_upper': str(upper),
            'scope': 'same qualified real regular A21 16-second LIN-path profile and activation; no automatic H18/profile/float transfer',
            'conditional_not_marginal_precision': True,
            'sufficient_for_endpoint_Schur_margin': False}


def conditional_aw_process_coercivity():
    """Rational positive loss on the conditional AW zero-action block only.

    Uses the literal integrated OU column/noise, its proved polynomial defect
    and existing isotropic regular-profile AW ceiling. This is not a bound on
    remaining signed work or a complete-word Schur margin. See app:aw-conditional-loss.
    """
    from .lin_path_certificate import small_x_source_defect
    from .aw_covariance_ceiling import ceiling
    eps, _, gram0, _ = small_x_source_defect()
    inverse_norm = max(sum(abs(v) for v in row) for row in inverse(gram0))
    assert inverse_norm == 81060
    # sigma^2*x >= (.05)^2*(.004/12), x<=.3, and ||D_h^-1 F e_a||^2<=4.
    shorted_noise = (1-eps)*F(1, 1200000)/(4*F(13, 10)**2*inverse_norm)
    c = shorted_noise/(ceiling(eps)+shorted_noise)
    assert F(1, 10**14) < c < 1
    return {'shorted_process_noise_lower': str(shorted_noise),
            'conditional_covariance_upper': str(ceiling(eps)),
            'mean_loss_fraction_lower': str(c),
            'covariance_loss_fraction_lower': str(2*c-c*c),
            'simple_loss_fraction_lower': '1/100000000000000',
            'scope': 'base process action on conditional AW mean/covariance directions; existing real isotropic regular profile, h in [.004,.006], tau in [.02,12], sigma in [.05,4]',
            'mixed_block_or_signed_work_absorption': False,
            'uniform_complete_word_gap': False}


def planar_pitch_prediction_calculus(step, angle, gyro_density, bias_density):
    """Literal real-operation pitch/BG process, not a uniform word-gap test.

    On the fixed-input planar stratum the polynomial rotation/integral terms
    annihilate e_y. F_E,Q_E have zero MEKF-root derivative. The NORMALIZED
    small-angle quaternion is not an exact exponential: its mean derivative
    leaves the rank-one ds port derived in app:planar-even-process-work.
    Rounded-coefficient/FMA/libm defects remain separate.
    """
    h, x, qg, qb = map(F, (step, angle, gyro_density, bias_density))
    if h <= 0 or abs(x) >= F(1, 100) or min(qg, qb) <= 0:
        raise ValueError('positive process densities and literal small-angle branch required')
    w = 1-x*x/8+x**4/384
    v = x/2-x**3/48+x**5/3840
    dw = -x/4+x**3/96
    dv = F(1, 2)-x*x/16+x**4/768
    norm2 = w*w+v*v
    angle_derivative = 2*(w*dv-v*dw)/norm2
    defect = x**6*(1920-80*x*x+x**4)/(14745600*norm2)
    if angle_derivative != 1-defect or not 0 <= defect <= x**6/7680:
        raise ArithmeticError('literal normalized quaternion derivative bound failed')
    Q = [[qg*h+qb*h**3/3, qb*h*h/2], [qb*h*h/2, qb*h]]
    ldlt(Q)
    source_reader = h*h*defect*defect*inverse(Q)[0][0]
    source_upper = h*defect*defect/qg
    if source_reader > source_upper:
        raise ArithmeticError('linked pitch process charge bound failed')
    return {'pitch_BG_transition': [[F(1), h], [F(0), F(1)]],
            'pitch_BG_process': Q, 'quaternion_norm_squared': norm2,
            'literal_angle_derivative': angle_derivative,
            'mean_covariance_transition_defect': -h*defect,
            'auxiliary_charge_per_squared_BG_variation': source_reader,
            'auxiliary_charge_upper': source_upper,
            'even_process_covariance_port_zero_fixed_input': True,
            'odd_process_port_zero': False,
            'uniform_complete_gap_verified': False}


def scalar_aw_innovation_reader(P, D, h, noise, aw_index, dh=None, dnoise=0):
    """Actual scalar correction's AA decrement and SAME Fisher reader (CR6).

    Fixed h/R covariance part satisfies dI_cov^2 <= I*(2 P_aa-I)*L_P.
    Literal row/noise variations are retained as one signed coefficient port;
    they are not included in that bound or classified as external by this tool.
    Scalar odd acc/mag/S channels use the same operation's P,h,R and tangent.
    """
    ldlt(P)
    n, a, noise = len(P), aw_index, F(noise)
    if (not 0 <= a < n or D != transpose(D) or noise <= 0
            or len(h) != 1 or len(h[0]) != n):
        raise ValueError('regular scalar row/noise and symmetric same-operation tangent required')
    dh = zeros(1, n) if dh is None else dh
    J = inverse(P)
    s = product(h, P, transpose(h))[0][0]+noise
    b = product(P, transpose(h))[a][0]
    db_cov = product(D, transpose(h))[a][0]
    ds_cov = product(h, D, transpose(h))[0][0]
    db_port = product(P, transpose(dh))[a][0]
    ds_port = 2*product(dh, P, transpose(h))[0][0]+F(dnoise)
    I = b*b/s
    di_cov = 2*b*db_cov/s-b*b*ds_cov/(s*s)
    di_port = 2*b*db_port/s-b*b*ds_port/(s*s)
    loss = (2*product(h, D, J, D, transpose(h))[0][0]/s
            - ds_cov*ds_cov/(s*s))
    K = add([[2*x for x in row] for row in J],
            [[x/s for x in row] for row in product(transpose(h), h)], -1)
    reader = [[2*F(i == a)-b*h[0][i]/s for i in range(n)]]
    reader_norm = product(reader, inverse(K), transpose(reader))[0][0]
    assert reader_norm == 2*P[a][a]-I
    assert di_cov*di_cov <= I*reader_norm*loss
    return {'innovation': s, 'AA_information_decrement': I,
            'AA_decrement_covariance_tangent': di_cov,
            'AA_decrement_coefficient_port': di_port,
            'AA_decrement_full_tangent': di_cov+di_port,
            'same_operation_Fisher_loss': loss,
            'sharp_reader_norm_squared': reader_norm,
            'sharp_charge_coefficient': I*reader_norm,
            'coefficient_port_absorption_verified': False}


def scalar_aw_face_fisher_balance(P, D, target, aw_index, dtarget=0,
                                 zero_gap_active=False):
    """Exact active scalar face in the EXISTING Fisher metric (OF2).

    Shipping use: AW_y in the planar odd block, fixed inherited target and
    source/private state, regular A21 branch. A supplied target tangent stays
    in the same signed gap; this tool does not classify its causal origin.
    No replacement covariance law,
    independent cross-covariance bound, or reachability assertion is supplied.
    Held effective-noise ports and off-parity tangents need the full causal
    derivative instead of this scalar specialization.
    """
    ldlt(P)
    n, a, target, dtarget = len(P), aw_index, F(target), F(dtarget)
    if (not 0 <= a < n or D != transpose(D) or target < P[a][a]
            or (target == P[a][a] and not zero_gap_active)):
        raise ValueError('symmetric tangent and strictly active scalar AW face required; zero gap needs its directional branch')
    o = [i for i in range(n) if i != a]
    B = [[P[i][j] for j in o] for i in o]
    dB = [[D[i][j] for j in o] for i in o]
    T = product([[P[a][j] for j in o]], inverse(B))
    dT = product(add([[D[a][j] for j in o]], product(T, dB), -1), inverse(B))
    beta = product(T, B, transpose(T))[0][0]
    d_beta = (product(T, dB, transpose(T))[0][0]
              + 2*product(dT, B, transpose(T))[0][0])
    C, Cnext, dC = P[a][a]-beta, target-beta, D[a][a]-d_beta
    Pnext, Dnext = [row[:] for row in P], [row[:] for row in D]
    Pnext[a][a], Dnext[a][a] = target, dtarget
    J, Jnext = inverse(P), inverse(Pnext)
    before, after = trace(product(J, D, J, D)), trace(product(Jnext, Dnext, Jnext, Dnext))
    positive = dC*dC/(C*C)+2*(1/C-1/Cnext)*product(dT, B, transpose(dT))[0][0]
    adverse = ((d_beta-dtarget)/Cnext)**2
    assert before-after == positive-adverse and positive >= 0
    delta, ratio, q, mu = target-P[a][a], C/Cnext, dC/C, D[a][a]-dtarget
    retained = deficit_charge = None
    if delta > 0:
        retained = (1-ratio**2)*(q+ratio*mu/(Cnext*(1-ratio**2)))**2
        retained += 2*(1/C-1/Cnext)*product(dT, B, transpose(dT))[0][0]
        deficit_charge = mu**2/(delta*(2*C+delta))
        assert before-after == retained-deficit_charge
    return {'P_next': Pnext, 'D_next': Dnext, 'B': B, 'T': T,
            'C': C, 'C_next': Cnext, 'beta': beta, 'd_beta': d_beta,
            'dC': dC, 'dT': dT, 'positive_face_action': positive,
            'regression_reader': (d_beta-dtarget)/Cnext, 'adverse_face_work': adverse,
            'AA_target_relative_tangent': mu, 'retained_coupled_square': retained,
            'deficit_reader_charge': deficit_charge,
            'gap': before-after}


def actual_aw_process_short():
    """CR25: short the FULL integrated process onto the original AW row.

    The normalized inverse AA entry is 8, rather than the unrelated full
    inverse row bound. This is an existing-process comparison, not a second
    covariance recursion. Its use below retains the unsplit receipt square.
    """
    from .lin_path_certificate import small_x_source_defect
    eps, _, b0, _ = small_x_source_defect()
    assert inverse(b0)[3][3] == 8
    xmin = F(1, 3000)
    full = F('.05')**2*xmin*(1-eps)/(8*(1+xmin)**2)
    q = full/2  # Q-q ee' >= Q/2; no inverse of a singular residual noise.
    assert q > F('5.2e-8')
    return {'full_short_lower': str(full), 'allocated_short': str(q),
            'allocated_short_lower': '0.000000052',
            'normalized_inverse_AA_entry': '8',
            'scope': 'existing regular real isotropic profile: h in [.004,.006], tau in [.02,12], sigma in [.05,4]; literal polynomial Q and PSD repair',
            'stationary_Q_AA_identity_used': False,
            'full_integrated_cross_covariance_retained': True,
            'uniform_signed_word_margin': None}


def residual_process_receipt_payment(X, process, D, aw_index):
    """CR44: price a NONZERO marginal receipt with the full residual Q.

    X=F P F', process=Q_ship-q_aw ee', D=F dP F' are the SAME
    preceding prediction already retained in CR37's L_rem. No process loss
    is added. The rank-one split is only an algebraic comparison of that
    full addition; its complementary half and transverse square remain.
    It is not an independently reachable covariance or a new recursion.
    """
    n, a = len(X), aw_index
    if not 0 <= a < n or D != transpose(D):
        raise ValueError('same symmetric prediction tangent and AW selector required')
    ldlt(X)
    ldlt(process)
    v, mu = X[a][a], D[a][a]
    w = [[X[i][a]] for i in range(n)]
    d = product(transpose(w), inverse(process), w)[0][0]
    allocation = [[x/(2*d) for x in row] for row in product(w, transpose(w))]
    remainder = add(process, allocation, -1)
    # w w'/d <= process by Cauchy in the actual process precision.
    assert is_psd(add(remainder, [[x/2 for x in row] for row in process], -1))
    ldlt(remainder)
    middle, terminal = add(X, allocation), add(X, process)
    J, Jmid, Jend = inverse(X), inverse(middle), inverse(terminal)
    fisher = lambda j: trace(product(j, D, j, D))
    paid = (4*d+v)/(v*(2*d+v)**2)
    s = v/(2*d+v)
    # ||(I-zz') X^-1/2 D X^-1/2 z||^2, evaluated without square roots.
    transverse = product([[D[a][j] for j in range(n)]], J,
                         [[D[i][a]] for i in range(n)])[0][0]/v-mu*mu/(v*v)
    assert transverse >= 0
    rank_loss = fisher(J)-fisher(Jmid)
    full_remainder = fisher(Jmid)-fisher(Jend)
    assert rank_loss == paid*mu*mu+2*s*transverse
    assert full_remainder >= 0
    return {'marginal_before_process': v, 'same_prefix_receipt': mu,
            'full_process_directional_inverse_charge': d,
            'paid_receipt_coefficient': paid,
            'rank_one_allocation': allocation,
            'retained_process_covariance': remainder,
            'retained_transverse_Fisher_loss': 2*s*transverse,
            'retained_full_process_Fisher_loss': full_remainder,
            'full_residual_process_Fisher_loss': rank_loss+full_remainder,
            'uniform_paid_coefficient': None,
            'uniform_receipt_domination_verified': False}


def prediction_face_fisher_balance(P, f, noise, D, target, aw_index, q_aw,
                                   zero_gap_active=False):
    """CR26: SAME process/face tangent, including a nonzero AA receipt.

    The literal process is split algebraically as Q-q ee' then q ee'.
    No process loss is added twice and no intermediate covariance is stored
    or used by the estimator. The actual row/cross blocks stay in f and Q.
    A zero-gap active formula is a linear extension valid only on mu<=0;
    callers must retain that cone, not apply it to all root coordinates.
    """
    a, q = aw_index, F(q_aw)
    n = len(P)
    if not 0 <= a < n or q <= 0:
        raise ValueError('positive allocated AW process short required')
    if any(f[a][j] for j in range(n) if j != a):
        raise ValueError('literal autonomous AW prediction row required')
    allocated = zeros(n, n)
    allocated[a][a] = q
    if not is_psd(add(noise, allocated, -2)):
        raise ValueError('actual full process must supply at least twice the allocated AW short')
    px, dx = product(f, P, transpose(f)), product(f, D, transpose(f))
    pm = add(px, noise)
    intermediate = add(pm, allocated, -1)
    ldlt(px)
    ldlt(intermediate)
    face = scalar_aw_face_fisher_balance(pm, dx, target, a,
                                        zero_gap_active=zero_gap_active)
    def fisher(p, d):
        j = inverse(p)
        return trace(product(j, d, j, d))
    residual = fisher(px, dx)-fisher(intermediate, dx)
    assert residual >= 0 and fisher(px, dx) == fisher(P, D)
    c0, A = face['C']-q, face['C_next']
    dc, mu, dt, B = face['dC'], face['AA_target_relative_tangent'], face['dT'], face['B']
    assert c0 > q and A > c0
    regression = 2*(1/c0-1/A)*product(dt, B, transpose(dt))[0][0]
    joint = residual+(dc/c0)**2-((dc-mu)/A)**2+regression
    ratio = c0/A
    denominator = (A-c0)*(A+c0)
    coupled_row = dc/c0+ratio*mu/(A*(1-ratio**2))
    retained = (1-ratio**2)*coupled_row**2+regression
    charge = mu**2/denominator
    direct = fisher(P, D)-fisher(face['P_next'], face['D_next'])
    assert direct == joint == residual+retained-charge
    payment = residual_process_receipt_payment(px, add(noise, allocated, -1), dx, a)
    assert payment['full_residual_process_Fisher_loss'] == residual
    unpaid = 1/denominator-payment['paid_receipt_coefficient']
    assert direct == (payment['retained_full_process_Fisher_loss']
                      + payment['retained_transverse_Fisher_loss']+retained-unpaid*mu*mu)
    return {**face, 'C_before_allocated_noise': c0,
            'allocated_process_short': q, 'remaining_process_Fisher_loss': residual,
            'effective_gap': A-c0, 'receipt_denominator': denominator,
            'effective_ratio': ratio, 'coupled_receipt_row': coupled_row,
            'coupled_conditional_action': retained,
            'joint_receipt_charge': charge, 'joint_gap': direct,
            'residual_process_receipt_payment': payment,
            'receipt_weight_after_residual_payment': unpaid,
            'zero_gap_active_cone': 'mu<=0' if zero_gap_active else None,
            'uniform_signed_word_margin_verified': False}


def scalar_acc_conditional_fisher(P, D, h, noise, aw_index):
    """CR18--CR20: actual AW-observing row, one shared covariance tangent.

    The planar odd acc row has h_a=1. This is an exact completion of its
    existing Fisher loss, not an independent bound on the floor receipt.
    Row/noise/auxiliary variations outside the CR12 partial fibre remain
    in the full correlated gap. No uniform innovation bound is asserted.
    """
    n, a, noise = len(P), aw_index, F(noise)
    if (not 0 <= a < n or D != transpose(D) or noise <= 0
            or len(h) != 1 or len(h[0]) != n or h[0][a] != 1):
        raise ValueError('actual scalar AW-observing acc row and symmetric tangent required')
    ldlt(P)
    o = [i for i in range(n) if i != a]
    B, dB = [[P[i][j] for j in o] for i in o], [[D[i][j] for j in o] for i in o]
    Binv = inverse(B)
    T = product([[P[a][j] for j in o]], Binv)
    dT = product(add([[D[a][j] for j in o]], product(T, dB), -1), Binv)
    beta = product(T, B, transpose(T))[0][0]
    d_beta = product(T, dB, transpose(T))[0][0]+2*product(dT, B, transpose(T))[0][0]
    C, dC = P[a][a]-beta, D[a][a]-d_beta
    k = add([[h[0][j] for j in o]], T)
    r = product(k, B, transpose(k))[0][0]
    s = noise+C+r
    u = add(product(dB, transpose(k)), product(B, transpose(dT)))
    t = product(dT, B, transpose(k))[0][0]/C
    v = dC+C*t
    Ku = add([[2*x/s for x in row] for row in Binv],
             [[x/(s*s) for x in row] for row in product(transpose(k), k)], -1)
    ldlt(Ku)
    residual = add(u, [[x*v/(2*s-r) for x in row]
                      for row in product(B, transpose(k))], -1)
    eta = 2*C*(2*s-r-C)/(s*(2*s-r))
    omega = 1-eta
    assert 0 < eta < 1 and omega > 0
    square = product(transpose(residual), Ku, residual)[0][0]
    loss = scalar_aw_innovation_reader(P, D, h, noise, a)['same_operation_Fisher_loss']
    assert loss == square+eta*(v/C)**2
    rho = noise/(noise+C)
    return {'B': B, 'T': T, 'C': C, 'dC': dC, 'dT': dT,
            'alignment': k, 'alignment_variance': r, 'innovation': s,
            'u': u, 'conditional_residual': v, 'linked_regression_reader': t,
            'positive_matrix': Ku, 'positive_residual': residual,
            'positive_square': square, 'paid_fraction': eta,
            'remaining_weight': omega, 'same_operation_Fisher_loss': loss,
            'C_after': rho*C, 'T_after': add([[rho*x for x in row] for row in k],
                                            [[h[0][j] for j in o]], -1),
            'uniform_paid_fraction_verified': False}


def covariance_word_signed_matrix(boundaries, prefixes, active_faces, aw_index,
                                  acc_face_pairs=(), zero_gap_active_faces=()):
    """Verify OF4 from supplied SAME-word covariance/tangent prefix maps.

    prefixes[i][j] is the entire actual D_i for root coordinate j; after a
    face it includes AA deletion, not a frozen congruence product. Ordinary
    steps must be the proved fixed-operand process/Joseph/congruence steps.
    This small exact-algebra helper does not construct or qualify a shipping
    word, an inherited root lift, magnetic service, or a uniform margin.
    """
    if len(boundaries) < 2 or len(prefixes) != len(boundaries):
        raise ValueError('linked covariance and tangent boundaries required')
    d = len(prefixes[0])
    if not d or any(len(ds) != d for ds in prefixes):
        raise ValueError('one common root-coordinate basis required')
    faces = set(active_faces)
    zero_faces = set(zero_gap_active_faces)
    if not zero_faces <= faces:
        raise ValueError('zero-gap active branches must be actual supplied face maps')
    if any(i < 0 or i >= len(boundaries)-1 for i in faces):
        raise ValueError('active face outside word')
    grams = []
    for P, ds in zip(boundaries, prefixes):
        ldlt(P)
        if any(D != transpose(D) or len(D) != len(P) for D in ds):
            raise ValueError('symmetric same-dimension covariance tangents required')
        J = inverse(P)
        JD = [product(J, D) for D in ds]
        grams.append([[trace(product(x, y)) for y in JD] for x in JD])
    positive, readers, ordinary = zeros(d, d), [], []
    receipt_action, receipt_rows, conditional_rows = zeros(d, d), [], []
    for i in range(len(boundaries)-1):
        if i not in faces:
            loss = add(grams[i], grams[i+1], -1)
            if not is_psd(loss):
                raise ValueError('ordinary step lacks the stated Fisher base-loss qualification')
            positive = add(positive, loss)
            receipt_action = add(receipt_action, loss)
            ordinary.append(loss)
            continue
        fs = [scalar_aw_face_fisher_balance(boundaries[i], D,
              boundaries[i+1][aw_index][aw_index], aw_index,
              zero_gap_active=i in zero_faces) for D in prefixes[i]]
        if any(f['P_next'] != boundaries[i+1] or f['D_next'] != prefixes[i+1][j]
               for j, f in enumerate(fs)):
            raise ValueError('face must retain cross covariance and its actual target-fixed derivative')
        f0 = fs[0]
        face_action = [[f['dC']*g['dC']/f0['C']**2
                       + 2*(1/f0['C']-1/f0['C_next'])*product(f['dT'], f0['B'], transpose(g['dT']))[0][0]
                       for g in fs] for f in fs]
        positive = add(positive, face_action)
        readers.append([f['regression_reader'] for f in fs])
        ratio = f0['C']/f0['C_next']
        face_unsplit_positive = [
            [(1-ratio**2)*f['dC']*g['dC']/f0['C']**2
             + 2*(1/f0['C']-1/f0['C_next'])
             * product(f['dT'], f0['B'], transpose(g['dT']))[0][0]
             for g in fs] for f in fs]
        receipt_action = add(receipt_action, face_unsplit_positive)
        receipt_rows.append([f['AA_target_relative_tangent']/f0['C_next'] for f in fs])
        conditional_rows.append([ratio*f['dC']/f0['C'] for f in fs])
    adverse = product(transpose(readers), readers) if readers else zeros(d, d)
    gap = add(positive, adverse, -1)
    assert gap == add(grams[0], grams[-1], -1)
    if receipt_rows:
        coupled = add(add(receipt_action,
                          product(transpose(conditional_rows), receipt_rows)),
                      product(transpose(receipt_rows), conditional_rows))
        coupled = add(coupled, product(transpose(receipt_rows), receipt_rows), -1)
        assert coupled == gap
    # Pair a floor only with its actual subsequent acc substep. Between
    # them CR7 S corrections and qualified block resets preserve C,dC.
    # Validate that constraint on every SAME prefix; never free the receipt.
    paired_action, paired_readers = [row[:] for row in positive], [row[:] for row in readers]
    weights, paired_corrections = [F(1) for _ in readers], set()
    pair_data = []
    for face, corr, h, noise in acc_face_pairs:
        if face not in faces or not face < corr < len(boundaries)-1 or corr in paired_corrections:
            raise ValueError('disjoint actual floor-to-acc chronology required')
        paired_corrections.add(corr)
        slot = sorted(faces).index(face)
        ps = [scalar_acc_conditional_fisher(boundaries[corr], D, h, noise, aw_index)
              for D in prefixes[corr]]
        fs = [scalar_aw_face_fisher_balance(boundaries[face], D,
              boundaries[face+1][aw_index][aw_index], aw_index,
              zero_gap_active=face in zero_faces) for D in prefixes[face]]
        if any(p['C'] != f['C_next'] or p['dC'] != -f['d_beta'] for p, f in zip(ps, fs)):
            raise ValueError('same conditional covariance/tangent must survive floor-to-acc interval')
        p0, eta, omega = ps[0], ps[0]['paid_fraction'], ps[0]['remaining_weight']
        lp = add(grams[corr], grams[corr+1], -1)
        paid = [[product(transpose(p['positive_residual']), p0['positive_matrix'],
                         q['positive_residual'])[0][0]
                 + eta*(p['linked_regression_reader']-f['regression_reader'])
                      *(q['linked_regression_reader']-g['regression_reader'])
                 for q, g in zip(ps, fs)] for p, f in zip(ps, fs)]
        if lp != paid:
            raise ValueError('acc loss must be its complete fixed-operand same-prefix Fisher loss')
        positive_part = [[product(transpose(p['positive_residual']), p0['positive_matrix'],
                                  q['positive_residual'])[0][0]
                          + eta/omega*p['linked_regression_reader']*q['linked_regression_reader']
                          for q in ps] for p in ps]
        paired_action = add(add(paired_action, lp, -1), positive_part)
        paired_readers[slot] = [f['regression_reader']+eta/omega*p['linked_regression_reader']
                                for p, f in zip(ps, fs)]
        weights[slot] = omega
        pair_data.append({'face': face, 'correction': corr, 'paid_fraction': eta,
                          'remaining_weight': omega, 'alignment_variance': p0['alignment_variance'],
                          'innovation': p0['innovation']})
    paired_adverse = zeros(d, d)
    for reader, weight in zip(paired_readers, weights):
        paired_adverse = add(paired_adverse, product(transpose([reader]), [reader]), weight)
    assert add(paired_action, paired_adverse, -1) == gap
    return {'root_metric': grams[0], 'positive_action': positive,
            'ordinary_losses': ordinary, 'face_reader': readers,
            'unsplit_receipt_action': receipt_action,
            'normalized_receipt_rows': receipt_rows,
            'conditional_receipt_cross_rows': conditional_rows,
            'adverse_face_work': adverse, 'signed_gap': gap,
            'acc_paired_positive_action': paired_action,
            'acc_paired_face_reader': paired_readers,
            'acc_paired_reader_weights': weights, 'acc_face_payments': pair_data}


def covariance_partial_word(P, root_tangents, operations, aw_index):
    """Construct all prefixes from the CR12 fixed-operand partial recursion.

    Operands must come from ONE qualified nominal shipping history. F/Q are
    the full odd prediction including integrated cross blocks and repairs;
    scalar h/R are the applied odd acc/S/mag channels. This routine does not
    select coefficients, qualify an inherited image, or assert reachability.
    Unlike supplied-prefix checking, every preceding AA deletion is propagated.
    General coefficient/target/mean/auxiliary variations require CR2's full
    lift and are not silently set to zero by this routine.
    """
    boundaries, prefixes, faces = [[row[:] for row in P]], [root_tangents], []
    zero_faces, branch_guards = [], []
    a, n = aw_index, len(P)
    if not 0 <= a < n or not root_tangents:
        raise ValueError('nonempty common covariance root and AW selector required')
    for event in operations:
        P, ds = boundaries[-1], prefixes[-1]
        kind = event['kind']
        if kind == 'prediction':
            f, q = event['F'], event['Q']
            if any(f[a][j] for j in range(n) if j != a) or not is_psd(q):
                raise ValueError('literal autonomous AW row and PSD applied process required')
            Pnext = add(product(f, P, transpose(f)), q)
            dn = [product(f, D, transpose(f)) for D in ds]
        elif kind == 'correction':
            h, noise = event['h'], F(event['noise'])
            out = correction(P, h, [[noise]], [[F(0)]])
            A = add(identity(n), product(out['K'], h), -1)
            Pnext, dn = out['C'], [product(A, D, transpose(A)) for D in ds]
        elif kind == 'reset':
            g = event['G']
            Pnext, dn = product(g, P, transpose(g)), [
                product(g, D, transpose(g)) for D in ds]
        elif kind == 'floor':
            target = F(event['target'])
            Pnext, dn = [row[:] for row in P], [[row[:] for row in D] for D in ds]
            if target > P[a][a]:
                faces.append(len(boundaries)-1)
                Pnext[a][a] = target
                for D in dn:
                    D[a][a] = F(0)
            elif target == P[a][a]:
                branch = event.get('zero_gap_branch')
                if branch not in (None, 'active', 'inactive'):
                    raise ValueError('zero-gap directional branch must be active or inactive')
                if branch is None and any(D[a][a] for D in ds):
                    raise ValueError('zero-gap crossing needs its actual directional map')
                if branch is not None:
                    branch_guards.append({'operation': len(boundaries)-1,
                                          'receipt_row': [D[a][a] for D in ds],
                                          'relation': '<=0' if branch == 'active' else '>=0'})
                if branch == 'active':
                    faces.append(len(boundaries)-1)
                    zero_faces.append(len(boundaries)-1)
                    for D in dn:
                        D[a][a] = F(0)
        else:
            raise ValueError('unsupported covariance partial operation')
        boundaries.append(Pnext)
        prefixes.append(dn)
    pairs = []
    for face in faces:
        for j in range(face+1, len(operations)):
            event = operations[j]
            if event['kind'] == 'correction':
                if event['h'][0][a] == 1:
                    pairs.append((face, j, event['h'], event['noise']))
                    break
                if event['h'][0][a] != 0:
                    break
            elif event['kind'] == 'reset':
                g = event['G']
                if (g[a] != identity(n)[a] or any(g[i][a] for i in range(n) if i != a)):
                    break
            else:
                break  # no pairing across prediction, another floor or a proof reset
    out = covariance_word_signed_matrix(boundaries, prefixes, faces, a, pairs, zero_faces)
    # Allocate each actual preceding prediction once, retaining every
    # completed q/mu square on the common root. Do not release coordinates.
    joint_action = zeros(len(root_tangents), len(root_tangents))
    joint_readers, joint_weights, joint_data, used = [], [], [], set()
    receipt_process_payments = []
    q = F(actual_aw_process_short()['allocated_short'])
    for face in faces:
        if face == 0 or operations[face-1]['kind'] != 'prediction':
            continue  # proof boundary is not a process or a floor event
        pred = face-1
        allocation = zeros(n, n)
        allocation[a][a] = 2*q
        if not is_psd(add(operations[pred]['Q'], allocation, -1)):
            continue  # outside CR25's profile: retain original exact OF4
        used.update((pred, face))
        ps = [prediction_face_fisher_balance(
            boundaries[pred], operations[pred]['F'], operations[pred]['Q'], D,
            operations[face]['target'], a, q, face in zero_faces) for D in prefixes[pred]]
        p0 = ps[0]
        if any(p['D_next'] != prefixes[face+1][j] or p['P_next'] != boundaries[face+1]
               for j, p in enumerate(ps)):
            raise ValueError('process and face must use their complete common prefix')
        # Polarize the residual PSD addition by its actual two Fisher Grams.
        px = product(operations[pred]['F'], boundaries[pred], transpose(operations[pred]['F']))
        intermediate = [row[:] for row in boundaries[face]]
        intermediate[a][a] -= q
        jx, jr = inverse(px), inverse(intermediate)
        dx = prefixes[face]
        residual = [[trace(product(jx, d, jx, e))-trace(product(jr, d, jr, e))
                     for e in dx] for d in dx]
        positive = [[(1-p0['effective_ratio']**2)*p['coupled_receipt_row']*r['coupled_receipt_row']
                     + 2*(1/p0['C_before_allocated_noise']-1/p0['C_next'])
                     * product(p['dT'], p0['B'], transpose(r['dT']))[0][0]
                     for r in ps] for p in ps]
        joint_action = add(joint_action, add(residual, positive))
        joint_readers.append([p['AA_target_relative_tangent'] for p in ps])
        joint_weights.append(1/p0['receipt_denominator'])
        receipt_process_payments.append(p0['residual_process_receipt_payment']['paid_receipt_coefficient'])
        joint_data.append({'prediction': pred, 'face': face,
                           'allocated_process_short': q,
                           'effective_gap': p0['effective_gap'],
                           'receipt_denominator': p0['receipt_denominator']})
    for i in range(len(operations)):
        if i in used:
            continue
        if i in faces:
            slot = faces.index(i)
            fs = [scalar_aw_face_fisher_balance(boundaries[i], D, operations[i]['target'], a,
                                                zero_gap_active=i in zero_faces) for D in prefixes[i]]
            positive = [[p['dC']*r['dC']/fs[0]['C']**2
                         + 2*(1/fs[0]['C']-1/fs[0]['C_next'])
                         * product(p['dT'], fs[0]['B'], transpose(r['dT']))[0][0]
                         for r in fs] for p in fs]
            joint_action = add(joint_action, positive)
            joint_readers.append(out['face_reader'][slot])
            joint_weights.append(F(1))
            receipt_process_payments.append(F(0))  # no invented preceding process
        else:
            ji, jn = inverse(boundaries[i]), inverse(boundaries[i+1])
            loss = [[trace(product(ji, d, ji, e))-trace(product(jn, dn, jn, en))
                     for e, en in zip(prefixes[i], prefixes[i+1])]
                    for d, dn in zip(prefixes[i], prefixes[i+1])]
            joint_action = add(joint_action, loss)
    adverse = zeros(len(root_tangents), len(root_tangents))
    for row, weight in zip(joint_readers, joint_weights):
        adverse = add(adverse, product(transpose([row]), [row]), weight)
    assert add(joint_action, adverse, -1) == out['signed_gap']
    paid_action = zeros(len(root_tangents), len(root_tangents))
    for row, payment in zip(joint_readers, receipt_process_payments):
        paid_action = add(paid_action, product(transpose([row]), [row]), payment)
    residual_action = add(joint_action, paid_action, -1)
    assert is_psd(residual_action)
    residual_weights = [weight-payment for weight, payment
                        in zip(joint_weights, receipt_process_payments)]
    residual_charge = zeros(len(root_tangents), len(root_tangents))
    for row, weight in zip(joint_readers, residual_weights):
        residual_charge = add(residual_charge, product(transpose([row]), [row]), weight)
    assert add(residual_action, residual_charge, -1) == out['signed_gap']
    out.update({'boundaries': boundaries, 'prefixes': prefixes, 'active_faces': faces,
                'zero_gap_branch_guards': branch_guards,
                'process_face_positive_action': joint_action,
                'process_face_receipt_rows': joint_readers,
                'process_face_receipt_weights': joint_weights,
                'process_face_payments': joint_data,
                'process_face_receipt_process_paid_weights': receipt_process_payments,
                'process_face_retained_positive_action': residual_action,
                'process_face_receipt_weights_after_process_payment': residual_weights,
                'inherited_causal_image_qualified': False})
    return out


def receipt_kernel_nuisance_margin():
    """CR16: signed WHOLE-word reserve on the actual nuisance/receipt kernel.

    All face work is nonnegative on that kernel, so the first process loss
    can be retained. Zero-gap events also require zero target-relative AA
    tangent: active-receipt rows alone do not imply that crossing condition.
    No reserve is asserted off this intersection; its actual rank,
    origin and 17-second regular A21 activation are separate qualifications.
    """
    from .lin_path_certificate import small_x_source_defect
    from .nuisance_upper_certificate import bounds
    from .root_covariance_certificate import process_floors
    eps, _, normalized, _ = small_x_source_defect()
    inverse_bound = max(sum(abs(x) for x in row) for row in inverse(normalized))
    hmin, hmax = F('.004'), F('.006')
    lin = ((1-eps)*F('.05')**2*(hmin/12)*hmin**6
           / ((1+hmax/F('.02'))**2*inverse_bound))
    q = min(lin, process_floors()[1])/F('1.01')**2
    m = max(bounds()[-1])
    delta = q/(m+q)
    reserve = delta*(2-delta)
    assert F('1.37935e-37') < reserve < F('1.37936e-37')
    return {'q_n': str(q), 'm_n': str(m), 'signed_kernel_reserve': str(reserve),
            'signed_kernel_reserve_lower': '1/10000000000000000000000000000000000000',
            'scope': 'CR12 pure odd covariance causal image, nuisance-supported root and zero target-relative receipt at every active or zero-gap face; nonempty word includes first prediction at root after 17 s regular default A21',
            'zero_gap_target_relative_tangent_required': True,
            'actual_kernel_dimension_verified': False,
            'whole_odd_gap_verified': False,
            'full_homogeneous_gap_verified': False}


def coupled_receipt_schur(K, L, Qfull, embedding):
    """CR17a exact common-root elimination; ordinary inverses only.

    L consists of independent ACTUAL receipt rows; embedding*L retains all
    faces. Qfull contains their actual conditional cross rows. Supplying this
    algebra neither proves those origins nor K>0 at a uniform root reserve.
    """
    ldlt(K)
    Ki = inverse(K)
    Z = product(L, Ki, transpose(L))
    ldlt(Z)
    Zi = inverse(Z)
    Y = product(Ki, transpose(L), Zi)
    shorted = add(Ki, product(Y, L, Ki), -1)
    Q = product(transpose(embedding), Qfull)
    remainder = add(add(add(add(Zi, product(Q, Y)),
                            product(transpose(Y), transpose(Q))),
                        product(transpose(embedding), embedding), -1),
                    product(Q, shorted, transpose(Q)), -1)
    return {'receipt_remainder': remainder, 'receipt_lift': Y,
            'kernel_inverse': shorted, 'conditional_rows': Q,
            'uniform_actual_remainder_verified': False}


def odd_covariance_gap_scope_check():
    """One exact D-class implication check, NOT a shipping execution (OF7).

    It even retains the fixed-root AW marginal tangent D_aa=0. The correction
    precedes prediction and its pending floor, as the shipping chronology can.
    F=I and Q=I/1000 are a stated formal relaxation, not the OU process.
    """
    P = [[F(1), F(0), F(1)], [F(0), F(1), F(1)], [F(1), F(1), F(3)]]
    H, R = [[F(0), F(1), F(0)]], [[F(1)]]
    S = add(product(H, P, transpose(H)), R)
    A = add(identity(3), product(P, transpose(H), inverse(S), H), -1)
    Pc = product(A, P)
    Q = [[F(i == j, 1000) for j in range(3)] for i in range(3)]
    Pp = add(Pc, Q)
    Pn = [row[:] for row in Pp]
    Pn[2][2] = F(3)
    # Coordinates (dB11,dB12,dB22,dT1,dT2), inherited fixed P_aa=3.
    ds = []
    for X, u in (([[F(1), F(0)], [F(0), F(0)]], [[F(0), F(0)]]),
                 ([[F(0), F(1)], [F(1), F(0)]], [[F(0), F(0)]]),
                 ([[F(0), F(0)], [F(0), F(1)]], [[F(0), F(0)]]),
                 (zeros(2, 2), [[F(1), F(0)]]),
                 (zeros(2, 2), [[F(0), F(1)]])):
        b = add(product(X, [[F(1)], [F(1)]]), transpose(u))
        ds.append([X[0]+b[0], X[1]+b[1], [b[0][0], b[1][0], F(0)]])
    dc = [product(A, D, transpose(A)) for D in ds]
    dn = [scalar_aw_face_fisher_balance(Pp, D, F(3), 2)['D_next'] for D in dc]
    word = covariance_word_signed_matrix([P, Pc, Pp, Pn], [ds, dc, dc, dn], [2], 2)
    ldlt(word['positive_action'])
    reader = product(word['face_reader'], inverse(word['positive_action']), transpose(word['face_reader']))[0][0]
    direction = [[F(-10)], [F(3)], [F(1)], [F(0)], [F(0)]]
    gap = product(transpose(direction), word['signed_gap'], direction)[0][0]
    assert reader > 1 and gap == -F(691428110973, 567390082009)
    return {'classification': 'D_SUFFICIENT_BOUND_FAILURE', 'shipping_counterexample': False,
            'root_AW_marginal_tangent_zero': True, 'strict_positive_process': True,
            'ordinary_and_face_chronology_retained': True,
            'relative_reader_threshold': '1', 'relative_reader_exact': str(reader),
            'negative_direction': ['-10', '3', '1', '0', '0'], 'signed_gap': str(gap),
            'relaxation': 'formal three-state covariance word; F=I,Q=I/1000, not shipping-reached or service-admitted'}


def certificate():
    return {
        'qualification': 'OU3_MEASUREMENT_FRAME_V2',
        'result_type': 'PROVED — analytical moving-frame/AW-shear identities and qualified planar magnetic-loss implication',
        'scope': 'regular real-operation complete 21-state CoG profile; no estimator modification',
        'positive_magnetic_loss_scope': 'nominal planar invariant stratum, even mean and parity-block covariance; not all nominal perturbations',
        'profile': 'planar wrapper sigma_a=.2, sigma_g=.00135, sigma_m=.8, b0=1e-10, adaptive S; no AtomS3R transfer',
        'frame': 'T=diag(R,R,I,I,I,I,R), R=actual nominal world-to-body rotation',
        'actual_P_gain_Joseph_congruence': True,
        'held_mask_congruence': True,
        'mag_row': 'Hbar_mag=[-[Bref]x,0,0,0,0,0,0]',
        'acc_row': 'Hbar_acc=[-[aw-g]x,0,0,0,0,I,I]',
        'attitude_row_derivative': 'exactly zero in the moving frame for both measurements',
        'remaining_acc_row_derivative': 'dHbar_acc=[-[daw]x,0,0,0,0,0,0]',
        'remaining_mag_row_derivative': 'dHbar_mag=[-[dBref]x,0,0,0,0,0,0]',
        'anisotropic_noise_port': 'dRbar=R^T dRnoise R-[omega]Rbar+Rbar[omega]',
        'isotropic_fixed_noise_rotation_port': '0',
        'covariance_frame_derivative': 'dPbar=T^T dP T-Omega Pbar+Pbar Omega',
        'shear_frame_connection': 'etabar=T^T eta-Pbar Omega Jbar ebar',
        'mag_mean_connection_port': '-Kbar([omega] rbar+Hbar Omega ebar)',
        'gain_information_identity': 'Kbar^T Cbar^-1 Kbar=Rbar^-1-Sbar^-1',
        'mag_connection_energy': 'v^T (Rbar^-1-Sbar^-1) v, v=[omega]rbar+Hbar Omega ebar',
        'acc_covariance_port_bound': '||qP||_C^2 <= 4 tr(C dH^T (R^-1-S^-1) dH) <= 8 lambda_max(C_theta,theta)||daw||^2/rmin',
        'coefficient_8_is_uniform_algebraic_constant': True,
        'AW_shear_proof': 'app:aw-shear-loss',
        'AW_shear': 'L=I-E_aw[aw]x E_theta^T, L^-1=I+E_aw[aw]x E_theta^T',
        'AW_shear_acc_row': '[+[g]x,0,0,0,0,I,I]',
        'AW_row_covariance_port_is_exact_coboundary': True,
        'AW_shear_connection_derivative_retained': True,
        'AW_shear_correction_balance': 'Delta W=-||H eta-m||_(S^-1)^2-lambda L_P+||m||_(R^-1)^2; m is linked, not independent noise',
        'AW_shear_post_correction_jump': 'L(aw+)L(aw)^-1=I+N(K_aw r); dJump=N(dK_aw r+K_aw dr)',
        'AW_shear_uniform_port_absorption_verified': False,
        'AW_frame_word_proof': 'app:aw-frame-word',
        'internal_frame_connections_cancel_in_complete_word': True,
        'frame_cancellation_requires_generated_covariance_suffix_score': True,
        'AW_frame_storage_difference': 'Phi_N-Phi_0; Phi=2 h^T daw+daw^T Q daw, Q positive definite',
        'AW_frame_endpoint_precision': 'J_aw,aw is conditional precision, not inverse marginal AW covariance',
        'AW_frame_endpoint_absorption_test': 'R_c>0 and Q_N^-1-Z R_c^-1 Z^T>=0; actual linked root maps, uniform margin OPEN',
        'AW_frame_endpoint_absorption_verified': False,
        'actual_AW_sync_conditional_precision_nonincrease': True,
        'actual_AW_sync_connection_square_nonincrease': True,
        'conditional_AW_joint_storage_decomposition': True,
        'conditional_AW_storage_and_mixed_reader_frame_invariant': True,
        'non_AW_measurement_conditional_storage_invariant': True,
        'non_AW_measurement_endpoint_square_matrix_decreases': True,
        'linked_OU_conditional_precision_and_connection_balance': True,
        'qualified_AW_conditional_precision_ceiling': qualified_aw_precision_ceiling(),
        'conditional_AW_process_coercivity': conditional_aw_process_coercivity(),
        'conditional_AW_loss_proof': 'app:aw-conditional-loss',
        'odd_covariance_complete_gap': {
            'proof': 'app:odd-covariance-signed-gap',
            'result_type': 'PROVED analytical restriction and exact Schur obstruction; uniform gap OPEN',
            'scope': 'regular central-planar A21 continuation after refinement, fixed delivered/private inherited history and target; actual causal root image must be retained',
            'ambient_symmetric_dimension': 45,
            'non_AW_generated_self_ports_zero': True,
            'E_L_and_source_score_priced_separately': False,
            'actual_AW_face_derivative_retained': True,
            'signed_matrix': 'lambda*(H_W-R_W^T R_W), evaluated on actual covariance prefixes',
            'remaining_operator': 'stacked d(T B T^T)/(C+Delta) at actual active AW_y faces',
            'magnetic_and_S_loss': 'same-operation directional Fisher Grams, transported through actual prefix including AW deletion',
            'relative_Schur_threshold': 'I-R_W*(H_W-c*H_0)^-1*R_W^T >= 0, H_W-c*H_0 > 0',
            'uniform_relative_reader_margin': None,
            'uniform_odd_covariance_gap': None,
            'full_even_odd_source_cross_blocks_discarded': False,
            'scope_check': odd_covariance_gap_scope_check(),
        },
        'causal_aw_deficit_gap': {
            'proof': 'app:correlated-complete-gap, CR4--CR43',
            'result_type': 'joint actual-process/nonzero-receipt inequality and scalar directional branches; full AW absorption OPEN',
            'scope': 'same qualified regular planar odd-covariance fibre; full correlated cross/target ports retained',
            'actual_process_AA_row': 'p_minus=phi^2*p+actual_q_aa',
            'queued_process_target_lag_retained': True,
            'stationary_Q_aa_identity_assumed_on_all_source_branches': False,
            'fixed_private_causal_receipt': 'mu=-sum actual_prediction_weights*dI_a after an actual active fixed-target face',
            'same_correction_sharp_reader': '(dI_cov)^2 <= I_a*(2P_aa-I_a)*L_P',
            'signed_row_noise_target_tuner_ports_retained': True,
            'active_face_completion': 'q^2-(C/Cplus*q-mu/Cplus)^2+regression_loss',
            'negative_completed_reader': 'mu^2/[Delta*(2C+Delta)], Delta>0 only',
            'non_AW_payment': 'I_a=T*(B-Bplus)*T^T; dI_a and dC share actual dB,dT',
            'common_tangent_unsplit_matrix': 'H_mu+Q_mu^T L_mu+L_mu^T Q_mu-L_mu^T L_mu; every row from the same actual prefix',
            'common_tangent_matrix_identity_verified': True,
            'partial_origin_zero_energy_kernel_test_verified': True,
            'partial_origin_uniform_parameterization_verified': False,
            'nuisance_receipt_kernel': receipt_kernel_nuisance_margin(),
            'full_cross_remainder': 'CR17: A_c-X_c^T N_c^-1 X_c >= 0, with all root metric cross entries retained',
            'receipt_only_remainder': 'CR17a: Z^-1+QY+Y^T Q^T-E^T E-Q K_sharp Q^T >= 0; actual dependent receipts retained',
            'receipt_only_inverse_uniformly_qualified': False,
            'receipt_only_remainder_uniformly_verified': False,
            'actual_process_face_nonzero_receipt_payment': {
                'proof': 'CR25--CR28',
                'process_short': actual_aw_process_short(),
                'common_root_joint_signed_inequality_verified': True,
                'joint_completion': 'L_rem+(1-a_star^2)*(dC/C_star+a_star*mu/[A*(1-a_star^2)])^2+2*(1/C_star-1/A)*dT B dT^T-mu^2/[(Delta+q_aw)*(A+C_star)]',
                'same_complete_matrix': 'G_OO/lambda=H_pmu-L_pmu^T D_pmu L_pmu; preceding process allocated once',
                'retained_receipt_square_discarded': False,
                'scalar_zero_gap_branch_accounting_verified': True,
                'zero_gap_active_guard': 'actual prefix mu<=0: AA deletion; mu>=0: identity; maps agree at mu=0',
                'actual_causal_image_or_rank_assumed': False,
                'original_CR17a_matrix_sign_changed': False,
                'relative_reader_threshold': 'D_pmu^-1-L_pmu*(H_pmu-c_O H_0)^-1*L_pmu^T>=0; inverse qualification and cone feasibility required',
                'uniform_receipt_remainder_verified': False,
                'uniform_signed_word_margin': None,
                'odd_AW_blocker_resolved': False,
            },
            'constrained_dangerous_direction': {
                'proof': 'CR29--CR35',
                'scope': 'one fixed qualified word on a justified closed linear history-tangent slice of CR12, intersected with its actual scalar face cones; no uniform image parameterization assumed',
                'root_normalized_Euler_equation': '(H-L^T D L-r H_0)y+U^T nu-V^T kappa=0; actual guards and complementarity retained',
                'fixed_word_zero_reserve_inverse_justified': True,
                'common_positive_reserve_inverse_justified': False,
                'feasible_receipt_eigenmode': 'A_F t=chi t; A_F=D^1/2 L_F H_F^-1 L_F^T D^1/2; y=H_F^-1 L_F^T D^1/2 t/chi; actual image and all cone guards required',
                'nonpositive_gap_threshold': 'chi>=1 on a reconstructed feasible mode',
                'balanced_shipping_root': 'D0=H_F^-1 sum_i D_ii mu_i prefix_i^*(e_a e_a^T); mu_i=e_a^T prefix_i(D0)e_a, with actual regression rows for unpaired faces',
                'integrated_dB_dT_dC_differentiation_retained': True,
                'same_transported_S_acc_magnetic_columns_retained': True,
                'zero_signed_gap_implies_zero_Fisher_action': False,
                'balanced_mode_has_nonzero_receipt_and_positive_first_process_loss': True,
                'complete_zero_action_kernel_excludes_balanced_mode': False,
                'uniform_normalized_history_image_compactness_assumed': False,
                'failure_class': 'D_SUFFICIENT_PROOF_INFERENCE_FAILURE',
                'shipping_balance_eigenmode_excluded_or_exhibited': False,
                'uniform_receipt_margin': None,
                'odd_AW_blocker_resolved': False,
                'complete_homogeneous_gap_proved': False,
            },
            'joint_process_face_acc_receipt_reduction': {
                'proof': 'CR36--CR40',
                'scope': 'actual prediction/active scalar floor/applied acc triples on CR12; only qualified non-AW S/block resets between floor and acc; all other operations unchanged',
                'same_receipt_reader_retained': True,
                'regression_transport': 'ell=U_T B_acc k^T; z=dT_face ell; kappa=ell^T B_face^-1 ell; 0<=kappa<=r',
                'same_root_positive_squares_retained': True,
                'reduced_inverse_receipt_weight': 'dbar^-1=A^2-C_star^2+eta*A^2/[omega+eta*kappa/(2*gamma*A^2)]; gamma=1/C_star-1/A>0',
                'paired_receipt_weight_strictly_reduced': True,
                'zero_kappa_division_used': False,
                'conditional_B_reader_constrained_to_zero': False,
                'complete_signed_matrix': 'G_OO/lambda=Hbar-L_pmu^T Dbar L_pmu; Hbar=H_pmu-L_pmu^T (D_pmu-Dbar) L_pmu',
                'full_integrated_process_and_directional_service_retained': True,
                'zero_reserve_inverse_justified_on_fixed_slice': True,
                'threshold_one_stationary_mode_invariant': True,
                'uniform_eigenmode_exclusion_verified': False,
                'weight_reduction_alone_implies_positive_margin': False,
                'uniform_receipt_margin': None,
                'odd_AW_blocker_resolved': False,
                'complete_homogeneous_gap_proved': False,
            },
            'full_residual_process_receipt_payment': {
                'proof': 'CR44--CR45',
                'scope': 'CR12 actual fixed-target partial fibre, preceding full integrated prediction allocated once as in CR25; nonzero and scalar zero-gap receipts retain their actual guards',
                'actual_operands': "X=F*P*F^T; Qr=Q_ship-q_aw*e_a*e_a^T; v=X_aa; d=e_a^T*X*Qr^-1*X*e_a",
                'paid_receipt_coefficient': 'alpha=(4*d+v)/(v*(2*d+v)^2)>0',
                'signed_substitution': 'L_rem=alpha*mu^2+retained_transverse_loss+retained_full_process_loss; replace L_rem once in the SAME Hbar',
                'full_process_remaining_lower': 'Qr-rank_one_allocation >= Qr/2 > 0',
                'full_integrated_Q_ao_and_all_prefix_deletions_retained': True,
                'nonzero_receipt_lower_inequality_verified': True,
                'uniform_paid_coefficient': None,
                'receipt_threshold_crossed': False,
                'full_endogenous_Schur_cost_paid': False,
                'complete_homogeneous_gap_proved': False,
            },
            'full_endogenous_receipt_schur_reduction': {
                'proof': 'CR41--CR43',
                'scope': 'justified fixed linear part of the existing inherited endogenous lift; scalar cones retain actual image and prefix guard feasibility',
                'positive_pivot': 'K_c=lambda*Hbar-c*J_r_OO; K_0>0 on each injective fixed-word odd slice, common c>0 not qualified',
                'dependent_receipts': 'L=E*L0; L0 consists of independent actual receipt rows',
                'receipt_corner': 'S_c=(L0*K_c^-1*L0^T)^-1-lambda*E^T*Dbar*E',
                'full_remainder': 'R_W(c)=[[S_c,Y_c^T*X_c],[X_c^T*Y_c,A_c-X_c^T*K_c_sharp*X_c]]',
                'full_schur_when_receipt_corner_positive': 'A_c-X_c^T*N_c^-1*X_c=T_c-B_c^T*S_c^-1*B_c; N_c^-1=K_c_sharp+Y_c*S_c^-1*Y_c^T',
                'all_metric_generated_source_reset_inherited_gauge_cross_work_retained': True,
                'restricted_kernel_margin_used_to_price_cross_work': False,
                'odd_margin_assumed_to_form_remainder': False,
                'uniform_positive_reserve_inverse_qualified': False,
                'actual_feasible_receipt_corner_sign_decided': False,
                'full_remainder_uniformly_nonnegative': False,
                'odd_AW_blocker_resolved': False,
                'complete_homogeneous_gap_proved': False,
            },
            'actual_acc_floor_payment': {
                'proof': 'CR18--CR24',
                'scope': 'CR12 regular planar pure covariance fibre; actual active floor then applied acc in the same cycle, with only qualified non-AW S/block resets between',
                'same_prefix_boundary_constraint': 'dC_acc_prefix=-d_beta_face; C_acc_prefix=Cplus_face',
                'same_acc_Fisher_completion': 'L_P=utilde^T K_u utilde+eta*(v/C)^2; v/C=t-b',
                'paid_fraction': 'eta=2*C*(2*s-r-C)/(s*(2*s-r)); 1>eta>C/s>0',
                'remaining_weight': 'omega=1-eta=[C*r+(R+r)*(2*R+r)]/[s*(2*s-r)]>0',
                'linked_remaining_reader': 'b+(eta/omega)*t; b=d_beta_face/Cplus, t=dT_acc B_acc (h_o+T_acc)^T/Cplus',
                'complete_same_storage_matrix': 'G_OO/lambda=H_acc-R_acc^T Omega R_acc; paired acc loss removed once',
                'relative_reader_threshold': 'Omega^-1-R_acc*(H_acc-c_O*H_0)^-1*R_acc^T>=0; H_acc-c_O*H_0>0',
                'uniform_paid_fraction': None,
                'uniform_relative_reader_margin': None,
                'conditional_B_covariance_reader_retained': True,
                'unpaired_faces_and_zero_gap_directional_maps_retained': True,
                'kernel_signed_reserve_lower': '1/10000000000000000000000000000000000000',
                'kernel_scope': 'actual CR12 image, nuisance-supported root, ker R_acc and zero-gap compatibility; first prediction included after 17 s qualified default regular A21',
                'actual_kernel_nontriviality_verified': False,
                'odd_AW_blocker_resolved': False,
            },
            'released_coordinate_Schur_relative_charge': '(2alpha-I)/(2alpha-I-2beta)>1 for beta>0',
            'released_coordinate_failure_class': 'D_SUFFICIENT_BOUND_FAILURE',
            'released_direction_proved_shipping_reachable': False,
            'uniform_actual_linked_deficit_reader_margin': None,
            'uniform_full_homogeneous_gap': None,
            'full_word_OPEN_dependencies_discharged': [],
            'new_storage_constructed': False,
        },
        'planar_even_process_work_proof': 'app:planar-even-process-work',
        'planar_fixed_input_even_dF_dQ_zero': True,
        'planar_fixed_input_even_ds_identically_zero': False,
        'planar_quaternion_polynomial_derivative_defect_retained': True,
        'planar_even_auxiliary_charge_BG_coefficient_upper_real': '1/1000000000000000000000000000',
        'planar_even_auxiliary_charge_scope': 'h<=.006, |h omega_y|<.01, q_g>=1e-6, fixed source/private state; original coordinates before separate projection/AW/reset operations; not a bound for full word work',
        'AW_shear_floor_face_uses_original_recovered_marginal': True,
        'AW_shear_preserves_qualified_planar_magnetic_loss': True,
        'planar_acc_physical_mismatch': 'm=J_y[t daw+omega e_aw-omega t D(t)(a_phys-g)-omega nu], D(t)=(R_y(-t)-I)/t',
        'planar_acc_curvature_charge': '||m0||_(R^-1)^2 <= (e^T P^-1 e) tr(R^-1 J_y B_v P B_v^T J_y^T)',
        'planar_acc_charge_retains_pitch_AW_cross_covariance': True,
        'planar_acc_charge_requires_nominal_AW_or_BA_box': False,
        'nominal_aw_cap_needed_for_dH_coefficient': False,
        'world_frame_is_full_filter_contraction': False,
        'all_time_nominal_innovation_bound_verified': False,
        'uniform_complete_word_gap': None,
        'planar_pitch_domain': planar_pitch_domain_bounds(),
        'planar_magnetic_loss': planar_magnetic_loss_margin(F(3, 500), 6, 100),
        'S_correction_fixed_noise_exact_loss': True,
        'remaining_feedback': ['linked acc mismatch and applied-increment AW-frame jump', 'reference and noise variation',
            'nominal-frame shear connection', 'prediction and reset endpoint-frame derivatives',
            'held bias/source ports', 'AW faces and projection', 'source/gauge/precision qualifications'],
        'same_record_physical_variation': 'dR_nominal=dP=0, hence omega=dPbar=0; not a nominal-root tangent',
        'default_AW_map_unchanged_in_frame': True,
        'planar_12_plus_9_parity_preserved': True,
        'structures_preserved': ['all 21 states and actual covariance/gain/Joseph', 'literal held BA mask',
            'actual nominal rotation and world AW state', 'reference/noise and covariance correlations',
            'injection/reset and inherited frame/clock chronology'],
        'relaxations_introduced': ['regular real-operation derivatives; finite precision and branch crossings separate',
            'CoG profile stated; nonzero lever-arm Jacobian needs its own retained port'],
        'all_time_planar_magnetic_service_verified': False,
        'theorem_closed': False,
    }


def planar_pitch_domain_bounds():
    """Uniform covariance bound on a regular planar service segment.

    Qualified real-operation implication, not planar magnetic admission.
    Two prior actually applied events are supplied by existing service/liveness
    on IMU-rooted windows strictly preceding the correction being analyzed. No initial P upper bound or reset is used.
    Exact source settings, parity and reference cone are compulsory.
    """
    import struct
    def represented(x):
        return F(struct.unpack('f', struct.pack('f', x))[0])
    u, a = F(1, 2**24), F(11, 5000)
    bmin = 75*(1-u)*(1-a*a/2)
    bmax = 75*(1+u)
    r = represented(.8)**2/bmin**2
    qg, qb = represented(.00135)**2, represented(1e-10)
    floor = represented(1e-12)
    # Pick the last IMU roots <=t-3-hmax and <=t-1-hmax. Each lags
    # its target by <=hmax. All events used precede t; no future correction.
    hmax = F(3, 500)
    delta_min, delta_max, tail = 1-hmax, 3+hmax, 1+2*hmax
    ratio = tail/delta_min
    # The exact positive pitch/BG block is unchanged by the real PSD
    # spectral repair. Charging an extra eps I per step is a safe outer;
    # it is NOT a bound on floating-point factorization error.
    floor_charge = 3*floor*(F(4)/F('.004')+2)
    pitch = ((1+ratio)**2+ratio**2)*r + qg*(tail**2/delta_min+tail)
    pitch += qb*(tail**2*delta_max+tail**3)/3+floor_charge
    bias = 2*r/delta_min**2+qg/delta_min+qb*(delta_max/3+tail)+floor_charge
    assert pitch < F(3, 5000)
    assert bias < F(1, 4000)
    reference_angle = (1+u)*a/((1-u)*(1-a*a/2))  # angle <= tan(angle)
    tilt = F(11, 105)  # 6 degrees = pi/30 < 11/105
    relative_min = (1-(tilt+reference_angle)**2/2)/(1+u)-75*u/bmin
    relative_max = 75*(1+u)/bmin
    kappa = max(1-relative_min, relative_max-1)
    tmax = bmax*bmax*F(3, 5000)/represented(.8)**2
    assert kappa < F(3, 500)
    assert tmax < 6
    mean_row_floor = bmin*bmin/(7*represented(.8)**2)
    covariance_row_floor = F(8, 49)*bmin*bmin/represented(.8)**2
    assert mean_row_floor > 1250
    assert covariance_row_floor > 1400
    return {
        'pitch_variance_upper_exact': str(pitch),
        'event_separation_interval': [str(delta_min), str(delta_max)],
        'last_event_age_upper': str(tail),
        'reader_entry_time_seconds': str(3+2*hmax),
        'BG_process_density': str(qb),
        'BG_density_source': 'ProxyStartupFusionConfig::b0=1e-10f, forwarded by initialize_ext; core standalone default differs',
        'pitch_variance_strict_upper': '3/5000',
        'BG_y_variance_upper_exact': str(bias),
        'BG_y_variance_strict_upper': '1/4000',
        'reference_norm_lower': str(bmin),
        'effective_pitch_measurement_variance_upper': str(r),
        'process_floor_reader_charge': str(floor_charge),
        'radial_fraction_abs_upper_exact': str(kappa),
        'radial_fraction_abs_strict_upper': '3/500',
        'magnetic_signal_to_noise_upper_exact': str(tmax),
        'magnetic_signal_to_noise_strict_upper': '6',
        'mean_pitch_loss_coefficient_lower_exact': str(mean_row_floor),
        'mean_pitch_loss_coefficient_strict_lower': '1250',
        'covariance_pitch_row_loss_coefficient_lower_exact': str(covariance_row_floor),
        'covariance_pitch_row_loss_coefficient_strict_lower': '1400',
        'coupled_magnetic_pitch_loss_lower': '1125*eta_pitch^2+1260*(dP J dP)_pitch,pitch',
        'coercivity_scope': 'planar pitch row only; not full transverse storage',
        'covariance_requires_nominal_tilt_bound': False,
        'radial_bound_requires_retained_tilt_chart': True,
        'covariance_scope': 'after 3+2*hmax seconds of regular planar actually-applied service with qualified reference cone; no hard reset in reader window',
        'tilt_scope': 'existing 6-degree retained nominal/central-physical comparison chart',
        'reference_cone_scope': 'qualified cone a=11/5000; float reference accumulation not certified',
        'service_assumed_not_proved_for_planar_witness': True,
        'initial_covariance_upper_assumed': False,
        'every_prefix_retained_chart_verified': False,
    }


def planar_magnetic_loss_margin(kappa_cap, signal_to_noise_cap, comparison_energy_cap,
                                covariance_weight=F(1), retained_fraction=F(9, 10)):
    """Rigorous uniform coupled loss fraction on the stated planar cell.

    Not a point Jacobian or empirical factor. The proof uses exact rank-one
    magnetic structure and the SAME correction's Fisher covariance loss.
    It prices both the cross term and the square of covariance-induced work.
    """
    k, t, E, weight, c = map(F, (kappa_cap, signal_to_noise_cap,
                                comparison_energy_cap, covariance_weight, retained_fraction))
    if min(k, t, E) < 0 or weight <= 0 or not 0 < c < 1:
        raise ValueError('valid nonnegative domain and positive storage weight required')
    ratio = (1+t)/(2+t)
    mean_loss = 1-2*k-k*k*t
    square_charge = k*k*E*t*ratio
    cross_squared = k*k*E*(1+k*t)**2*ratio
    mean_spare = mean_loss-c
    covariance_spare = weight*(1-c)-square_charge
    determinant = mean_spare*covariance_spare-cross_squared
    verified = mean_spare > 0 and covariance_spare > 0 and determinant > 0
    return {'qualification': 'PLANAR_MAGNETIC_LINKED_LOCAL_LOSS',
            'retained_fraction': str(c) if verified else None,
            'mean_loss_lower': str(mean_loss),
            'covariance_square_charge_upper': str(square_charge),
            'cross_coefficient_squared_upper': str(cross_squared),
            'strict_two_by_two_determinant': str(determinant),
            'verified_on_stated_domain': verified,
            'failure_class': None if verified else 'D_SUFFICIENT_BOUND_FAILURE',
            'domain': {'abs_kappa_max': str(k), 'signal_to_noise_max': str(t),
                       'comparison_energy_max': str(E), 'covariance_weight': str(weight)},
            'inequality': 'Wplus<=W-c*(L_eta+lambda*L_cov) on the planar fixed-reference/noise magnetic substep before injection',
            'strict_complete_word_gap_claimed': False,
            'activation_domain_forward_invariant': False}


def rank_one_magnetic_balance(P, h, noise_variance, e, eta, dP, kappa, weight=F(1)):
    """Exact planar magnetic coupled balance in arbitrary active-even dimension.

    h is the literal scalar pitch observation after rotating the 3-D isotropic
    magnetic row. kappa is its linked radial innovation fraction. This tests
    the analytical formula; arbitrary input matrices are not shipping cells.
    """
    ldlt(P)
    noise_variance, kappa, weight = map(F, (noise_variance, kappa, weight))
    if noise_variance <= 0 or weight <= 0 or dP != transpose(dP):
        raise ValueError('positive noise/weight and symmetric covariance tangent required')
    J = inverse(P)
    R = [[noise_variance]]
    S = add(product(h, P, transpose(h)), R)
    K = product(P, transpose(h), inverse(S))
    B = product(K, h)
    A = add(identity(len(P)), B, -1)
    C = product(A, P)
    Jc = inverse(C)
    dC = product(A, dP, transpose(A))
    de = add(eta, product(dP, J, e))
    eta_plus = add(product(A, eta), product(B, de), -kappa)
    Lmean = product(transpose(eta), transpose(h), inverse(S), h, eta)[0][0]
    Lcov = trace(product(J, dP, J, dP))-trace(product(Jc, dC, Jc, dC))
    t = product(h, P, transpose(h))[0][0]/noise_variance
    x, z = product(h, eta)[0][0], product(h, dP, J, e)[0][0]
    a = 1+2*kappa-kappa*kappa*t
    work = -2*kappa*(1-kappa*t)*x*z/S[0][0]+kappa*kappa*t*z*z/S[0][0]
    before = product(transpose(eta), J, eta)[0][0]+weight*trace(product(J, dP, J, dP))
    after = product(transpose(eta_plus), Jc, eta_plus)[0][0]+weight*trace(product(Jc, dC, Jc, dC))
    assert after == before-a*Lmean-weight*Lcov+work
    energy = product(transpose(e), J, e)[0][0]
    assert z*z/S[0][0] <= energy*(1+t)/(2+t)*Lcov
    return {'initial': before, 'final': after, 'mean_loss': Lmean,
            'covariance_loss': Lcov, 'signed_work': work,
            'comparison_energy': energy, 'signal_to_noise': t}
