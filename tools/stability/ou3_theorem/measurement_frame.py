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
