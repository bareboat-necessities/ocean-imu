"""Literal CoG measurement feedback in a moving nominal-world frame.

All 21 states and the same P/K/Joseph operands remain. This is an exact change
of frame, including its derivative, not an invariant-EKF replacement. Post-
injection frame transport, reference/noise changes and branch jumps remain.
"""
from __future__ import annotations
from fractions import Fraction as F

from .information_shear_word import zeros, trace
from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, ldlt, transpose
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
