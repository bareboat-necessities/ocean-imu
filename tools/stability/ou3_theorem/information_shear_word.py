"""Exact information-shear word action, without a uniform-gap assertion.

Analytical identities are proved in app:information-shear-word. Rational
operands test identities, not reachable cells. Prefix maps must be derivatives
of one complete causal shipping history, including auxiliary-state feedback.
"""
from __future__ import annotations
from fractions import Fraction as F

from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, transpose
from .planar_complete_word_storage import product


def zeros(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def trace(a):
    return sum((row[i] for i, row in enumerate(a)), F(0))


def symmetric_basis(n):
    """Unscaled symmetric coordinates; metric retains both off-diagonal terms."""
    result = []
    for i in range(n):
        for j in range(i, n):
            a = zeros(n, n)
            a[i][j] = a[j][i] = F(1)
            result.append(a)
    return result


def joint_metric(P, weight):
    """eta' P^-1 eta + weight tr(P^-1 dP P^-1 dP)."""
    ldlt(P)
    weight = F(weight)
    if weight <= 0:
        raise ValueError('positive covariance weight required')
    n, J = len(P), inverse(P)
    basis = symmetric_basis(n)
    out = zeros(n+len(basis), n+len(basis))
    for i in range(n):
        out[i][:n] = J[i]
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            out[n+i][n+j] = weight*trace(product(J, a, J, b))
    return out


def joint_congruence(B):
    """Base tangent map (eta,dP) -> (B eta,B dP B')."""
    n = len(B)
    if any(len(row) != n for row in B):
        raise ValueError('square base map required')
    basis = symmetric_basis(n)
    out = zeros(n+len(basis), n+len(basis))
    for i in range(n):
        out[i][:n] = B[i]
    indices = [(i, j) for i in range(n) for j in range(i, n)]
    for col, a in enumerate(basis):
        image = product(B, a, transpose(B))
        for row, (i, j) in enumerate(indices):
            out[n+row][n+col] = image[i][j]
    return out


def linked_action_balance(base_maps, metrics, prefix_maps):
    """Exact root quadratic identity with state-dependent ports retained.

    v_i=T_i a, p_i=(T_(i+1)-L_i T_i)a. The ONE root a includes inherited
    auxiliary/source coordinates or their justified lift. T_i may be
    rectangular. No Markov closure on mean/P alone is assumed. Zero-action
    kernel is intersection ker(Q_i T_i), Q_i=J_i-L_i'J_(i+1)L_i >= 0.
    Signed port work is NOT an exogenous forcing bound.
    """
    if not base_maps or len(metrics) != len(base_maps)+1 or len(prefix_maps) != len(metrics):
        raise ValueError('nonempty paired maps and boundary metrics required')
    for metric in metrics:
        ldlt(metric)
    r = len(prefix_maps[0][0])
    action, work = zeros(r, r), zeros(r, r)
    local_losses, ports = [], []
    for i, L in enumerate(base_maps):
        T, Tnext, Jnext = prefix_maps[i], prefix_maps[i+1], metrics[i+1]
        Q = add(metrics[i], product(transpose(L), Jnext, L), -1)
        if not is_psd(Q):
            raise ValueError('D_SUFFICIENT_BOUND_FAILURE: selected base action is not PSD')
        LT = product(L, T)
        E = add(Tnext, LT, -1)
        cross = product(transpose(LT), Jnext, E)
        action = add(action, product(transpose(T), Q, T))
        work = add(work, add(add(cross, transpose(cross)), product(transpose(E), Jnext, E)))
        local_losses.append(Q)
        ports.append(E)
    initial = product(transpose(prefix_maps[0]), metrics[0], prefix_maps[0])
    final = product(transpose(prefix_maps[-1]), metrics[-1], prefix_maps[-1])
    gap = add(action, work, -1)
    if gap != add(initial, final, -1):
        raise ArithmeticError('linked word identity failed')
    return {'action': action, 'signed_port_work': work, 'gap': gap,
            'initial': initial, 'terminal': final, 'local_losses': local_losses,
            'root_port_maps': ports}


def precision_projector(J, basis):
    """Fixed-rank basis only. Rank loss fails at the actual Gram inverse."""
    ldlt(J)
    gram = product(transpose(basis), J, basis)
    ldlt(gram)
    return product(basis, inverse(gram), transpose(basis), J)


def projector_differential(J, basis, dJ, dbasis):
    """Exact dPi, retaining precision and physical-basis transport."""
    gram = product(transpose(basis), J, basis)
    ldlt(gram)
    U = inverse(gram)
    dgram = add(add(product(transpose(dbasis), J, basis),
                     product(transpose(basis), dJ, basis)),
                 product(transpose(basis), J, dbasis))
    out = add(product(dbasis, U, transpose(basis), J),
              product(basis, U, transpose(dbasis), J))
    out = add(out, product(basis, U, transpose(basis), dJ))
    return add(out, product(basis, U, dgram, U, transpose(basis), J), -1)


def congruence_shear(P, G, e, mismatch, dP, dG, de, dmismatch):
    """P+=GPG', e+=Ge+t; t is the LITERAL chart/injection discrepancy.

    G is shipping covariance reset, not an asserted exact finite SO(3) map.
    No bound on t or dt is inferred. G must be invertible on this branch.
    """
    C = product(G, P, transpose(G))
    ldlt(P)
    ldlt(C)
    dC = add(add(product(G, dP, transpose(G)), product(dG, P, transpose(G))),
             product(G, P, transpose(dG)))
    ep = add(product(G, e), mismatch)
    dep = add(add(product(G, de), product(dG, e)), dmismatch)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(dep, product(dC, inverse(C), ep), -1)
    linked = add(product(G, eta), product(G, P, transpose(dG), inverse(transpose(G)), inverse(P), e), -1)
    linked = add(add(linked, dmismatch), product(dC, inverse(C), mismatch), -1)
    assert direct == linked
    return {'eta_plus': direct, 'covariance_plus': C, 'covariance_tangent_plus': dC,
            'base': product(G, eta), 'port': add(linked, product(G, eta), -1)}


def covariance_increment_shear(P, increment, e, dP, dincrement, de):
    """Default AW: actual PSD increment, including its face/target derivative."""
    if not is_psd(increment):
        raise ValueError('actual PSD increment required')
    C = add(P, increment)
    ldlt(P)
    ldlt(C)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(de, product(add(dP, dincrement), inverse(C), e), -1)
    linked = add(eta, product(dP, add(inverse(P), inverse(C), -1), e))
    linked = add(linked, product(dincrement, inverse(C), e), -1)
    assert linked == direct
    return {'eta_plus': direct, 'base': eta, 'port': add(linked, eta, -1)}


def mean_shift_shear(P, e, shift, dP, de, dshift):
    """Bias projection with unchanged P: no covariance-metric contraction assumed."""
    ldlt(P)
    eta = add(de, product(dP, inverse(P), e), -1)
    direct = add(add(de, dshift), product(dP, inverse(P), add(e, shift)), -1)
    linked = add(add(eta, dshift), product(dP, inverse(P), shift), -1)
    assert direct == linked
    return {'eta_plus': direct, 'port': add(linked, eta, -1)}


def word_score_normal_form(base_maps, covariances, comparisons, gain_scores,
                           covariance_ports, mean_ports, dP0, de0):
    """Exact causal tangent identity, NOT a frozen-factor approximation.

    D+=B D B'+U, de+=B de+B D d+v. On an optimal additive correction
    d=H' S^-1 r; elsewhere d=0. U and v are the actual linked remainder
    derivatives (including all auxiliary feedback), never independent inputs.
    The routine verifies the normal form against direct tangent propagation.
    """
    N, n = len(base_maps), len(dP0)
    if not N or any(len(x) != N for x in (gain_scores, covariance_ports, mean_ports)):
        raise ValueError('one score and two linked ports per operation required')
    if len(covariances) != N+1 or len(comparisons) != N+1:
        raise ValueError('actual boundary covariance and comparison required')
    for P in covariances:
        ldlt(P)
    J = [inverse(P) for P in covariances]
    scores = [zeros(n, 1) for _ in range(N+1)]
    local = []
    for i, B in enumerate(base_maps):
        local.append(add(add(product(J[i], comparisons[i]), gain_scores[i]),
                         product(transpose(B), J[i+1], comparisons[i+1]), -1))
    for i in range(N-1, -1, -1):
        scores[i] = add(local[i], product(transpose(base_maps[i]), scores[i+1]))
    # suffix[i] = B_(N-1)...B_i; it is only the dissipative BASE transport.
    suffix = [identity(n) for _ in range(N+1)]
    for i in range(N-1, -1, -1):
        suffix[i] = product(suffix[i+1], base_maps[i])
    eta0 = add(de0, product(dP0, J[0], comparisons[0]), -1)
    eta = add(product(suffix[0], eta0), product(suffix[0], dP0, scores[0]))
    Dnormal = product(suffix[0], dP0, transpose(suffix[0]))
    retained_ports = []
    for i, (U, v) in enumerate(zip(covariance_ports, mean_ports)):
        g = add(v, product(U, J[i+1], comparisons[i+1]), -1)
        port = add(g, product(U, scores[i+1]))
        retained_ports.append(port)
        eta = add(eta, product(suffix[i+1], port))
        Dnormal = add(Dnormal, product(suffix[i+1], U, transpose(suffix[i+1])))
    D, de = dP0, de0
    for B, d, U, v in zip(base_maps, gain_scores, covariance_ports, mean_ports):
        de = add(add(product(B, de), product(B, D, d)), v)
        D = add(product(B, D, transpose(B)), U)
    direct = add(de, product(D, J[-1], comparisons[-1]), -1)
    if D != Dnormal or direct != eta:
        raise ArithmeticError('complete-word score identity failed')
    return {'local_scores': local, 'suffix_scores': scores, 'base_suffixes': suffix,
            'retained_ports': retained_ports, 'eta_terminal': eta,
            'covariance_tangent_terminal': Dnormal}


def _range_solution(G, q):
    """Exact consistent solution, free coordinates zero; no discarded kernel."""
    n = len(G)
    a = [[F(x) for x in G[i]]+[F(q[i][0])] for i in range(n)]
    row, pivots = 0, []
    for col in range(n):
        pivot = next((i for i in range(row, n) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [x/scale for x in a[row]]
        for i in range(n):
            if i != row:
                scale = a[i][col]
                a[i] = [x-scale*y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
    if any(not any(r[:n]) and r[n] for r in a):
        raise ValueError('score outside action range: retain its kernel forcing')
    z = zeros(n, 1)
    for i, col in enumerate(pivots):
        z[col][0] = a[i][n]
    return z


def score_loss_charge(P, C, M, X, score):
    """For G=J-M' C^-1 M>=0 and G z=score, charge chi=score' z.

    ||M X score||_(C^-1)^2 <= chi/2 * Fisher_loss(X).
    A singular G is allowed, but an out-of-range score fails closed.
    This does not certify a uniform chi or full word coercivity.
    """
    ldlt(P)
    ldlt(C)
    if X != transpose(X):
        raise ValueError('symmetric covariance variation required')
    J, Jnext = inverse(P), inverse(C)
    G = add(J, product(transpose(M), Jnext, M), -1)
    if not is_psd(G):
        raise ValueError('D_SUFFICIENT_BOUND_FAILURE: base mean action is not PSD')
    z = _range_solution(G, score)
    chi = product(transpose(score), z)[0][0]
    DX = product(M, X, transpose(M))
    loss = trace(product(J, X, J, X))-trace(product(Jnext, DX, Jnext, DX))
    image = product(M, X, score)
    charge = product(transpose(image), Jnext, image)[0][0]
    if not (chi >= 0 and loss >= 0 and charge <= chi*loss/2):
        raise ArithmeticError('linked score/Fisher-loss inequality failed')
    return {'mean_action': G, 'score_charge': chi, 'covariance_loss': loss,
            'image_energy': charge, 'upper_bound': chi*loss/2}


def comparison_word_budget(base_maps, covariances, comparisons, corrections):
    """Exact inherited comparison-energy budget for the loss-generated score.

    corrections maps an index to its actual optimal (H,R,r) operands. It must
    be a separate additive substep; actual chart/injection shifts follow it.
    Other steps retain s=e_next-B e and its signed cross-plus-square work.
    No all-time activation, source bound, or observability is inferred.
    """
    N, n = len(base_maps), len(covariances[0])
    if not N or len(covariances) != N+1 or len(comparisons) != N+1:
        raise ValueError('paired nonempty word required')
    if any(i not in range(N) for i in corrections):
        raise ValueError('correction index outside word')
    for P in covariances:
        ldlt(P)
    J = [inverse(P) for P in covariances]
    prefix, loss_score, source_score = identity(n), zeros(n, 1), zeros(n, 1)
    loss, nis, supply = F(0), F(0), F(0)
    for i, B in enumerate(base_maps):
        e, ep, P, C = comparisons[i], comparisons[i+1], covariances[i], covariances[i+1]
        Q = add(J[i], product(transpose(B), J[i+1], B), -1)
        if not is_psd(Q):
            raise ValueError('D_SUFFICIENT_BOUND_FAILURE: comparison base action is not PSD')
        if i in corrections:
            H, R, r = corrections[i]
            ldlt(R)
            S = add(product(H, P, transpose(H)), R)
            K = product(P, transpose(H), inverse(S))
            A = add(identity(n), product(K, H), -1)
            if B != A or C != product(A, P) or ep != add(e, product(K, r)):
                raise ValueError('literal optimal additive operands required; retain masks/chart shifts separately')
            mismatch = add(r, product(H, e))
            nis += product(transpose(r), inverse(S), r)[0][0]
            supply += product(transpose(mismatch), inverse(R), mismatch)[0][0]
        else:
            s = add(ep, product(B, e), -1)
            loss += product(transpose(e), Q, e)[0][0]
            supply += 2*product(transpose(e), transpose(B), J[i+1], s)[0][0]
            supply += product(transpose(s), J[i+1], s)[0][0]
            loss_score = add(loss_score, product(transpose(prefix), Q, e))
            source_score = add(source_score, product(transpose(prefix), transpose(B), J[i+1], s), -1)
        prefix = product(B, prefix)
    initial = product(transpose(comparisons[0]), J[0], comparisons[0])[0][0]
    terminal = product(transpose(comparisons[-1]), J[-1], comparisons[-1])[0][0]
    G = add(J[0], product(transpose(prefix), J[-1], prefix), -1)
    z = _range_solution(G, loss_score)
    chi = product(transpose(loss_score), z)[0][0]
    if terminal-initial != -loss-nis+supply or not 0 <= chi <= loss:
        raise ArithmeticError('linked comparison-energy/score budget failed')
    return {'initial_energy': initial, 'terminal_energy': terminal,
            'nonmeasurement_loss': loss, 'innovation_dissipation': nis,
            'signed_supply': supply, 'loss_generated_score': loss_score,
            'source_chart_score': source_score, 'loss_score_charge': chi,
            'root_plus_supply_budget': initial+supply-nis}


def process_source_score(P, B, Q, e, s):
    """Combined process/source score on the ACTUAL PSD-noise support.

    No OU prior is imposed on physical truth. Q t=s is checked exactly;
    unsupported held/chart forcing must be retained outside this identity.
    All values are real-operation operands from one causal history.
    """
    ldlt(P)
    if Q != transpose(Q) or not is_psd(Q):
        raise ValueError('actual PSD process covariance required')
    try:
        t = _range_solution(Q, s)
    except ValueError as exc:
        raise ValueError('physical defect outside process-noise range; retain explicit unsupported forcing') from exc
    C = add(product(B, P, transpose(B)), Q)
    ldlt(C)
    J, Jnext = inverse(P), inverse(C)
    ep = add(product(B, e), s)
    G = add(J, product(transpose(B), Jnext, B), -1)
    score = add(product(J, e), product(transpose(B), Jnext, ep), -1)
    witness = add(e, product(P, transpose(B), t), -1)
    source_action = product(transpose(s), t)[0][0]
    charge = product(transpose(score), witness)[0][0]
    comparison_loss = product(transpose(e), G, e)[0][0]
    source_score_charge = source_action-product(transpose(s), Jnext, s)[0][0]
    signed_cross = product(transpose(e), transpose(B), Jnext, s)[0][0]
    change_budget = (product(transpose(e), J, e)[0][0]+source_action
                     - product(transpose(ep), Jnext, ep)[0][0])
    if not is_psd(G) or product(G, witness) != score or charge != change_budget or charge < 0:
        raise ArithmeticError('combined process/source score identity failed')
    if (source_score_charge < 0 or signed_cross**2 > comparison_loss*source_score_charge
            or charge != comparison_loss+source_score_charge-2*signed_cross):
        raise ArithmeticError('linked comparison/source score decomposition failed')
    return {'covariance_next': C, 'comparison_next': ep, 'mean_loss': G,
            'combined_score': score, 'range_witness': witness,
            'physical_process_action': source_action, 'combined_score_charge': charge,
            'comparison_loss_energy': comparison_loss,
            'source_only_score_charge': source_score_charge,
            'signed_comparison_source_cross': signed_cross}


def source_qualified_word_score(base_maps, covariances, comparisons, corrections,
                                covered_process_indices):
    """Linked total process score budget; uncovered scores are NEVER dropped.

    covered_process_indices selects actual prediction/PSD-addition substeps.
    Corrections use the reached optimal additive branch; chart shifts and
    unsupported held forcing are separate, retained operations. The inherited
    source budget is an identity, not a uniform numerical cap.
    """
    old = comparison_word_budget(base_maps, covariances, comparisons, corrections)
    N, n = len(base_maps), len(covariances[0])
    covered = set(covered_process_indices)
    if covered.intersection(corrections) or any(i not in range(N) for i in covered):
        raise ValueError('disjoint valid process and correction indices required')
    prefix, score, uncovered_score = identity(n), zeros(n, 1), zeros(n, 1)
    budget = source_action = uncovered_change = F(0)
    comparison_loss = source_score_charge = signed_cross = F(0)
    local_charges = []
    for i, B in enumerate(base_maps):
        P, C, e, ep = covariances[i], covariances[i+1], comparisons[i], comparisons[i+1]
        J, Jnext = inverse(P), inverse(C)
        if i in covered:
            Q = add(C, product(B, P, transpose(B)), -1)
            s = add(ep, product(B, e), -1)
            local = process_source_score(P, B, Q, e, s)
            score = add(score, product(transpose(prefix), local['combined_score']))
            budget += local['combined_score_charge']
            source_action += local['physical_process_action']
            comparison_loss += local['comparison_loss_energy']
            source_score_charge += local['source_only_score_charge']
            signed_cross += local['signed_comparison_source_cross']
            local_charges.append(local['combined_score_charge'])
        elif i not in corrections:
            a = add(product(J, e), product(transpose(B), Jnext, ep), -1)
            uncovered_score = add(uncovered_score, product(transpose(prefix), a))
            uncovered_change += (product(transpose(ep), Jnext, ep)[0][0]
                                  - product(transpose(e), J, e)[0][0])
        prefix = product(B, prefix)
    G = add(inverse(covariances[0]), product(transpose(prefix), inverse(covariances[-1]), prefix), -1)
    witness = _range_solution(G, score)
    charge = product(transpose(score), witness)[0][0]
    # Only correction defect action; old signed_supply also includes process
    # cross work, which is already included in the new combined local charges.
    correction_action = F(0)
    for i, (H, R, r) in corrections.items():
        delta = add(r, product(H, comparisons[i]))
        correction_action += product(transpose(delta), inverse(R), delta)[0][0]
    telescope = (old['initial_energy']-old['terminal_energy']+source_action+
                 correction_action-old['innovation_dissipation']+uncovered_change)
    if budget != telescope or not 0 <= charge <= budget:
        raise ArithmeticError('source-qualified word score budget failed')
    if (budget != comparison_loss+source_score_charge-2*signed_cross
            or signed_cross**2 > comparison_loss*source_score_charge):
        raise ArithmeticError('grouped signed score budget failed')
    if not is_psd(add([[budget*x for x in row] for row in G], product(score, transpose(score)), -1)):
        raise ArithmeticError('directional score Gram bound failed')
    total = add(score, uncovered_score)
    if total != add(old['loss_generated_score'], old['source_chart_score']):
        raise ArithmeticError('uncovered score was lost')
    return {'covered_score': score, 'uncovered_score': uncovered_score,
            'actual_total_score': total, 'mean_action': G,
            'covered_score_charge': charge, 'local_score_charges': local_charges,
            'linked_budget': budget, 'physical_process_action': source_action,
            'correction_defect_action': correction_action,
            'uncovered_signed_energy_change': uncovered_change,
            'comparison_loss_energy': comparison_loss,
            'source_only_score_charge': source_score_charge,
            'signed_comparison_source_cross': signed_cross,
            'uniform_budget_verified': False}


def held_bias_boundary_score(Pb, bias_initial, bias_final, duration,
                             slow_amplitude, slow_rate, precision_ceiling):
    """Exact unsupported score on the reached, entirely held BA stratum.

    Structural hypotheses are proved from shipping in app:held-score-boundary:
    P=diag(P_X,Pb), B=diag(B_X,I), constant held nominal BA and Pb. This
    helper checks the physical endpoint bounds and the ACTUAL Pb comparison.
    It does not remove BA residuals/dR_eff or cover an A21/release word.
    """
    ldlt(Pb)
    duration, slow_amplitude, slow_rate, precision_ceiling = map(
        F, (duration, slow_amplitude, slow_rate, precision_ceiling))
    if min(duration, slow_amplitude, slow_rate) < 0 or precision_ceiling <= 0:
        raise ValueError('nonnegative physical bounds and positive precision ceiling required')
    Jb = inverse(Pb)
    if not is_psd(add([[precision_ceiling*x for x in row] for row in identity(len(Pb))], Jb, -1)):
        raise ValueError('actual held precision exceeds the supplied ceiling')
    for b in (bias_initial, bias_final):
        if product(transpose(b), b)[0][0] > slow_amplitude**2:
            raise ValueError('physical SLOW endpoint exceeds amplitude contract')
    delta = add(bias_final, bias_initial, -1)
    cap_squared = min(slow_rate**2*duration**2, 4*slow_amplitude**2)
    if product(transpose(delta), delta)[0][0] > cap_squared:
        raise ValueError('physical SLOW increment exceeds same-history contract')
    q = product(Jb, delta)
    energy = product(transpose(q), Pb, q)[0][0]
    upper = precision_ceiling*cap_squared
    if energy > upper:
        raise ArithmeticError('held unsupported boundary score cap failed')
    return {'held_score': q, 'held_score_energy': energy,
            'uniform_boundary_upper': upper,
            'active_score_component_identically_zero': True,
            'held_bias_remains_in_residual_and_effective_noise': True,
            'release_or_active_BA_covered': False}


def process_augmented_shear(P, B, Q, e, s, D, de, dB, dQ, ds, weight):
    """Exact same-operation generated process-port balance, including dB/dQ/ds.

    An algebraic augmentation factors the literal derivative; it is NOT an
    independent stochastic-source model or a replacement observer. Q must be
    SPD here. H18 unsupported coordinates and nonsmooth faces stay explicit.
    Derivatives must be the outputs of the same causal auxiliary/source lift.
    """
    ldlt(P)
    ldlt(Q)
    if weight <= 0 or D != transpose(D) or dQ != transpose(dQ):
        raise ValueError('positive weight and symmetric covariance derivatives required')
    n, J, Qinv = len(P), inverse(P), inverse(Q)
    C = add(product(B, P, transpose(B)), Q)
    Jnext = inverse(C)
    ep = add(product(B, e), s)
    DP = add(add(add(product(B, D, transpose(B)), product(dB, P, transpose(B))),
                 product(B, P, transpose(dB))), dQ)
    dep = add(add(product(B, de), product(dB, e)), ds)
    eta = add(de, product(D, J, e), -1)
    etap = add(dep, product(DP, Jnext, ep), -1)
    Sigma, X = zeros(2*n, 2*n), zeros(2*n, 2*n)
    cross = product(P, transpose(dB))
    for i in range(n):
        for j in range(n):
            Sigma[i][j], Sigma[n+i][n+j] = P[i][j], Q[i][j]
            X[i][j], X[n+i][n+j] = D[i][j], dQ[i][j]
            X[i][n+j], X[n+j][i] = cross[i][j], cross[i][j]
    Sinv, channel = inverse(Sigma), [B[i]+identity(n)[i] for i in range(n)]
    z = e+s
    nu = add(ds, product(dQ, Qinv, s), -1)
    aeta = add(eta, product(P, transpose(dB), Qinv, s), -1)+nu
    action = add(Sinv, product(transpose(channel), Jnext, channel), -1)
    score = product(action, z)
    feedback = product(channel, X, score)
    base_eta = product(channel, aeta)
    if etap != add(base_eta, feedback) or DP != product(channel, X, transpose(channel)):
        raise ArithmeticError('literal process derivative augmentation failed')
    auxiliary_work = -2*product(transpose(eta), transpose(dB), Qinv, s)[0][0]
    auxiliary_work += product(transpose(s), Qinv, dB, P, transpose(dB), Qinv, s)[0][0]
    auxiliary_work += product(transpose(nu), Qinv, nu)[0][0]
    auxiliary_work += weight*(2*trace(product(Qinv, dB, P, transpose(dB)))+
                              trace(product(Qinv, dQ, Qinv, dQ)))
    mean_loss = product(transpose(aeta), action, aeta)[0][0]
    Fisher_loss = trace(product(Sinv, X, Sinv, X))-trace(product(Jnext, DP, Jnext, DP))
    mixed = 2*product(transpose(base_eta), Jnext, feedback)[0][0]
    square = product(transpose(feedback), Jnext, feedback)[0][0]
    chi = product(transpose(z), action, z)[0][0]
    initial = product(transpose(eta), J, eta)[0][0]+weight*trace(product(J, D, J, D))
    terminal = product(transpose(etap), Jnext, etap)[0][0]+weight*trace(product(Jnext, DP, Jnext, DP))
    if terminal-initial != auxiliary_work-mean_loss-weight*Fisher_loss+mixed+square:
        raise ArithmeticError('signed augmented process balance failed')
    if min(mean_loss, Fisher_loss, chi, square) < 0 or square > chi*Fisher_loss/2:
        raise ArithmeticError('linked augmented process Fisher charge failed')
    if chi != process_source_score(P, B, Q, e, s)['combined_score_charge']:
        raise ArithmeticError('process/source charge mismatch')
    return {'covariance_tangent_next': DP, 'shear_next': etap,
            'signed_auxiliary_work': auxiliary_work, 'hidden_mean_loss': mean_loss,
            'augmented_Fisher_loss': Fisher_loss, 'signed_feedback_cross': mixed,
            'feedback_square': square, 'combined_process_score_charge': chi,
            'feedback_square_upper': chi*Fisher_loss/2,
            'actual_storage_change': terminal-initial,
            'uniform_absorption_verified': False}


def conditional_mixed_coefficients(C, T, q, weight):
    """Exact conditional root mean/Fisher Schur operator, not a uniform bound.

    C is actual conditional root covariance, T=L' J_N L, q=E_aw' q_word.
    Positive C^-1-T comes from the qualified first process loss, retained
    through subsequent dissipative base maps. Full signed ports are separate.
    Unscaled symmetric coordinates retain BOTH off-diagonal trace terms.
    """
    ldlt(C)
    ldlt(T)
    if weight <= 0:
        raise ValueError('positive covariance weight required')
    Jc = inverse(C)
    A = add(Jc, T, -1)
    ldlt(A)
    K = add(T, product(T, inverse(A), T))
    basis = symmetric_basis(len(C))
    G = [[trace(add(product(Jc, X, Jc, Y), product(T, X, T, Y), -1))
          for Y in basis] for X in basis]
    ldlt(G)
    columns = [product(X, q) for X in basis]
    E = [[col[i][0] for col in columns] for i in range(len(C))]
    reader = product(E, inverse(G), transpose(E))
    margin = add([[weight*x for x in row] for row in inverse(K)], reader, -1)
    operator = add([[weight*x for x in row] for row in G],
                   product(transpose(E), K, E), -1)
    # Two Schur complements of [[lambda G, E'],[E,K^-1]].
    assert is_psd(operator) == is_psd(margin)
    return {'mean_loss': A, 'Fisher_loss_Gram': G, 'score_map': E,
            'completed_score_weight': K, 'directional_score_reader': reader,
            'three_row_margin': margin, 'net_conditional_Fisher_operator': operator}


def conditional_word_mixed_bound(P, PN, M, q, eta0, D0, etaN, DN,
                                 selector, weight, endpoint_adjustment):
    """Sharp linked-packet bound for the full actual suffix-word operands.

    etaN,DN MUST be the complete tangent, including generated U*q suffixes.
    The aggregate packets can depend on ALL root/source/auxiliary coordinates;
    completion is a pointwise identity, not minimization over independent noise.
    For input-frame operands endpoint_adjustment is Phi_0-Phi_N-Gamma_0+Gamma_N;
    for already transformed operands it is only -Gamma_0+Gamma_N. Never count
    frame work twice. Callers retain it explicitly (zero for an unquotiented
    identity). No shipping margin/domain is inferred here.
    """
    ldlt(P)
    ldlt(PN)
    if D0 != transpose(D0) or DN != transpose(DN):
        raise ValueError('symmetric actual covariance tangents required')
    J, JN = inverse(P), inverse(PN)
    C = inverse(product(transpose(selector), J, selector))
    u = product(C, transpose(selector), J, eta0)
    Y = product(C, transpose(selector), J, D0, J, selector, C)
    L = product(M, selector)
    T, qa = product(transpose(L), JN, L), product(transpose(selector), q)
    coeff = conditional_mixed_coefficients(C, T, qa, weight)
    A, S = coeff['mean_loss'], coeff['net_conditional_Fisher_operator']
    try:
        ldlt(S)
    except (ValueError, ArithmeticError, ZeroDivisionError) as exc:
        raise ValueError('D_SUFFICIENT_BOUND_FAILURE: conditional signed Fisher margin is not positive; no kernel deletion') from exc
    t = add(etaN, product(L, add(u, product(Y, qa))), -1)
    Z = add(DN, product(L, Y, transpose(L)), -1)
    b = product(transpose(L), JN, t)
    v = product(transpose(L), JN, Z, JN, L)
    hmean = product(add(identity(len(C)), product(T, inverse(A))), b)
    basis = symmetric_basis(len(C))
    h = [[product(transpose(qa), X, hmean)[0][0]+weight*trace(product(X, v))]
         for X in basis]
    y = [[Y[i][j]] for i in range(len(C)) for j in range(i, len(C))]
    shifted_mean = add(u, product(inverse(A), add(product(T, Y, qa), b)), -1)
    shifted_covariance = add(y, product(inverse(S), h), -1)
    mean_square = product(transpose(shifted_mean), A, shifted_mean)[0][0]
    covariance_square = product(transpose(shifted_covariance), S, shifted_covariance)[0][0]
    initial = product(transpose(eta0), J, eta0)[0][0]+weight*trace(product(J, D0, J, D0))
    terminal = product(transpose(etaN), JN, etaN)[0][0]+weight*trace(product(JN, DN, JN, DN))
    Jc = inverse(C)
    conditional = product(transpose(u), Jc, u)[0][0]+weight*trace(product(Jc, Y, Jc, Y))
    complement = initial-conditional
    packet_self = product(transpose(t), JN, t)[0][0]+weight*trace(product(JN, Z, JN, Z))
    coupling = product(transpose(b), inverse(A), b)[0][0]
    coupling += product(transpose(h), inverse(S), h)[0][0]
    packet = packet_self+coupling
    net_conditional = product(transpose(u), A, u)[0][0]
    net_conditional -= 2*product(transpose(u), T, Y, qa)[0][0]
    net_conditional -= product(transpose(qa), Y, T, Y, qa)[0][0]
    net_conditional += weight*product(transpose(y), coeff['Fisher_loss_Gram'], y)[0][0]
    lower = complement-packet+endpoint_adjustment
    half_lower = net_conditional/2+complement-packet_self-2*coupling+endpoint_adjustment
    gap = initial-terminal+endpoint_adjustment
    assert complement >= 0 and min(mean_square, covariance_square) >= 0
    assert gap == mean_square+covariance_square+lower
    assert net_conditional >= 0 and gap >= half_lower
    return {**coeff, 'conditional_covariance': C, 'conditional_root_mean': u,
            'conditional_root_covariance_tangent': Y, 'aggregate_mean_packet': t,
            'aggregate_covariance_packet': Z, 'linked_packet_charge': packet,
            'packet_self_energy': packet_self, 'linked_mixed_coupling_charge': coupling,
            'net_conditional_root_form': net_conditional, 'half_loss_lower_bound': half_lower,
            'complementary_root_storage': complement,
            'retained_mean_square': mean_square, 'retained_covariance_square': covariance_square,
            'endpoint_adjustment': endpoint_adjustment, 'actual_signed_gap': gap,
            'linked_lower_bound': lower, 'uniform_margin_verified': False}


def bordered_comparison_storage(P, e, kappa, D, de, dkappa):
    """Exact bordered Fisher metric; no uniform equivalence is asserted.

    kappa is proof bookkeeping, never an estimator state. Its Schur slack
    must be positive. Every covariance entry and the scalar tangent remain.
    """
    ldlt(P)
    n, J = len(P), inverse(P)
    M = [P[i][:]+[e[i][0]] for i in range(n)]+[[x[0] for x in e]+[kappa]]
    dM = [D[i][:]+[de[i][0]] for i in range(n)]+[[x[0] for x in de]+[dkappa]]
    V = product(transpose(e), J, e)[0][0]
    c = kappa-V
    if c <= 0:
        raise ValueError('positive bordered comparison slack required')
    ldlt(M)
    eta = add(de, product(D, J, e), -1)
    dc = (dkappa-2*product(transpose(e), J, de)[0][0]
          +product(transpose(e), J, D, J, e)[0][0])
    fisher = trace(product(inverse(M), dM, inverse(M), dM))
    split = trace(product(J, D, J, D))+2*product(transpose(eta), J, eta)[0][0]/c+(dc/c)**2
    if fisher != split:
        raise ArithmeticError('bordered Fisher decomposition failed')
    return {'matrix': M, 'tangent': dM, 'slack': c, 'slack_tangent': dc,
            'Fisher_storage': fisher, 'mean_weight': 2/c,
            'uniform_coercivity_verified': False}


def bordered_process(P, e, kappa, B, Q, s):
    """Literal supported process as a PSD addition on the bordered matrix."""
    n = len(P)
    old = bordered_comparison_storage(P, e, kappa, zeros(n, n), zeros(n, 1), F(0))
    p = process_source_score(P, B, Q, e, s)
    kp = kappa+p['physical_process_action']
    new = bordered_comparison_storage(p['covariance_next'], p['comparison_next'], kp,
                                     zeros(n, n), zeros(n, 1), F(0))
    A = [row[:]+[F(0)] for row in B]+[[F(0)]*n+[F(1)]]
    U = [Q[i][:]+[s[i][0]] for i in range(n)]
    U += [[x[0] for x in s]+[p['physical_process_action']]]
    if not is_psd(U) or new['matrix'] != add(product(A, old['matrix'], transpose(A)), U):
        raise ArithmeticError('bordered process addition failed')
    if new['slack']-old['slack'] != p['combined_score_charge']:
        raise ArithmeticError('bordered slack/score identity failed')
    return {'P_next': p['covariance_next'], 'e_next': p['comparison_next'],
            'kappa_next': kp, 'matrix_next': new['matrix'], 'base': A,
            'addition': U, 'slack_next': new['slack'],
            'slack_increment': p['combined_score_charge']}


def bordered_correction(P, e, kappa, H, R, delta):
    """Optimal additive correction as one joint covariance Schur complement.

    delta=r+H e is the ACTUAL comparison defect, including -S_physical.
    Fixed operands yield Fisher loss; dH,dR,d(delta) require the full lift
    derivative. Arbitrary held masks are not optimal corrections.
    """
    n, m = len(P), len(R)
    ldlt(R)
    old = bordered_comparison_storage(P, e, kappa, zeros(n, n), zeros(n, 1), F(0))
    S = add(product(H, P, transpose(H)), R)
    K = product(P, transpose(H), inverse(S))
    r = add(delta, product(H, e), -1)
    C = add(P, product(K, H, P), -1)
    ep = add(e, product(K, r))
    kp = (kappa+product(transpose(delta), inverse(R), delta)[0][0]
          -product(transpose(r), inverse(S), r)[0][0])
    new = bordered_comparison_storage(C, ep, kp, zeros(n, n), zeros(n, 1), F(0))
    Sigma, L = zeros(n+1+m, n+1+m), identity(n+1+m)
    for i in range(n+1):
        Sigma[i][:n+1] = old['matrix'][i][:]
    for i in range(m):
        Sigma[n+1+i][n+1:] = R[i][:]
        L[n+1+i][:n] = H[i][:]
    t = product(inverse(R), delta)
    L[n][n+1:] = [-row[0] for row in t]
    Y = product(L, Sigma, transpose(L))
    cross = [row[n+1:] for row in Y[:n+1]]
    regression = product(cross, inverse(S))
    schur = add([row[:n+1] for row in Y[:n+1]], product(regression, transpose(cross)), -1)
    if schur != new['matrix'] or new['slack'] != old['slack']:
        raise ArithmeticError('bordered correction Schur/slack identity failed')
    return {'P_next': C, 'e_next': ep, 'kappa_next': kp,
            'matrix_next': new['matrix'], 'slack_next': new['slack'],
            'input': Sigma, 'lift': L, 'joint': Y,
            'innovation': S, 'regression': regression}


def certificate():
    return {
        'qualification': 'OU3_INFORMATION_SHEAR_WORD_V2',
        'result_type': 'PROVED — analytical balance and kernel characterization; uniform coercivity OPEN',
        'proof_appendix_label': 'app:information-shear-word',
        'governing_strategy': 'symmetry/kernel -> justified quotient/tube -> linked dissipation -> forcing -> minimal activation domain',
        'general_theorem_assumes_existing_magnetic_service': True,
        'general_theorem_waits_for_planar_admission': False,
        'storage': 'W=eta^T J eta+lambda tr(J dP J dP), eta=de-dP J e, lambda>0',
        'word_balance': 'W_N-W_0=-sum v_i^T Q_i v_i+sum[2(L_i v_i)^T J_(i+1) p_i+p_i^T J_(i+1) p_i]',
        'port_scope': 'literal linked state-dependent/source ports, not independent disturbances; inherited auxiliary tangents remain causal',
        'covariance_zero_loss': '(I-U^T U) Z=0, Z=P^-1/2 dP P^-1/2',
        'optimal_correction_zero_loss': 'H eta=0 and H dP=0; reached held BA uses active block and retained held ports',
        'word_action_kernel': 'intersection_i ker(Q_i T_i) for ACTUAL complete causal prefix tangents T_i',
        'quotient_balance': 'Wperp_N-Wperp_0=-A_W+F_W+gauge_energy_0-gauge_energy_N',
        'moving_projector_differential_retained': True,
        'physical_compatibility_equals_homogeneous_kernel': False,
        'source_uniform_coercivity_requirement': 'B^T[G_action-G_port-G_gauge,0+G_gauge,N]B >= c B^T Jperp,0 B, c>0, on a justified transverse lift and every admitted activated word',
        'relative_Schur_obligation': 'for G-cJ=[[A,B],[B^T,N]], N>0 and A-B N^-1 B^T>=0 on justified service-coordinate strata; actual cross precision retained',
        'uniform_c': None,
        'word_score_proof': 'app:linked-word-score',
        'word_score_note': 'docs/ou3-linked-word-score.md',
        'optimal_additive_local_score': 'J e+H^T S^-1 r-A^T C^-1(e+K r)=0',
        'word_score_normal_form': 'eta_N=M eta_0+M dP_0 q_0+sum M_(N,i+1)[v_i-U_i J_(i+1)e_(i+1)+U_i q_(i+1)]',
        'generated_covariance_suffix_score_retained': True,
        'score_loss_charge': 'G z=q implies ||M X q||_(J_N)^2 <= (q^T z/2) L_P(X), G=J_0-M^T J_N M>=0',
        'score_outside_action_range': 'retain kernel component as forcing; never discard it by pseudoinverse',
        'uniform_suffix_score_charge': None,
        'conditional_mixed_word_proof': 'app:conditional-mixed-word',
        'conditional_mixed_mean_Fisher_completion': True,
        'conditional_directional_score_test': 'lambda K^-1-E_q G_F^-1 E_q^T>0; actual C,T,q, not scalar chi/c',
        'conditional_generated_packets': 't=eta_N-L(u+Y q_a), Z=dP_N-L Y L^T; complete actual tangents required',
        'conditional_linked_packet_bound': 'gap=mean_square+covariance_square+W_complement-PacketCharge+endpoint_adjustment',
        'conditional_half_loss_bound': 'gap>=net_conditional/2+W_complement-packet_self-2 linked_mixed_charge+endpoint_adjustment',
        'conditional_packet_independence_assumed': False,
        'conditional_AW_state_dependent_frame_invariance': True,
        'conditional_frame_work_retained_in_complementary_packets': True,
        'conditional_packet_absorption_proof': 'app:conditional-packet-absorption',
        'complementary_base_packet_Schur_absorption_verified': True,
        'fixed_half_charge_domination_is_necessary': False,
        'isotropic_score_cancellation_extends_to_anisotropic_loss': False,
        'process_source_score_proof': 'app:process-source-score',
        'supported_process_source_score_identity_verified': True,
        'combined_process_word_directional_budget_verified': True,
        'unsupported_source_score_discarded': False,
        'causal_generated_process_augmented_balance_verified': True,
        'generated_process_square_Fisher_charge_verified': True,
        'generated_process_cross_absorbed': False,
        'combined_process_score_budget_AW_frame_invariant': True,
        'uniform_combined_process_source_budget': None,
        'held_score_boundary_proof': 'app:held-score-boundary',
        'held_BA_unsupported_score_uniform_boundary_bound': True,
        'held_BA_unsupported_score_active_component_zero': True,
        'held_BA_residual_and_effective_noise_ports_retained': True,
        'active_process_budget_minimal_domain': 'B_W=E_L+S_q-2 C; C^2<=E_L S_q. With bounded S_q, uniform B_W iff uniform comparison loss energy E_L; neither cap is supplied for general activated A21 words.',
        'uniform_conditional_score_margin_verified': False,
        'uniform_generated_packet_absorption_verified': False,
        'bordered_comparison_storage_proof': 'app:bordered-comparison-storage',
        'bordered_process_correction_score_absorption': 'PROVED for fixed process/row/noise/comparison-defect operands; causal variations remain signed ports',
        'bordered_slack_balance': 'c_N=c_0+B_W on a covered process/correction word; no physical or estimator restart',
        'bordered_storage_uniform_coercivity_verified': False,
        'bordered_storage_full_endogenous_work_absorbed': False,
        'comparison_word_budget': 'chi_loss <= E_loss = V_0-V_N+Supply_W-sum NIS <= V_0+Supply_W-sum NIS',
        'comparison_supply': 'sum corrections ||r+H e||_(R^-1)^2 plus actual prediction/reset/projection signed cross-plus-square work',
        'physical_S_reset_by_pseudo_measurement': False,
        'realized_NIS_is_observability_lower_bound': False,
        'minimal_obligations_note': 'docs/ou3-information-shear-strategy.md',
        'full_literal_derivative_domain_certified': False,
        'nonlinear_remainder_and_precision_transfer_certified': False,
        'structures_preserved': ['complete 21-state mean and covariance', 'literal Joseph/masks and effective held noise',
            'prediction mismatch and dF/dQ', 'dH/dR/dr', 'injection/reset and bias projection',
            'default AW pending-target/face differential', 'inherited reference/gates/clocks',
            'same-history gauge/source and moving projectors'],
        'relaxations_introduced': ['regular real-operation tangent identities; nonsmooth/arithmetic closure remains separate',
            'dissipative base plus exact signed ports is an identity, not a surrogate execution'],
        'radius_solve_performed': False,
        'all_time_planar_magnetic_service_verified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
