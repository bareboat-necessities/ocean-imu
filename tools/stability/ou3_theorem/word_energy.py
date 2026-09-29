"""Exact complete-word covariance-energy algebra and information-lifting audit.

For M mapping the root error to the current error, each prediction contributes
M' (P^-1 - F' Ppred^-1 F) M and each correction contributes
Mprior' H' S^-1 H Mprior. A congruent nonsingular reset contributes zero.
Their sum D satisfies M_end' P_end^-1 M_end + D = P_root^-1.
Only a full-state D >= delta P_root^-1 yields rho <= 1-delta. A principal
submatrix of D cannot be embedded as an independent full-state observation.

This is an exact algebra layer. Supplying one finite word does not certify the
shipping source envelope, nonlinear error, or target arithmetic.
"""
from fractions import Fraction as F
from .lin_path_certificate import inverse
from .matrix_certificates import add, congruence, encoded, identity, is_psd, ldlt, matmul, transpose


def _matrix(a):
    return [[F(v) for v in row] for row in a]


def word_identity(root_covariance, events):
    p=_matrix(root_covariance); ldlt(p)
    root_inverse=inverse(p); m=identity(len(p)); loss=[[F(0) for _ in p] for _ in p]
    for event in events:
        before_inverse=inverse(p)
        kind=event['kind']
        if kind=='prediction':
            f,q=_matrix(event['F']),_matrix(event['Q'])
            if not is_psd(q):
                raise ValueError('process increment must be PSD')
            p=add(congruence(p,transpose(f)),q)
            ldlt(p)
            decrement=add(before_inverse,congruence(inverse(p),f),F(-1))
            loss=add(loss,congruence(decrement,m))
            m=matmul(f,m)
        elif kind=='correction':
            h,r=_matrix(event['H']),_matrix(event['R']); ldlt(r)
            s=add(congruence(p,transpose(h)),r)
            k=matmul(matmul(p,transpose(h)),inverse(s))
            a=add(identity(len(p)),matmul(k,h),F(-1))
            decrement=congruence(inverse(s),h)
            loss=add(loss,congruence(decrement,m))
            # Literal Joseph form; no information-form assumption about Q order.
            p=add(congruence(p,transpose(a)),congruence(r,transpose(k)))
            m=matmul(a,m)
        elif kind=='reset':
            g=_matrix(event['G']); inverse(g)
            p=congruence(p,transpose(g)); m=matmul(g,m)
        else:
            raise ValueError('unknown energy event')
        if not is_psd(loss):
            raise ArithmeticError('energy loss is not PSD')
    endpoint=congruence(inverse(p),m)
    if add(endpoint,loss)!=root_inverse:
        raise ArithmeticError('complete-word energy identity failed')
    return {'algebra_verified':True,'root_precision':root_inverse,
            'endpoint_energy':endpoint,'loss':loss,'end_transport':m,
            'end_covariance':p,'source_uniform_verified':False}


def full_loss_margin(word, delta):
    """Exact finite-word check; never a source-uniform certificate."""
    delta=F(delta)
    if not 0<delta<1:
        raise ValueError('strict loss margin must lie in (0,1)')
    return is_psd(add(word['loss'],word['root_precision'],-delta))



def block_schur_loss_certificate(loss, split):
    """Exact nuisance-eliminated loss and completed-square congruence.

    Partition a symmetric loss D as [[A,B],[B',N]].  If N>0, then
      D = T' diag(S,N) T,  S=A-B N^-1 B',
      T=[[I,0],[N^-1 B',I]].
    Thus S is exactly the AG loss remaining after optimal nuisance
    cancellation.  This is an algebraic reduction only: contraction still
    requires comparison of the complete D with root precision.
    """
    d=_matrix(loss)
    n=len(d)
    if not 0<split<n or d!=transpose(d):
        raise ValueError('symmetric loss and interior split required')
    a=[row[:split] for row in d[:split]]
    b=[row[split:] for row in d[:split]]
    nn=[row[split:] for row in d[split:]]
    ldlt(nn)
    ni=inverse(nn)
    schur=add(a,matmul(matmul(b,ni),transpose(b)),F(-1))
    ldlt(schur)
    t=identity(n)
    nib=matmul(ni,transpose(b))
    for i,row in enumerate(nib):
        for j,v in enumerate(row):
            t[split+i][j]=v
    diag=[[F(0) for _ in range(n)] for _ in range(n)]
    for i,row in enumerate(schur):
        diag[i][:split]=row
    for i,row in enumerate(nn):
        diag[split+i][split:]=row
    if congruence(diag,t)!=d:
        raise ArithmeticError('block completed-square identity failed')
    return {'verified':True,'schur':schur,'nuisance':nn,
            'ag_after_optimal_nuisance_cancellation_spd':True,
            'full_contraction_implied_without_root_comparison':False}


def block_generalized_loss_implication(root_precision, loss, split, delta):
    """Exact sufficient test D >= delta J_root after block-Schur reduction.

    The Schur factorization is used to retain cross-coordinate cancellation;
    the final comparison remains the required full generalized loss inequality.
    """
    cert=block_schur_loss_certificate(loss,split)
    j=_matrix(root_precision); ldlt(j)
    delta=F(delta)
    if not 0<delta<1:
        raise ValueError('delta in (0,1) required')
    ok=is_psd(add(_matrix(loss),j,-delta))
    return {**cert,'delta':str(delta),'full_generalized_loss_margin_verified':ok,
            'rho0_upper':str(1-delta) if ok else None}


def block_schur_self_test():
    """Exact correlated example: principal AG information overstates usable loss."""
    d=[[F(5),F(2),F(2)],[F(2),F(4),F(1)],[F(2),F(1),F(2)]]
    cert=block_schur_loss_certificate(d,2)
    # N=2 and B=(2,1)' so S=A-BB'/2.
    expected=[[F(3),F(1)],[F(1),F(7,2)]]
    if cert['schur']!=expected:
        raise ArithmeticError('unexpected Schur complement')
    return {'verified':True,'principal_AG_loss':encoded([r[:2] for r in d[:2]]),
            'nuisance_eliminated_AG_loss':encoded(expected),
            'cross_cancellation_retained':True}



def relative_schur_contraction_certificate(root_precision, loss, split, delta):
    """Exact block-Schur characterization of D >= delta J in relative metric.

    Partition both D and J.  D-delta J is positive definite iff its nuisance
    block N-delta J_nn is positive definite and the corresponding AG Schur
    complement
      A-delta J_aa
      -(B-delta J_an)(N-delta J_nn)^-1(B-delta J_an)'
    is positive definite.
    This is the source-uniform target that preserves root cross precision;
    no square root, scalar precision cap, or independent principal-block
    information is introduced.
    """
    d,j=_matrix(loss),_matrix(root_precision)
    n=len(d); delta=F(delta)
    if len(j)!=n or not 0<split<n or not 0<delta<1:
        raise ValueError('matching matrices, interior split, delta in (0,1) required')
    e=add(d,j,-delta)
    aa=[row[:split] for row in e[:split]]
    an=[row[split:] for row in e[:split]]
    nn=[row[split:] for row in e[split:]]
    try:
        ldlt(nn)
    except Exception:
        return {'verified':False,'nuisance_relative_positive':False,
                'ag_relative_schur_positive':False}
    schur=add(aa,matmul(matmul(an,inverse(nn)),transpose(an)),F(-1))
    try:
        ldlt(schur); ok=True
    except Exception:
        ok=False
    return {'verified':ok,'nuisance_relative_positive':True,
            'ag_relative_schur_positive':ok,
            'relative_AG_schur':schur if ok else None,
            'equivalent_full_margin':is_psd(e)}



def prediction_relative_identity(root_covariance, transition, process_covariance):
    """Exact relative form of one prediction decrement.

    For P>0, nonsingular F and Q>=0, let P-=F P F'+Q and
      Dp=P^-1-F' P-^-1 F.
    Put C=F^-1 Q F^-T. Then P-=F(P+C)F' and therefore
      Dp=P^-1-(P+C)^-1
        =P^-1 C (P+C)^-1.
    More importantly, after root-energy congruence,
      P^(1/2) Dp P^(1/2)=I-(I+P^-1/2 C P^-1/2)^-1.
    Thus a useful relative prediction margin needs a LOWER bound on C
    relative to P, equivalently C>=epsilon P. A positive absolute Q floor
    alone cannot give source-uniform relative contraction when P is unbounded.
    """
    p,f,q=_matrix(root_covariance),_matrix(transition),_matrix(process_covariance)
    ldlt(p); fi=inverse(f)
    c=matmul(matmul(fi,q),transpose(fi))
    pplus=add(p,c)
    decrement=add(inverse(p),inverse(pplus),F(-1))
    # Original-coordinate prediction gives the same decrement.
    pred=add(congruence(p,transpose(f)),q)
    original=add(inverse(p),congruence(inverse(pred),f),F(-1))
    if decrement!=original:
        raise ArithmeticError('prediction relative identity failed')
    return {'verified':True,'pulled_back_process':c,'decrement':decrement,
            'relative_margin_requires':'F^-1 Q F^-T >= epsilon P_root',
            'absolute_process_floor_sufficient_without_P_upper':False}


def prediction_relative_margin(root_covariance, transition, process_covariance, epsilon):
    """Exact finite-word implication C>=epsilon P => Dp>=eps/(1+eps) P^-1."""
    p,f,q=_matrix(root_covariance),_matrix(transition),_matrix(process_covariance)
    eps=F(epsilon)
    if eps<=0: raise ValueError('positive epsilon required')
    c=matmul(matmul(inverse(f),q),transpose(inverse(f)))
    premise=is_psd(add(c,p,-eps))
    d=prediction_relative_identity(p,f,q)['decrement']
    delta=eps/(1+eps)
    conclusion=is_psd(add(d,inverse(p),-delta))
    if premise and not conclusion:
        raise ArithmeticError('relative prediction implication failed')
    return {'premise_verified':premise,'delta':str(delta),
            'conclusion_verified':conclusion}


def word_smoother_identity(root_covariance, events):
    """Exact fixed-point-smoother form of the complete-word contraction.

    For the optimal-gain frozen word, with C=Cov(x_0,x_k|y) propagated by
    C<-C F', C<-C (I-KH)', C<-C G' and Sigma_00 reduced by C H' S^-1 H C',
      M = C_end' P_0^-1   (exactly),
      M' P_end^-1 M = P_0^-1 (Sigma_00|y - Sigma_00|y,x_end) P_0^-1,
    with Sigma_00|y,x_end = Sigma_00|y - C_end P_end^-1 C_end'. Hence
    rho=lambda_max(P_0^-1/2 (Sigma_00|y - Sigma_00|y,x_end) P_0^-1/2).
    Sigma_00|y is the root covariance after the word's data (information
    mechanism); Sigma_00|y,x_end is what the terminal state leaves unknown
    (forgetting mechanism). Neither alone controls a mixed word.
    """
    p0 = _matrix(root_covariance)
    ldlt(p0)
    p, c, s00, m = p0, p0, p0, identity(len(p0))
    for event in events:
        kind = event['kind']
        if kind == 'prediction':
            f, q = _matrix(event['F']), _matrix(event['Q'])
            if not is_psd(q):
                raise ValueError('process increment must be PSD')
            p, c, m = add(congruence(p, transpose(f)), q), matmul(c, transpose(f)), matmul(f, m)
        elif kind == 'correction':
            h, r = _matrix(event['H']), _matrix(event['R'])
            ldlt(r)
            s = add(congruence(p, transpose(h)), r)
            si = inverse(s)
            k = matmul(matmul(p, transpose(h)), si)
            a = add(identity(len(p)), matmul(k, h), F(-1))
            ch = matmul(c, transpose(h))
            s00 = add(s00, matmul(matmul(ch, si), transpose(ch)), F(-1))
            c = matmul(c, transpose(a))
            p = add(congruence(p, transpose(a)), congruence(r, transpose(k)))
            m = matmul(a, m)
        elif kind == 'reset':
            g = _matrix(event['G'])
            inverse(g)
            p, c, m = congruence(p, transpose(g)), matmul(c, transpose(g)), matmul(g, m)
        else:
            raise ValueError('unknown energy event')
    p0i = inverse(p0)
    if m != matmul(transpose(c), p0i):
        raise ArithmeticError('closed-loop map is not the smoother cross covariance')
    forgetting = add(s00, matmul(matmul(c, inverse(p)), transpose(c)), F(-1))
    explained = add(s00, forgetting, F(-1))
    if congruence(inverse(p), m) != congruence(explained, p0i):
        raise ArithmeticError('smoother contraction identity failed')
    return {'algebra_verified': True, 'root_covariance': p0,
            'smoothed_root_covariance': s00, 'terminal_conditioned_root_covariance': forgetting,
            'explained_root_covariance': explained, 'end_covariance': p,
            'source_uniform_verified': False}


def smoother_margin_tests(word, delta):
    """Exact checks of the full and the two separated sufficient conditions.

    full:        Sigma_00|y - Sigma_00|y,x_end <= (1-delta) P_0  (<=> rho<=1-delta)
    information: Sigma_00|y <= (1-delta) P_0
    forgetting:  Sigma_00|y,x_end >= delta P_0
    Either separated test implies the full one; the converse fails.
    """
    delta = F(delta)
    if not 0 < delta < 1:
        raise ValueError('delta in (0,1) required')
    p0 = word['root_covariance']
    full = is_psd(add([[(1-delta)*v for v in row] for row in p0], word['explained_root_covariance'], F(-1)))
    info = is_psd(add([[(1-delta)*v for v in row] for row in p0], word['smoothed_root_covariance'], F(-1)))
    forget = is_psd(add(word['terminal_conditioned_root_covariance'], p0, -delta))
    if (info or forget) and not full:
        raise ArithmeticError('separated sufficient condition failed to imply the full test')
    return {'full': full, 'information_only': info, 'forgetting_only': forget}


def mixed_mechanism_example():
    """Exact two-state word: separated bounds certify nothing, the word gives 1/2.

    State 1 is observed with no process noise (information), state 2 is
    unobserved with process noise (forgetting). Carried 0.32-s shipping words
    show the same split: translation contracts by forgetting and slow bias
    directions by information, so any source-uniform contraction certificate
    must be one joint matrix inequality, not two scalar ones.
    """
    word = word_smoother_identity(identity(2), [
        {'kind': 'prediction', 'F': identity(2), 'Q': [[0, 0], [0, 1]]},
        {'kind': 'correction', 'H': [[1, 0]], 'R': [[1]]}])
    energy = word_identity(identity(2), [
        {'kind': 'prediction', 'F': identity(2), 'Q': [[0, 0], [0, 1]]},
        {'kind': 'correction', 'H': [[1, 0]], 'R': [[1]]}])['endpoint_energy']
    tests = smoother_margin_tests(word, F(1, 2))
    if energy != [[F(1, 2), 0], [0, F(1, 2)]] or tests != {'full': True, 'information_only': False, 'forgetting_only': False}:
        raise ArithmeticError('mixed-mechanism example failed')
    for delta in (F(1, 10**6),):
        weak = smoother_margin_tests(word, delta)
        if weak['information_only'] or weak['forgetting_only']:
            raise ArithmeticError('a separated bound unexpectedly certified the mixed word')
    return {'verified': True, 'rho': '1/2', 'information_only_margin': '0',
            'forgetting_only_margin': '0', 'joint_matrix_inequality_required': True}


def restricted_service_counterexample():
    # A genuine positive-noise Kalman correction with S=4 I. Its restricted
    # heading/bias loss is I_2, but the full loss has a cancellation direction.
    p=[[F(i==j,10) for j in range(3)] for i in range(3)]
    event={'kind':'correction','H':[[2,0,2],[0,2,0]],
           'R':[[F(16,5),0],[0,F(18,5)]]}
    word=word_identity(p,[event])
    d=word['loss']; e=[[F(1)],[F(0)],[F(-1)]]
    return {
        'qualification':'OU3_RESTRICTED_SERVICE_LIFTING_COUNTEREXAMPLE_V1',
        'verified':d==[[F(1),F(0),F(1)],[F(0),F(1),F(0)],[F(1),F(0),F(1)]],
        'full_loss':encoded(d),'restricted_heading_bias_loss':encoded([r[:2] for r in d[:2]]),
        'null_error':['1','0','-1'],
        'root_energy':str(congruence(word['root_precision'],e)[0][0]),
        'endpoint_energy':str(congruence(word['endpoint_energy'],e)[0][0]),
        'full_state_strict_margin':False,
        'scope':'falsifies the principal-block lifting implication, not a shipping marine-history counterexample',
        'required_replacement':'full transported loss with cross blocks and covariance-whitened coercivity',
    }
