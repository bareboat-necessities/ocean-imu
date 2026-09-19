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
