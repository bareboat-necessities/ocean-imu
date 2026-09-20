"""Matrix endpoint-action certificate for the interleaved A21 LIN covariance.

Enters V_next <= rho V + supply only as a covariance-comparison ingredient.
The raw S-observability normalization is conditional and is NOT a closed-loop
word loss or a full-state rho. Every bound here uses exact rational arithmetic.
"""
from fractions import Fraction as F

from .lin_path_certificate import (
    derivative, hermite_basis, inverse, rational_record, small_x_source_defect,
    sqrt_floor,
)
from .matrix_certificates import encoded, ldlt

SCALES = ('2.4','18','132','4')


def gram(order, time, scales):
    p=[derivative(q,order) for q in hermite_basis()]
    powers=(2,1,0,3)
    return [[scales[i]*scales[j]*time**(powers[i]+powers[j]+1-2*order)*
             sum((a*b/F(k+l+1) for k,a in enumerate(p[i])
                  for l,b in enumerate(p[j])),F(0))
             for j in range(4)] for i in range(4)]


def gram_ceiling(order, scales):
    """G(T) <= G(16)+r I for all T in [16,16.006].

    Each entry is a constant times a single integer power of T, hence monotone.
    The symmetric interval perturbation has norm bounded by its absolute row sum.
    """
    lo=gram(order,F(16),scales); hi=gram(order,F('16.006'),scales)
    r=max(sum(abs(x-y) for x,y in zip(a,b)) for a,b in zip(lo,hi))
    return [[lo[i][j]+(r if i==j else 0) for j in range(4)] for i in range(4)]


def action_matrix(scales=SCALES):
    scales=tuple(F(str(v)) for v in scales)
    if len(scales)!=4 or min(scales)<=0:
        raise ValueError('four positive coordinate scales required')
    g=[gram_ceiling(d,scales) for d in range(5)]
    h=F('.004'); sigma=F('.05'); r_s=F('.075')
    eps,defect,_,_=small_x_source_defect()
    def sampled(d,i,j):
        return 2/h*g[d][i][j]+2*h*g[d+1][i][j]
    # Process energy: (a'+lambda a)^2/(2 sigma^2 lambda).
    # The cross integral is a(T)^2/(2 sigma^2), regardless of lambda switches.
    a=[[2/(1-eps)*((12*g[4][i][j]+50*g[3][i][j]+
                    (scales[3]**2 if i==j==3 else 0))/(2*sigma**2)+
                    defect*sampled(3,i,j))+
        3/sigma**2*sampled(3,i,j)+sampled(0,i,j)/r_s**2
        for j in range(4)] for i in range(4)]
    ldlt(a)
    return a


def _mul(x,y):
    vals=[a*b for a in x for b in y]
    return min(vals),max(vals)


def _sub(x,y):
    return x[0]-y[1],x[1]-y[0]


def _scale(x,s):
    return _mul(x,(s,s))


def _sqrt_upper(q):
    lo=sqrt_floor(q)
    return lo if lo*lo==q else lo+F(1,10**24)


def translation_normalization(a, scales=SCALES):
    """Bound J_S >= mu A_N for raw neutral S rows at three actual event times.

    mu = det(O)_min^2 / (R_S,max * trace(adj(O)' A_N adj(O))_max).
    This is coordinate invariant if A and O receive the same congruence.
    Event times cover [0,.156], [8,8.156], [16,16.156] without sampling.
    """
    v,p,s,_=map(lambda x:F(str(x)),scales)
    an=[row[:3] for row in a[:3]]
    times=[(F(t),F(t)+F('.156')) for t in (0,8,16)]
    ceiling=F(0)
    for j in range(3):
        k,l=(j+1)%3,(j+2)%3
        tk,tl=times[k],times[l]
        diff=_sub(tk,tl)
        square_diff=_sub(_mul(tk,tk),_mul(tl,tl))
        cofactor=[_scale(diff,p*s),_scale(square_diff,-v*s/2),
                  _scale(_mul(_mul(tk,tl),diff),v*p/2)]
        x,y=F(8*k),F(8*l)
        center=[p*s*(x-y),-v*s*(x*x-y*y)/2,v*p*x*y*(x-y)/2]
        delta=[max(abs(z[0]-c),abs(z[1]-c)) for z,c in zip(cofactor,center)]
        e=sum((center[i]*an[i][k]*center[k] for i in range(3) for k in range(3)),F(0))
        r=sum((delta[i]*abs(an[i][k])*delta[k] for i in range(3) for k in range(3)),F(0))
        ceiling+=(_sqrt_upper(e)+_sqrt_upper(r))**2
    det_min=v*p*s/2*F('7.844')**2*F('15.844')
    return det_min**2/(F(100)**2*ceiling)


def certificate():
    a=action_matrix(); lower=inverse(a); _,pivots=ldlt(a)
    mu=translation_normalization(a)
    return {
        'qualification':'OU3_LIN_MATRIX_PATH_ACTION_REAL_ARITHMETIC_V1',
        'verified':True,
        'scope':'A21 regular prediction/correction/reset/PSD-sync words in the bound shipping default profile',
        'state_order':['v','p','S','a_w'], 'scales':list(SCALES),
        'window_s':['16','16.006'], 'tau_s':['0.02','12'],
        'precision_ceiling':encoded(a), 'covariance_lower_matrix':encoded(lower),
        'precision_ldlt_pivots':[str(p) for p in pivots],
        'comparison':'P >= E_LIN A_inverse E_LIN_transpose; a singular full-state lower bound',
        'cross_covariance':'retained through endpoint precision, not discarded',
        'normalization_raw_translation_lower':rational_record(mu),
        'normalization_scope':'raw S rows on neutral endpoint paths with a_w(T)=0; not a complete-word loss',
        'full_state_coercivity_verified':False,
        'constructive_full_A21_mu_rho_enclosure':False,
        'float32_covariance_factor_verified':False, 'theorem_closed':False,
    }
