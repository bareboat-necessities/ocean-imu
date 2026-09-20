"""Exact factor residual certificates, including singular nuisance spaces.

These prove implications for supplied rational matrices. A floating factor
needs a separately certified error envelope before it certifies a true word.
"""
from fractions import Fraction as F
from .matrix_certificates import add, encoded, identity, ldlt, matmul, transpose
from .lin_path_certificate import sqrt_floor


def _rational(a):
    b=[[F(x) for x in row] for row in a]
    if not b or not b[0] or any(len(r)!=len(b[0]) for r in b):
        raise ValueError('nonempty rectangular matrix required')
    return b


def _pivots(a):
    """Exact independent-column indices, with no numerical rank threshold."""
    b=[r[:] for r in a]; row=0; pivots=[]
    for j in range(len(b[0])):
        k=next((k for k in range(row,len(b)) if b[k][j]),None)
        if k is None:
            continue
        b[row],b[k]=b[k],b[row]
        scale=b[row][j]; b[row]=[x/scale for x in b[row]]
        for k in range(row+1,len(b)):
            scale=b[k][j]
            b[k]=[x-scale*y for x,y in zip(b[k],b[row])]
        pivots.append(j); row+=1
        if row==len(b):
            break
    return pivots


def _solve_spd(a,b):
    l,d=ldlt(a); n=len(a); x=[r[:] for r in b]
    for i in range(n):
        for j in range(len(b[0])):
            x[i][j]-=sum((l[i][k]*x[k][j] for k in range(i)),F(0))
    x=[[v/d[i] for v in row] for i,row in enumerate(x)]
    for i in range(n-1,-1,-1):
        for j in range(len(b[0])):
            x[i][j]-=sum((l[k][i]*x[k][j] for k in range(i+1,n)),F(0))
    return x


def eliminate_nuisance(factor,heading_columns):
    """min_n ||B_h h+B_n n||^2 = ||residual h||^2 exactly.

    Select an exact basis of range(B_n), solve only its SPD Gram system,
    and verify orthogonality to ALL nuisance columns. This is the generalized
    Schur complement even when B_n' B_n is singular. The returned minimizer
    sets unused nuisance coefficients to zero; no singular inverse is formed.
    """
    a=_rational(factor); heading=list(heading_columns); width=len(a[0])
    if not heading or len(set(heading))!=len(heading) or any(i<0 or i>=width for i in heading):
        raise ValueError('distinct valid heading columns required')
    nuisance=[i for i in range(width) if i not in heading]
    if not nuisance:
        raise ValueError('at least one nuisance column required')
    h=[[row[i] for i in heading] for row in a]
    n=[[row[i] for i in nuisance] for row in a]
    pivots=_pivots(n)
    minimizer=[[F(0) for _ in heading] for _ in nuisance]
    if pivots:
        basis=[[r[i] for i in pivots] for r in n]; bt=transpose(basis)
        x=_solve_spd(matmul(bt,basis),matmul(bt,h))
        residual=add(h,matmul(basis,x),F(-1))
        for j,i in enumerate(pivots):
            minimizer[i]=[-v for v in x[j]]
    else:
        residual=h
    if any(any(row) for row in matmul(transpose(n),residual)):
        raise ArithmeticError('nuisance range projection failed')
    return {'verified':True,'nuisance_rank':len(pivots),
            'nuisance_gram_singular':len(pivots)<len(nuisance),
            'residual_factor':residual,'reduced_information':matmul(transpose(residual),residual),
            'minimizer_map':minimizer,'full_factor_rank':len(_pivots(a)),
            'source_uniform_verified':False}


def _sqrt_upper(x):
    lo=sqrt_floor(x)
    return lo if lo*lo==x else lo+F(1,10**24)


def singular_floor(factor,right_inverse_candidate,*,factor_error_norm_upper):
    """Residual lower singular bound for a square rational factor.

    If ||B X-I|| <= e<1 then sigma_min(B)>=(1-e)/||X||.
    Subtract the supplied operator error eta to cover ||B_true-B||<=eta.
    Frobenius upper bounds and directed rational roots make this a verified
    matrix implication. The caller must prove eta for the ENTIRE word.
    """
    b,x=_rational(factor),_rational(right_inverse_candidate)
    n=len(b)
    if len(b[0])!=n or len(x)!=n or len(x[0])!=n:
        raise ValueError('square same-size factor and preconditioner required')
    eta=F(factor_error_norm_upper)
    if eta<0:
        raise ValueError('nonnegative factor error required')
    residual=add(matmul(b,x),identity(n),F(-1))
    e=_sqrt_upper(sum((v*v for row in residual for v in row),F(0)))
    xn=_sqrt_upper(sum((v*v for row in x for v in row),F(0)))
    floor=max(F(0),(1-e)/xn-eta) if xn and e<1 else F(0)
    return {'matrix_implication_verified':True,'positive_floor':floor>0,
            'sigma_lower':str(floor),'delta_lower':str(floor*floor),
            'residual_norm_upper':str(e),'factor_error_norm_upper':str(eta),
            'whole_word_error_verified':False,'source_uniform_verified':False}


def audit_certificate():
    cancel=eliminate_nuisance([[1,0,1,1],[0,1,0,0]],[0,1])
    separated=eliminate_nuisance([[1,0,1,1],[0,1,0,0],[0,0,1,1]],[0,1])
    return {'qualification':'OU3_EXACT_SINGULAR_NUISANCE_FACTOR_AUDIT_V1',
            'verified':True,'singular_nuisance_rank':cancel['nuisance_rank'],
            'cancelled_heading_information':encoded(cancel['reduced_information']),
            'separated_heading_information':encoded(separated['reduced_information']),
            'separated_full_rank':separated['full_factor_rank'],
            'full_state_dimension':4,'source_uniform_verified':False,
            'theorem_closed':False}
