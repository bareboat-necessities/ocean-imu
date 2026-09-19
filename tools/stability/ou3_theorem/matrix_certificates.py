"""Small exact-rational matrices used by proof certificates, not the estimator."""
from fractions import Fraction as F


def transpose(a):
    return [list(r) for r in zip(*a)]


def matmul(a, b):
    if len(a[0]) != len(b):
        raise ValueError('incompatible matrix shapes')
    return [[sum((x*y for x,y in zip(r,c)), F(0)) for c in zip(*b)] for r in a]


def add(a, b, weight=F(1)):
    if len(a) != len(b) or any(len(x) != len(y) for x,y in zip(a,b)):
        raise ValueError('incompatible matrix shapes')
    return [[x+weight*y for x,y in zip(r,s)] for r,s in zip(a,b)]


def congruence(a, b):
    """b' a b."""
    return matmul(transpose(b), matmul(a,b))


def identity(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]


def ldlt(a):
    """Exact SPD certificate. Pivots are not scalar eigenvalue floors."""
    n=len(a)
    if not n or any(len(r)!=n for r in a) or a != transpose(a):
        raise ValueError('nonempty symmetric square matrix required')
    l=identity(n); d=[]
    for j in range(n):
        p=a[j][j]-sum((l[j][k]**2*d[k] for k in range(j)),F(0))
        if p<=0:
            raise ValueError('matrix is not positive definite')
        d.append(p)
        for i in range(j+1,n):
            l[i][j]=(a[i][j]-sum((l[i][k]*l[j][k]*d[k] for k in range(j)),F(0)))/p
    return l,d


def encoded(a):
    return [[str(v) for v in row] for row in a]


def is_psd(a):
    """Exact semidefinite elimination, including singular Schur complements."""
    if not a or any(len(r)!=len(a) for r in a) or a != transpose(a):
        return False
    b=[list(r) for r in a]
    while b:
        p=b[0][0]
        if p<0:
            return False
        if p==0:
            if any(b[0][j] for j in range(1,len(b))):
                return False
            b=[r[1:] for r in b[1:]]
        else:
            b=[[b[i][j]-b[i][0]*b[0][j]/p for j in range(1,len(b))]
               for i in range(1,len(b))]
    return True
