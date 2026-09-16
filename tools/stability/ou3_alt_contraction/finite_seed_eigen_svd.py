"""Pinned Eigen 3.4.0 scalar 2x3 QR/Jacobi seed dependency graph.

The axis and loop predicate are computed from the input vectors. No solver
axis, reflector, or reciprocal can be supplied by a caller. The graph retains
IEEE infinities/NaNs in the Jacobi work matrix: those unused entries need not
be finite for Eigen to return its unchanged third Q column. Status flags and
floating-point traps are outside this named nontrapping scalar profile.

QR norm downdates and singular-value sorting are omitted from the dependency
slice: the second pivot has only one candidate; sorting touches columns 0/1
only. Both loops have statically bounded lengths. The unbounded Jacobi loop
is executed with cycle detection, not silently replaced by a fixed iteration
count. The source-uniform bound is proved separately by the exact rational
normalization/QR/two-sweep induction in finite_seed_svd_roundoff.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as START
from tools.stability.ou3_alt_contraction import finite_seed_svd_axis_reduction as REDUCTION

INF = 0x7f800000
NAN = 0x7fc00000
ZERO, ONE, TWO = 0, 0x3f800000, 0x40000000
MIN_NORMAL = 0x00800000
PRECISION = M.exact_bits(F(1, 2**22))
PROFILE = M.PROFILE + '-nontrapping-ieee-specials'


def special(b):
    return (b & INF) == INF


def isnan(b):
    return special(b) and (b & 0x7fffff) != 0


def zero(b):
    return (b & ~M.SIGN) == 0


def neg(b):
    return b ^ M.SIGN


def absolute(b):
    return b & ~M.SIGN


def less(a,b):
    if isnan(a) or isnan(b): return False
    if a == b or zero(a) and zero(b): return False
    if special(a): return bool(a & M.SIGN)
    if special(b): return not bool(b & M.SIGN)
    return M.value(a) < M.value(b)


def maximum(a,b):
    # Eigen numext::maxi (std::max on this profile) keeps its first operand
    # when the ordered comparison is false, including unordered NaNs.
    return b if less(a,b) else a


@dataclass(frozen=True)
class Op:
    kind: str
    args: tuple[int,...]
    result: int


class Arithmetic:
    def __init__(self): self.operations = []

    def _round(self,x,negative_zero=False):
        try: return M.round_bits(x,negative_zero=negative_zero)
        except OverflowError: return INF | (M.SIGN if x < 0 else 0)

    def op(self,kind,*args):
        if any(isnan(x) for x in args):
            out=NAN
        elif kind=='sqrt':
            x=args[0]
            if zero(x): out=x
            elif x & M.SIGN: out=NAN
            elif special(x): out=INF
            else: out=M.exact_bits(START.sqrt32(M.value(x)))
        else:
            a,b=args
            if kind=='add':
                if special(a) and special(b): out=a if a==b else NAN
                elif special(a): out=a
                elif special(b): out=b
                else: out=self._round(M.value(a)+M.value(b),a==M.SIGN and b==M.SIGN)
            elif kind=='mul':
                sign=(a ^ b) & M.SIGN
                if special(a) or special(b):
                    out=NAN if zero(a) or zero(b) else INF|sign
                else: out=self._round(M.value(a)*M.value(b),bool(sign))
            elif kind=='div':
                sign=(a ^ b) & M.SIGN
                if special(a) and special(b) or zero(a) and zero(b): out=NAN
                elif special(a) or zero(b): out=INF|sign
                elif special(b): out=sign
                else: out=self._round(M.value(a)/M.value(b),bool(sign))
            else: raise ValueError('unknown scalar Eigen operation')
        self.operations.append(Op(kind,tuple(args),out))
        return out

    def add(self,a,b): return self.op('add',a,b)
    def sub(self,a,b): return self.add(a,neg(b))
    def mul(self,a,b): return self.op('mul',a,b)
    def div(self,a,b): return self.op('div',a,b)
    def sqrt(self,a): return self.op('sqrt',a)

    def sum(self,values,*,fixed=False):
        values=tuple(values)
        if not values: return ZERO
        if len(values)==1: return values[0]
        if fixed:
            k=len(values)//2
            return self.add(self.sum(values[:k],fixed=True),self.sum(values[k:],fixed=True))
        out=values[0]
        for x in values[1:]: out=self.add(out,x)
        return out

    def norm(self,values,*,fixed=False):
        return self.sqrt(self.sum((self.mul(x,x) for x in values),fixed=fixed))


@dataclass(frozen=True)
class Reflector:
    essential: tuple[int,...]
    tau: int
    beta: int


def householder(a,x):
    tail_sq=a.sum(a.mul(z,z) for z in x[1:])
    if not less(MIN_NORMAL,tail_sq):
        return Reflector((ZERO,)*(len(x)-1),ZERO,x[0])
    beta=a.sqrt(a.add(a.mul(x[0],x[0]),tail_sq))
    if not less(x[0],ZERO): beta=neg(beta)
    denominator=a.sub(x[0],beta)
    essential=tuple(a.div(z,denominator) for z in x[1:])
    tau=a.div(a.sub(beta,x[0]),beta)
    return Reflector(essential,tau,beta)


def apply_householder(a,h,column):
    if zero(h.tau): return tuple(column)
    temporary=a.sum(a.mul(v,x) for v,x in zip(h.essential,column[1:]))
    temporary=a.add(temporary,column[0])
    first=a.sub(column[0],a.mul(h.tau,temporary))
    tail=tuple(a.sub(x,a.mul(a.mul(h.tau,v),temporary))
               for v,x in zip(h.essential,column[1:]))
    return (first,*tail)


@dataclass(frozen=True)
class QR:
    scale: int
    first_pivot: int
    reflectors: tuple[Reflector,Reflector]
    q: tuple[tuple[int,...],...]
    work: tuple[tuple[int,...],...]


def qr_precondition(a,v0,v1):
    scale=max((*map(absolute,v0),*map(absolute,v1)))
    if special(scale): raise ValueError('nonfinite Eigen seed input')
    if zero(scale): scale=ONE
    cols=[tuple(a.div(x,scale) for x in v) for v in (v0,v1)]
    norms=tuple(a.norm(c,fixed=True) for c in cols)
    pivot=1 if less(norms[0],norms[1]) else 0
    if pivot: cols.reverse()
    h0=householder(a,cols[0])
    c1=apply_householder(a,h0,cols[1])
    h1=householder(a,c1[1:])
    # The rectangular preconditioner writes adjoint(upper(R[0:2,0:2])).
    work=((h0.beta,ZERO),(c1[0],h1.beta))
    # HouseholderSequence::evalTo initializes I and applies H1 then H0.
    columns=[]
    for j in range(3):
        c=tuple(ONE if i==j else ZERO for i in range(3))
        c=(c[0],*apply_householder(a,h1,c[1:]))
        columns.append(apply_householder(a,h0,c))
    q=tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))
    return QR(scale,pivot,(h0,h1),q,work)


def plane(a,x,y,c,s):
    return (a.add(a.mul(c,x),a.mul(s,y)),
            a.add(a.mul(neg(s),x),a.mul(c,y)))


def left(a,m,p,q,c,s):
    out=[list(row) for row in m]
    for j in range(len(m[0])):
        out[p][j],out[q][j]=plane(a,m[p][j],m[q][j],c,s)
    return tuple(map(tuple,out))


def right(a,m,p,q,c,s):
    out=[list(row) for row in m]
    for i in range(len(m)):
        out[i][p],out[i][q]=plane(a,m[i][p],m[i][q],c,neg(s))
    return tuple(map(tuple,out))


def rotations(a,work):
    # Eigen's only loop pair is p=1,q=0: preserve this reversed extraction.
    m=((work[1][1],work[1][0]),(work[0][1],work[0][0]))
    t=a.add(m[0][0],m[1][1]); d=a.sub(m[1][0],m[0][1])
    if less(absolute(d),MIN_NORMAL):
        c0,s0=ONE,ZERO
    else:
        u=a.div(t,d)
        temporary=a.sqrt(a.add(ONE,a.mul(u,u)))
        s0=a.div(ONE,temporary); c0=a.div(u,temporary)
    m=left(a,m,0,1,c0,s0)
    x,y,z=m[0][0],m[0][1],m[1][1]
    denominator=a.mul(TWO,absolute(y))
    if less(denominator,MIN_NORMAL):
        cr,sr=ONE,ZERO
    else:
        tau=a.div(a.sub(x,z),denominator)
        w=a.sqrt(a.add(a.mul(tau,tau),ONE))
        t=a.div(ONE,a.add(tau,w) if less(ZERO,tau) else a.sub(tau,w))
        sign=ONE if less(ZERO,t) else neg(ONE)
        n=a.div(ONE,a.sqrt(a.add(a.mul(t,t),ONE)))
        sr=a.mul(a.mul(a.mul(neg(sign),a.div(y,absolute(y))),absolute(t)),n)
        cr=n
    # j_left = rot1 * j_right.transpose().
    cl=a.sub(a.mul(c0,cr),a.mul(s0,neg(sr)))
    sl=a.add(a.mul(c0,neg(sr)),a.mul(s0,cr))
    return cl,sl,cr,sr


@dataclass(frozen=True)
class Result:
    v0: tuple[int,...]
    v1: tuple[int,...]
    qr: QR
    work_history: tuple
    v_history: tuple
    thresholds: tuple[int,...]
    axis: tuple[int,...]
    operations: tuple[Op,...]
    profile: str=PROFILE

    @property
    def sweep_count(self): return len(self.work_history)-1

    def validate(self):
        expected=solve_bits(self.v0,self.v1)
        if self!=expected: raise ValueError('Eigen scalar result detached from input vectors')


def solve_bits(v0,v1):
    v0,v1=tuple(v0),tuple(v1)
    if len(v0)!=3 or len(v1)!=3: raise ValueError('two three-vectors required')
    for x in (*v0,*v1): M.value(x)
    a=Arithmetic(); qr=qr_precondition(a,v0,v1)
    work,v=qr.work,qr.q
    maxdiag=maximum(absolute(work[0][0]),absolute(work[1][1]))
    history=[work]; vs=[v]; thresholds=[]; seen=set()
    while True:
        threshold=maximum(MIN_NORMAL,a.mul(PRECISION,maxdiag))
        thresholds.append(threshold)
        if not (less(threshold,absolute(work[1][0])) or less(threshold,absolute(work[0][1]))):
            break
        state=(work,maxdiag)
        if state in seen:
            raise ArithmeticError('exact scalar Eigen Jacobi work-state cycle')
        seen.add(state)
        cl,sl,cr,sr=rotations(a,work)
        work=left(a,work,1,0,cl,sl)
        work=right(a,work,1,0,cr,sr)
        v=right(a,v,1,0,cr,sr)
        maxdiag=maximum(maxdiag,maximum(absolute(work[1][1]),absolute(work[0][0])))
        history.append(work); vs.append(v)
    axis=tuple(row[2] for row in v)
    if axis!=tuple(row[2] for row in qr.q): raise AssertionError('Jacobi corrupted the invariant axis column')
    return Result(v0,v1,qr,tuple(history),tuple(vs),tuple(thresholds),axis,tuple(a.operations))


def solve(v0,v1):
    return solve_bits(tuple(M.exact_bits(x) for x in v0),tuple(M.exact_bits(x) for x in v1))


def in_source_domain(v0,v1):
    """The normalized near-antiparallel branch, before QR scaling."""
    from tools.stability.ou3_alt_contraction import finite_seed_svd_roundoff as R
    if len(v0)!=3 or tuple(v1)!=(0,0,1): return False
    n2=sum(F(x)**2 for x in v0)
    delta=F(1,10**6)
    return 1-delta<n2<1+delta and v0[2]<R.B32.add(-1,R.B32.rn32(F(1,100000)))


def source_uniform_certificate():
    from tools.stability.ou3_alt_contraction import finite_seed_svd_roundoff as R
    theorem=R.build()
    assert theorem['universal_termination_promoted']
    return {'qualification':'OU3_ALT_EIGEN_340_NEAR_ANTIPARALLEL_TOTALITY_V1',
            'reviewed_source':REDUCTION.build(),
            'arithmetic_family':theorem['arithmetic_family'],
            'maximum_Jacobi_sweeps':theorem['maximum_Jacobi_sweeps'],
            'QR_axis_norm2_upper':theorem['QR']['axis_norm2_upper'],
            'source_normalizations':theorem['source_and_order']['normalization'],
            'second_sweep_residual_upper':theorem['second_sweep']['offdiagonal_error'],
            'retained_termination_threshold_lower':theorem['second_sweep']['threshold_lower'],
            'both_pivots_rank_one_rank_two_and_tiny_tails_covered':True,
            'source_uniform_QR_and_Jacobi_totality_closed':True,
            'axis_is_computed_QR_column_two':True,
            'local_same_tree_FMA_contractions_covered':True,
            'target_compiler_correspondence_proved_here':False}


def witness(v0,v1):
    if not in_source_domain(v0,v1):
        raise ValueError('Eigen seed witness requires the normalized near-antiparallel source domain')
    result=solve(v0,v1)
    if result.sweep_count>2:
        raise AssertionError('scalar solver contradicts the uniform two-sweep theorem')
    return START.svd_witness(v0,v1,tuple(M.value(x) for x in result.axis)),result


def readiness():
    return {
        'pinned_Eigen_source': REDUCTION.build(),
        'profile': PROFILE,
        'axis_computed_from_scaled_source_vectors_without_solver_input_port': True,
        'both_QR_pivots_and_householder_zero_tail_branches_materialized': True,
        'Jacobi_2x2_operation_and_termination_predicate_graph_materialized': True,
        'discarded_Jacobi_nonfinite_arithmetic_retained': True,
        'per_input_result_is_universal_termination_certificate': False,
        'source_uniform_Jacobi_termination_closed': True,
        'source_uniform_totality_certificate':source_uniform_certificate(),
        'target_compiler_correspondence_closed': False,
        'storage_search_allowed': False,
    }
