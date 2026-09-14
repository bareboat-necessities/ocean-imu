"""Source-bound first-sample private Mahony update under the scalar profile.

The ordinary Eigen FromTwoVectors branch is evaluated from the SAME rounded
accelerometer packet, including both normalizations (acc/acc.norm(), then
FromTwoVectors' a.normalized()).  No final quaternion or reciprocal can be
supplied by the caller.  A too-small first sample preserves the observer; the
outer guard/LPF/statistics event still advances.

This is a conditional program identity, not startup admission.  The nearly
opposite-vector JacobiSVD branch is materialized with an explicit solver-axis
witness; target Eigen solver correspondence, sqrt/library/profile qualification
and universal source-domain bounds remain fail-closed.  It is not legitimate to
exclude that branch from COMPLETE-BRMM merely because the solver correspondence
is not yet certified.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


class UnqualifiedSeedBranch(ValueError):
    """An admitted source may take a branch not yet qualified by this module."""


@dataclass(frozen=True)
class SVDWitness:
    """Literal Eigen JacobiSVD axis output for the near-opposite branch.

    The axis is the actual binary32 ``matrixV().col(2)`` result.  We retain its
    nullspace/unit residuals rather than assuming an exact symbolic cross axis.
    This binds the branch topology and downstream quaternion arithmetic while
    leaving target Eigen solver correspondence/source-uniform residual bounds
    as explicit deployment obligations.
    """
    axis:tuple
    v0_dot:F
    v1_dot:F
    norm2_minus_one:F
    def __post_init__(self):
        if len(self.axis)!=3: raise ValueError('JacobiSVD axis must be vec3')
        axis=tuple(F(x) for x in self.axis)
        for x in axis: M.exact_bits(x)
        object.__setattr__(self,'axis',axis)
        for name in ('v0_dot','v1_dot','norm2_minus_one'):
            object.__setattr__(self,name,F(getattr(self,name)))


def svd_witness(v0,v1,axis):
    """Attach, but do not certify, the actual Eigen solver axis."""
    a=tuple(F(x) for x in axis)
    if len(a)!=3: raise ValueError('JacobiSVD axis must be vec3')
    for x in a: M.exact_bits(x)
    return SVDWitness(a,sum(x*y for x,y in zip(v0,a)),
                      sum(x*y for x,y in zip(v1,a)),sum(x*x for x in a)-1)


def rn(x):
    return M.value(M.round_bits(F(x)))


def sqrt32(x):
    """Exact RNE binary32 sqrt using rational squared-midpoint comparisons.

    Binary search finds the two adjacent nonnegative binary32 words bracketing
    sqrt(x).  Squaring their arithmetic midpoint preserves order, so the final
    comparison is exact even when sqrt(x) is irrational.  Zero/subnormals are
    retained.  This specifies correctly rounded sqrt; it does not assert that
    the deployed libm has this property.
    """
    x=F(x)
    if x<0: raise ValueError('negative startup sqrt radicand')
    M.exact_bits(x)
    if not x: return F(0)
    lo,hi=0,M.MAX_FINITE
    while lo<hi:
        mid=(lo+hi+1)//2
        if M.value(mid)**2<=x: lo=mid
        else: hi=mid-1
    y=M.value(lo)
    if y*y==x: return y
    if lo==M.MAX_FINITE: raise OverflowError('startup sqrt overflow')
    upper=M.value(lo+1); mid2=((y+upper)/2)**2
    return upper if x>mid2 or (x==mid2 and lo&1) else y


@dataclass(frozen=True)
class Operation:
    kind:str
    operands:tuple
    result:F
    residual:F

    def validate(self):
        for x in (*self.operands,self.result): M.exact_bits(x)
        if self.kind=='sqrt':
            if len(self.operands)!=1: raise ValueError('sqrt arity')
            x=self.operands[0]; b=M.exact_bits(self.result)
            if x<0 or self.result<0: raise ValueError('nonnegative sqrt relation required')
            y=self.result
            lower=(y+M.value(b-1))/2 if b else F(0)
            upper=(y+M.value(b+1))/2 if b<M.MAX_FINITE else M.pow2(128)
            even=not(b&1)
            if not ((x>lower*lower or (x==lower*lower and even)) and
                    (x<upper*upper or (x==upper*upper and even))):
                raise ValueError('startup sqrt outside exact squared midpoint cell')
            residual=y*y-x
        else:
            if len(self.operands)!=2: raise ValueError('binary arithmetic arity')
            a,b=self.operands
            if self.kind=='add': exact=a+b
            elif self.kind=='sub': exact=a-b
            elif self.kind=='mul': exact=a*b
            elif self.kind=='div': exact=a/b
            else: raise ValueError('unknown seed operation')
            if not M.rounding_cell_contains(exact,M.exact_bits(self.result)):
                raise ValueError('seed operation outside exact RNE cell')
            residual=self.result-exact
        if self.residual!=residual: raise ValueError('seed operation residual detached')


class Arithmetic:
    def __init__(self): self.operations=[]

    def op(self,kind,*args):
        args=tuple(F(x) for x in args)
        for x in args: M.exact_bits(x)
        if kind=='sqrt':
            out=sqrt32(args[0]); residual=out*out-args[0]
        else:
            a,b=args
            if kind=='add': exact=a+b
            elif kind=='sub': exact=a-b
            elif kind=='mul': exact=a*b
            elif kind=='div': exact=a/b
            else: raise ValueError('unknown seed operation')
            out=rn(exact); residual=out-exact
        self.operations.append(Operation(kind,args,out,residual))
        return out

    def dot3(self,a,b):
        p=tuple(self.op('mul',x,y) for x,y in zip(a,b))
        return self.op('add',p[0],self.op('add',p[1],p[2]))


@dataclass(frozen=True)
class Seed:
    acc:tuple
    norm:F
    normalized_down:tuple|None
    quaternion:tuple|None
    branch:str
    operations:tuple[Operation,...]
    svd:SVDWitness|None=None

    def validate(self):
        for op in self.operations: op.validate()
        if self!=seed(self.acc,svd=self.svd): raise ValueError('startup seed detached from same accelerometer')


def seed(acc,*,svd:SVDWitness|None=None):
    if len(acc)!=3: raise ValueError('startup seed requires vec3 accelerometer')
    acc=tuple(F(x) for x in acc)
    for x in acc: M.exact_bits(x)
    a=Arithmetic(); op=a.op
    norm=op('sqrt',a.dot3(acc,acc))
    if not norm>rn(F(1,1000)):
        if svd is not None: raise ValueError('small-accel branch consumes no JacobiSVD witness')
        return Seed(acc,norm,None,None,'acc-norm-too-small',tuple(a.operations),None)
    down=tuple(-op('div',x,norm) for x in acc)
    down_norm=op('sqrt',a.dot3(down,down))
    v0=tuple(op('div',x,down_norm) for x in down)
    z=(F(0),F(0),F(1))
    z_norm=op('sqrt',a.dot3(z,z))
    v1=tuple(op('div',x,z_norm) for x in z)
    c=a.dot3(v1,v0)
    cutoff=op('add',-1,rn(F(1,100000)))
    if c<cutoff:
        if svd is None:
            raise UnqualifiedSeedBranch('near-antiparallel Eigen JacobiSVD startup seed requires explicit solver witness')
        if not isinstance(svd,SVDWitness): raise TypeError('JacobiSVD solver witness required')
        if svd != svd_witness(v0,v1,svd.axis):
            raise ValueError('JacobiSVD witness residuals detached from same normalized vectors')
        cc=max(c,F(-1))
        w2=op('mul',op('add',1,cc),F(1,2))
        w=op('sqrt',w2)
        scale=op('sqrt',op('sub',1,w2))
        vec=tuple(op('mul',x,scale) for x in svd.axis)
        q=(w,*vec)
        return Seed(acc,norm,v0,q,'near-antiparallel-JacobiSVD',tuple(a.operations),svd)
    if svd is not None: raise ValueError('ordinary FromTwoVectors branch consumes no JacobiSVD witness')
    axis=tuple(op('sub',op('mul',v0[j],v1[k]),op('mul',v0[k],v1[j]))
               for j,k in ((1,2),(2,0),(0,1)))
    s=op('sqrt',op('mul',op('add',1,c),2))
    invs=op('div',1,s)
    vec=tuple(op('mul',x,invs) for x in axis)
    q=(op('mul',s,F(1,2)),*vec)
    return Seed(acc,norm,v0,q,'ordinary-FromTwoVectors',tuple(a.operations),None)


def step(state:V.State,cfg:V.Config,*,dt,gyro,acc,svd:SVDWitness|None=None,profile=M.PROFILE):
    if not isinstance(state,V.State) or not isinstance(cfg,V.Config):
        raise TypeError('private observer State/Config required')
    if profile!=M.PROFILE: raise ValueError('unqualified compiler/evaluation profile')
    M.exact_bits(dt)
    if F(dt)<=0: raise ValueError('positive finite startup dt required')
    if len(gyro)!=3 or len(acc)!=3: raise ValueError('three-dimensional IMU packet required')
    for x in (*gyro,*acc,*state.q,*state.integral,state.elapsed,state.up,
              cfg.two_kp,cfg.two_ki,cfg.gravity,cfg.settle_sec): M.exact_bits(x)
    if state.initialized:
        if svd is not None: raise ValueError('initialized Mahony step consumes no startup JacobiSVD witness')
        return M.step_initialized(state,cfg,dt=dt,gyro=gyro,acc=acc,profile=profile)
    if state!=V.State():
        raise ValueError('uninitialized private observer must retain literal reset state')
    initial=seed(acc,svd=svd)
    if initial.quaternion is None:
        return M.StepResult(V.Result(state,state.up),(),(),startup=initial)
    initialized=replace(state,q=initial.quaternion,initialized=True)
    result=M.step_initialized(initialized,cfg,dt=dt,gyro=gyro,acc=acc,profile=profile)
    return replace(result,startup=initial)


def readiness():
    return {
        'ordinary_FromTwoVectors_seed_computed_from_same_binary32_accelerometer':True,
        'both_shipping_accelerometer_normalizations_retained':True,
        'first_valid_sample_executes_Mahony_after_seeding_in_same_event':True,
        'too_small_first_accel_preserves_uninitialized_observer':True,
        'seed_quaternion_or_reciprocal_input_port_present':False,
        'near_antiparallel_JacobiSVD_branch_topology_materialized_with_solver_witness':True,
        'near_antiparallel_JacobiSVD_solver_correspondence_qualified':False,
        'target_sqrt_Eigen_and_compiler_correspondence_closed':False,
        'every_admitted_startup_history_covered':False,
        'source_uniform_complete_600_step_word_qualified':False,
        'storage_search_allowed':False,
        'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
