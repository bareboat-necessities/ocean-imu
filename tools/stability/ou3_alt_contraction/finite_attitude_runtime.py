"""Finite attitude/gyro-bias covariance runtime descriptor for ALT.

This follows the shipping with_gyro_bias=true path.  It removes arbitrary
``F_AA,Q_AA`` at the next composition layer by reconstructing them from angular
rate, Qbase and explicit branch witnesses.  Trigonometric values and Eigen LDLT
outcomes are still runtime/source/finite-precision obligations; their presence
is explicit and cannot be mistaken for a universal certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P

EPS_PSD = F(1,10**12)
ISO_TOL = F(1,10**9)
RATE_SMALL2 = F(1,10**14)  # (1e-7)^2


def skew(v): return M.skew(P.vec(v,3))


@dataclass(frozen=True)
class TrigWitness:
    """One nonzero-rate sin/cos evaluation at a declared time."""
    t: F
    sin_theta: F
    cos_theta: F
    inv_omega: F
    def __post_init__(self):
        t,s,c,inv=map(P.rational,(self.t,self.sin_theta,self.cos_theta,self.inv_omega))
        if t < 0 or inv <= 0: raise ValueError('valid time and positive inverse rate required')
        if s*s+c*c != 1: raise ValueError('trig witness must lie on unit circle')
        object.__setattr__(self,'t',t); object.__setattr__(self,'sin_theta',s); object.__setattr__(self,'cos_theta',c); object.__setattr__(self,'inv_omega',inv)


@dataclass(frozen=True)
class AngularRuntime:
    w: tuple
    h: F
    full: TrigWitness | None = None
    half: TrigWitness | None = None
    def __post_init__(self):
        w=tuple(P.vec(self.w,3)); h=P.rational(self.h)
        if h <= 0: raise ValueError('positive attitude step required')
        w2=M.dot(w,w); small=w2 < RATE_SMALL2
        if small:
            if self.full is not None or self.half is not None: raise ValueError('small-rate branch consumes no trig witnesses')
        else:
            if not isinstance(self.full,TrigWitness) or not isinstance(self.half,TrigWitness): raise ValueError('general branch requires full and half trig witnesses')
            if self.full.t != h or self.half.t != h/2: raise ValueError('trig witness times detached from same step')
            for q in (self.full,self.half):
                if q.inv_omega*q.inv_omega*w2 != 1: raise ValueError('inverse-rate witness detached from same angular rate')
        object.__setattr__(self,'w',w); object.__setattr__(self,'h',h)
    @property
    def small_rate(self): return M.dot(self.w,self.w) < RATE_SMALL2

    def rot_B(self,t, witness=None):
        t=P.rational(t); W=skew(self.w); W2=M.mm(W,W); I=M.eye(3)
        if self.small_rate:
            if witness is not None: raise ValueError('small-rate helper consumes no trig witness')
            R=M.plus(M.plus(I,M.scaled(W,-t)),M.scaled(W2,t*t/2))
            B=M.plus(M.plus(M.scaled(I,t),M.scaled(W,-t*t/2)),M.scaled(W2,t*t*t/6))
            return R,B
        q=witness
        if not isinstance(q,TrigWitness) or q.t != t: raise ValueError('same-time trig witness required')
        K=M.scaled(W,q.inv_omega); K2=M.mm(K,K); inv2=q.inv_omega*q.inv_omega
        R=M.plus(M.plus(I,M.scaled(K,-q.sin_theta)),M.scaled(K2,1-q.cos_theta))
        B=M.plus(M.plus(M.scaled(I,t),M.scaled(W,-(1-q.cos_theta)*inv2)),M.scaled(W2,(t-q.sin_theta*q.inv_omega)*inv2))
        return R,B

    def integral_B(self):
        W=skew(self.w); W2=M.mm(W,W); I=M.eye(3); h=self.h
        if self.small_rate:
            return M.plus(M.plus(M.scaled(I,h*h/2),M.scaled(W,-h**3/6)),M.scaled(W2,h**4/24))
        q=self.full; inv2=q.inv_omega*q.inv_omega
        return M.plus(M.plus(M.scaled(I,h*h/2),M.scaled(W,-(h-q.sin_theta*q.inv_omega)*inv2)),M.scaled(W2,(h*h/2+(q.cos_theta-1)*inv2)*inv2))

    def F_AA(self):
        R,B=self.rot_B(self.h,None if self.small_rate else self.full); out=M.eye(6)
        for i in range(3):
            for j in range(3): out[i][j]=R[i][j]; out[i][3+j]=B[i][j]
        return out


def is_isotropic_shipping(Q):
    Q=M.mat(Q,3,3); a,b,c=Q[0][0],Q[1][1],Q[2][2]; mean=(a+b+c)/3
    off=sum(abs(Q[i][j]) for i in range(3) for j in range(3) if i!=j)
    return abs(a-mean)+abs(b-mean)+abs(c-mean)+off <= ISO_TOL*(1+abs(mean))


def simpson_R(runtime:AngularRuntime,Q):
    Q=M.mat(Q,3,3); I=M.eye(3); Rm,_=runtime.rot_B(runtime.h/2,None if runtime.small_rate else runtime.half); R1,_=runtime.rot_B(runtime.h,None if runtime.small_rate else runtime.full)
    return M.scaled(M.plus(M.plus(Q,M.scaled(M.mm(M.mm(Rm,Q),M.transpose(Rm)),4)),M.mm(M.mm(R1,Q),M.transpose(R1))),runtime.h/6)


def simpson_B(runtime:AngularRuntime,Q):
    Q=M.mat(Q,3,3); _,Bm=runtime.rot_B(runtime.h/2,None if runtime.small_rate else runtime.half); _,B1=runtime.rot_B(runtime.h,None if runtime.small_rate else runtime.full)
    return M.scaled(M.plus(M.scaled(M.mm(M.mm(Bm,Q),M.transpose(Bm)),4),M.mm(M.mm(B1,Q),M.transpose(B1))),runtime.h/6)


def structured_Q_pre_hygiene(runtime:AngularRuntime,Qbase):
    Q=M.mat(Qbase,6,6); Qg=[r[:3] for r in Q[:3]]; Qbg=[r[3:] for r in Q[3:]]
    if is_isotropic_shipping(Qg): IR=M.scaled(M.eye(3),Qg[0][0]*runtime.h)
    else: IR=simpson_R(runtime,Qg)
    IBB=simpson_B(runtime,Qbg); Qtt=M.plus(IR,IBB); Qbb=M.scaled(Qbg,runtime.h); Qtb=M.mm(runtime.integral_B(),Qbg)
    out=M.zeros(6,6)
    for i in range(3):
        for j in range(3):
            out[i][j]=Qtt[i][j]; out[i][3+j]=Qtb[i][j]; out[3+i][j]=Qtb[j][i]; out[3+i][3+j]=Qbb[i][j]
    return M.scaled(M.plus(out,M.transpose(out)),F(1,2))


def psd_hygiene_6(S, *, first_ldlt_success, second_ldlt_success=None, eps=EPS_PSD):
    """Literal finite-value 6x6 project_psd_ou_iii control branches.

    Nonfinite replacement is a deployment-roundoff obligation and therefore not
    silently modeled here: exact rational S is finite.  LDLT success booleans
    remain explicit runtime witnesses.
    """
    S=M.mat(S,6,6); eps=P.rational(eps)
    if eps <= 0 or not isinstance(first_ldlt_success,bool): raise ValueError('positive eps and literal first LDLT outcome required')
    out=M.scaled(M.plus(S,M.transpose(S)),F(1,2))
    if first_ldlt_success:
        if second_ldlt_success is not None: raise ValueError('unused second LDLT outcome')
        return out
    lb=[]
    for i in range(6): lb.append(out[i][i]-sum(abs(out[i][j]) for j in range(6) if j!=i))
    min_lb=min(lb)
    if not (min_lb > eps):
        shift=eps-min_lb
        for i in range(6): out[i][i]+=shift
    if not isinstance(second_ldlt_success,bool): raise TypeError('failed first LDLT requires second outcome')
    if not second_ldlt_success:
        for i in range(6): out[i][i]+=10*eps
    return M.scaled(M.plus(out,M.transpose(out)),F(1,2))


def attitude_blocks(runtime:AngularRuntime,Qbase,*,use_exact_Q=True,
                    first_ldlt_success=True,second_ldlt_success=None):
    """Return literal shipping F_AA,Q_AA for a declared Q branch."""
    if not isinstance(use_exact_Q,bool): raise TypeError('literal exact/fast Q branch required')
    Q=M.mat(Qbase,6,6)
    if Q != M.transpose(Q): raise ValueError('symmetric Qbase required')
    if use_exact_Q:
        q=structured_Q_pre_hygiene(runtime,Q)
        q=psd_hygiene_6(q,first_ldlt_success=first_ldlt_success,second_ldlt_success=second_ldlt_success)
    else:
        if second_ldlt_success is not None: raise ValueError('fast Q branch has no PSD retry outcome')
        q=M.scaled(Q,runtime.h)
    return runtime.F_AA(),q


def readiness():
    return {
      'F_AA_constant_rate_structure_materialized':True,
      'structured_Q_AA_simpson_and_integral_B_materialized':True,
      'fast_Q_AA_branch_materialized':True,
      'isotropic_Qg_branch_predicate_materialized':True,
      'Q_AA_psd_hygiene_control_branches_materialized':True,
      'trig_values_same_runtime_source_attached':False,
      'Eigen_LDLT_outcomes_finite_precision_attached':False,
      'nonfinite_replacement_finite_precision_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
