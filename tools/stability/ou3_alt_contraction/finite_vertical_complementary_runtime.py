"""Finite private Mahony/vertical-acceleration runtime for ALT.

The deployed default WPE and sigma channel share one VerticalAccelComplementary
observer.  It consumes the same raw gyro/accelerometer sample before the MEKF
uses it.  This module materializes the initialized IMU-only Mahony recurrence,
including accel normalization, proportional/integral branches, quaternion Euler
integration/normalization and the levelled up-acceleration readout.

The float fast-inverse-sqrt results and first-sample FromTwoVectors seed remain
explicit runtime witnesses.  They are not silently replaced by exact sqrt or a
free attitude.  Euler outputs are omitted because VerticalAccelComplementary
computes them but does not consume them for its state or vertical output.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)
def dot(a,b): return sum((a[i]*b[i] for i in range(len(a))),F(0))


@dataclass(frozen=True)
class Config:
    two_kp:F
    two_ki:F
    gravity:F
    settle_sec:F
    def __post_init__(self):
        for n in ('two_kp','two_ki','gravity','settle_sec'): object.__setattr__(self,n,R(getattr(self,n)))
        if self.two_kp<0 or self.two_ki<0 or self.gravity<=0 or self.settle_sec<0: raise ValueError('invalid complementary configuration')


@dataclass(frozen=True)
class State:
    q:tuple=(F(1),F(0),F(0),F(0))
    integral:tuple=(F(0),F(0),F(0))
    initialized:bool=False
    elapsed:F=F(0)
    up:F=F(0)
    def __post_init__(self):
        q=tuple(R(x) for x in self.q); integ=tuple(R(x) for x in self.integral)
        if len(q)!=4 or len(integ)!=3: raise ValueError('q/integral dimensions')
        object.__setattr__(self,'q',q); object.__setattr__(self,'integral',integ); object.__setattr__(self,'elapsed',R(self.elapsed)); object.__setattr__(self,'up',R(self.up))
        if not isinstance(self.initialized,bool): raise TypeError('literal initialization flag required')
        if self.elapsed<0: raise ValueError('negative elapsed')


@dataclass(frozen=True)
class InvSqrtWitness:
    input_norm_sq:F
    reciprocal:F
    def __post_init__(self):
        n,r=R(self.input_norm_sq),R(self.reciprocal)
        if n<=0 or r<=0: raise ValueError('positive inverse-sqrt witness required')
        object.__setattr__(self,'input_norm_sq',n); object.__setattr__(self,'reciprocal',r)


@dataclass(frozen=True)
class SeedWitness:
    q:tuple
    def __post_init__(self):
        q=tuple(R(x) for x in self.q)
        if len(q)!=4: raise ValueError('seed quaternion dimension')
        object.__setattr__(self,'q',q)


@dataclass(frozen=True)
class Result:
    state:State
    vertical_accel:F


def step(s:State,cfg:Config,*,dt,gyro,acc,
         accel_invnorm:InvSqrtWitness|None=None,
         quat_invnorm:InvSqrtWitness|None=None,
         seed:SeedWitness|None=None):
    dt=R(dt); g=tuple(R(x) for x in gyro); a_raw=tuple(R(x) for x in acc)
    if dt<=0 or len(g)!=3 or len(a_raw)!=3: raise ValueError('positive dt and 3-vectors required')
    q=s.q; initialized=s.initialized
    acc_sq=dot(a_raw,a_raw)
    # VerticalAccelComplementary returns without mutation if the first accel norm
    # is <=1e-3.  Squared threshold is 1e-6.
    if not initialized:
        if acc_sq<=F(1,10**6):
            if any(x is not None for x in (accel_invnorm,quat_invnorm,seed)): raise ValueError('unseeded invalid-acc branch consumes no witnesses')
            return Result(s,s.up)
        if seed is None: raise ValueError('first valid sample requires FromTwoVectors seed witness')
        q=seed.q; initialized=True
    elif seed is not None:
        raise ValueError('initialized branch consumes no seed witness')

    # Mahony is fed -acc. Feedback is skipped only for exactly all-zero accel.
    ax,ay,az=(-a_raw[0],-a_raw[1],-a_raw[2])
    gx,gy,gz=g
    integ=list(s.integral)
    if not (ax==0 and ay==0 and az==0):
        if accel_invnorm is None or accel_invnorm.input_norm_sq!=acc_sq:
            raise ValueError('accelerometer inverse-sqrt witness detached from SAME raw sample')
        rn=accel_invnorm.reciprocal; ax*=rn; ay*=rn; az*=rn
        q0,q1,q2,q3=q
        hvx=q1*q3-q0*q2
        hvy=q0*q1+q2*q3
        hvz=F(1,2)*(q0*q0-q1*q1-q2*q2+q3*q3)
        hex_=ay*hvz-az*hvy; hey=az*hvx-ax*hvz; hez=ax*hvy-ay*hvx
        if cfg.two_ki>0:
            integ[0]+=cfg.two_ki*hex_*dt; integ[1]+=cfg.two_ki*hey*dt; integ[2]+=cfg.two_ki*hez*dt
            gx+=integ[0]; gy+=integ[1]; gz+=integ[2]
        else:
            integ=[F(0),F(0),F(0)]
        gx+=cfg.two_kp*hex_; gy+=cfg.two_kp*hey; gz+=cfg.two_kp*hez
    else:
        if accel_invnorm is not None: raise ValueError('zero-acc Mahony branch consumes no accel normalization witness')

    half=F(1,2)*dt; gx*=half; gy*=half; gz*=half
    q0,q1,q2,q3=q; qa,qb,qc=q0,q1,q2
    un=(q0+(-qb*gx-qc*gy-q3*gz),
        q1+(qa*gx+qc*gz-q3*gy),
        q2+(qa*gy-qb*gz+q3*gx),
        q3+(qa*gz+qb*gy-qc*gx))
    qsq=dot(un,un)
    if quat_invnorm is None or quat_invnorm.input_norm_sq!=qsq:
        raise ValueError('quaternion inverse-sqrt witness detached from SAME Euler step')
    qr=quat_invnorm.reciprocal
    qn=tuple(v*qr for v in un)
    w,x,y,z=qn
    down=(2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z)
    up=-(dot(down,a_raw)+cfg.gravity)
    nxt=State(qn,tuple(integ),initialized,s.elapsed+dt,up)
    return Result(nxt,up)


def readiness():
    return {
      'raw_gyro_accel_same_sample_inputs_retained':True,
      'mahony_feedback_and_quaternion_recurrence_materialized':True,
      'vertical_levelled_output_materialized':True,
      'first_sample_FromTwoVectors_seed_attached':False,
      'fast_inv_sqrt_binary32_attached':False,
      'raw_sensor_BRMM_disturbance_relation_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
