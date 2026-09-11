"""Exact finite real-arithmetic prediction graph for ALT joint24.

This module closes the deterministic *mean/error algebra* for one prediction
step.  It does not certify source reachability, covariance propagation,
finite precision, or a storage inequality.

The finite attitude error uses homogeneous quaternions.  If

  q_s = shipping step for (-omega_hat + e_bg) h,
  q_e ~ (2,c) for current Cayley error,
  q_n = shipping step for (-omega_hat) h,

then the relative attitude for that SAMPLED SHADOW is represented
exactly by

  q_plus ~ q_s * q_e * conj(q_n),
  c_plus = 2 v_plus / w_plus,

provided w_plus != 0.  Normalization of each step quaternion cancels in this
projective identity, so both deployed quat_from_delta_theta branches can be
retained as hard graph variables rather than linearized.

This is physical attitude only if an additional source identity establishes
q_s as the continuous physical increment. In general it does not. Use
finite_physical_prediction for the actual increment and retain its angular
defect; do not silently promote this sampled shadow to physical truth.

Translation uses the already-proved same-history q15 forcing recurrence, and
accelerometer-bias truth/error uses one shared physical driver w_b.  No source
coordinate is independently re-boxed here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from typing import Sequence

NSTATE=24


def rational(x):
    if isinstance(x,(float,bool)):
        raise TypeError('exact graph accepts integers, Fraction or rational strings')
    return F(x)


def vec(x,n):
    if len(x)!=n: raise ValueError(f'expected vector length {n}')
    return [rational(v) for v in x]


def quat_mul(a,b):
    aw,*av=a; bw,*bv=b
    ax,ay,az=av; bx,by,bz=bv
    return [
      aw*bw-ax*bx-ay*by-az*bz,
      aw*bx+bw*ax+ay*bz-az*by,
      aw*by+bw*ay+az*bx-ax*bz,
      aw*bz+bw*az+ax*by-ay*bx,
    ]


def quat_conj(q): return [q[0],-q[1],-q[2],-q[3]]


def relative_attitude_prediction(c, q_shadow, q_nominal):
    """Exact projective finite Cayley prediction.

    q_shadow and q_nominal are any nonzero scalar multiples of the deployed
    step quaternions. Their normalization therefore need not be repeated.
    """
    c=vec(c,3); qs=vec(q_shadow,4); qn=vec(q_nominal,4)
    qe=[F(2),*c]
    qp=quat_mul(quat_mul(qs,qe),quat_conj(qn))
    if qp[0]==0: raise ValueError('predicted relative attitude hits Cayley pole')
    return [F(2)*qp[i]/qp[0] for i in range(1,4)]


def translation_axis_coefficients(phi_va,phi_pa,phi_Sa,alpha,h):
    return tuple(map(rational,(phi_va,phi_pa,phi_Sa,alpha,h)))


def translation_prediction(error12,q15,coeff):
    """Exact [e_v,e_p,e_S,e_aw] xyz update from one correlated q15 witness.

    error12 block order: v_xyz,p_xyz,S_xyz,aw_xyz.
    q15 order: a0_xyz,a1_xyz,J0_xyz,J1_xyz,J2_xyz.
    """
    e=vec(error12,12); q=vec(q15,15)
    phi_va,phi_pa,phi_Sa,alpha,h=coeff
    out=[F(0)]*12
    h2=h*h
    for a in range(3):
        ev,ep,eS,ea=e[a],e[3+a],e[6+a],e[9+a]
        a0,a1,J0,J1,J2=q[a],q[3+a],q[6+a],q[9+a],q[12+a]
        out[a]=ev+phi_va*ea+(J0-phi_va*a0)
        out[3+a]=ep+h*ev+phi_pa*ea+(J1-phi_pa*a0)
        out[6+a]=eS+h*ep+h2*ev/F(2)+phi_Sa*ea+(J2-phi_Sa*a0)
        out[9+a]=alpha*ea+(a1-alpha*a0)
    return out


def bias_prediction(e_ba,beta,w_bias,phi_hat,phi_true):
    """One-history exact recurrence with the SAME physical driver in both rows."""
    e=vec(e_ba,3); b=vec(beta,3); w=vec(w_bias,3)
    ph=rational(phi_hat); pt=rational(phi_true)
    return ([ph*e[i]+(pt-ph)*b[i]+w[i] for i in range(3)],
            [pt*b[i]+w[i] for i in range(3)])


@dataclass(frozen=True)
class FinitePrediction:
    state_after:list
    attitude_denominator:F


def prediction(mode,z,*,q_shadow,q_nominal,q15,coeff,w_bias,phi_true,phi_hat=None):
    """Exact joint24 mean/error prediction conditional on hard source graphs."""
    if mode not in ('H','A'): raise ValueError('mode must be H or A')
    z=vec(z,24)
    cplus=relative_attitude_prediction(z[:3],q_shadow,q_nominal)
    qp=quat_mul(quat_mul(vec(q_shadow,4),[F(2),*z[:3]]),quat_conj(vec(q_nominal,4)))
    trans=translation_prediction(z[6:18],q15,coeff)
    out=[F(0)]*24
    out[:3]=cplus
    out[3:6]=z[3:6]
    out[6:18]=trans
    ph=F(1) if mode=='H' else rational(phi_hat)
    eba,beta=bias_prediction(z[18:21],z[21:24],w_bias,ph,phi_true)
    out[18:21]=eba;out[21:24]=beta
    return FinitePrediction(out,qp[0])


def readiness():
    return {
      'map_representation':'finite_physical_prediction_descriptor',
      'finite_attitude_prediction_identity':True,
      'finite_translation_q15_forcing_identity':True,
      'finite_shared_bias_driver_identity':True,
      'source_uniform_step_quaternion_graph_attached':False,
      'source_uniform_q15_same_history_relation_attached':False,
      'source_uniform_bias_parameter_root_attached':False,
      'covariance_successor_attached':False,
      'finite_precision_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
