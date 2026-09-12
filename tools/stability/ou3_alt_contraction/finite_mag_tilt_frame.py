"""Exact-real yaw-stripped boat quaternion used by startup MagAutoTuner.

Shipping's wrapper computes

  q_bw <- normalize(q_bw_in)
  c = R(q_bw)[0,0]
  s = R(q_bw)[1,0]
  yaw = atan2(s,c)
  q_tilt = normalize(AngleAxis(-yaw, z) * q_bw)

before each startup magnetic accumulation sample.  This module removes a free
``q_tilt`` operand from the ALT word.  The atan2/AngleAxis binary32 execution is
represented by an exact half-angle witness tied to the SAME (c,s); deployment
transcendental/rounding correspondence remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return M.rational(x)


@dataclass(frozen=True)
class SqrtWitness:
    radicand:F
    value:F
    def __post_init__(self):
        a,v=R(self.radicand),R(self.value)
        if a<0 or v<0 or v*v!=a: raise ValueError('exact nonnegative sqrt witness required')
        object.__setattr__(self,'radicand',a); object.__setattr__(self,'value',v)


@dataclass(frozen=True)
class YawHalfWitness:
    """Exact ideal-real AngleAxis half-angle tied to atan2 direction (c,s)."""
    c:F
    s:F
    radial:SqrtWitness
    cos_half:F
    sin_half:F
    def __post_init__(self):
        c,s,ch,sh=map(R,(self.c,self.s,self.cos_half,self.sin_half))
        if not isinstance(self.radial,SqrtWitness) or self.radial.radicand!=c*c+s*s:
            raise ValueError('yaw radial witness detached from same rotation-column direction')
        r=self.radial.value
        if r<=0: raise ValueError('positive heading radial required for nondegenerate yaw witness')
        cy,sy=c/r,s/r
        if ch<0 or ch*ch!=(1+cy)/2 or sh*sh!=(1-cy)/2 or 2*ch*sh!=sy:
            raise ValueError('half-angle witness detached from atan2 direction')
        object.__setattr__(self,'c',c); object.__setattr__(self,'s',s)
        object.__setattr__(self,'cos_half',ch); object.__setattr__(self,'sin_half',sh)


@dataclass(frozen=True)
class Result:
    q_input:tuple
    q_normalized:tuple
    heading_c:F
    heading_s:F
    q_tilt:tuple


def yaw_removed(q_bw_in,*,q_norm:SqrtWitness,yaw_half:YawHalfWitness|None):
    q=tuple(M.vec(q_bw_in,4)); q2=M.dot(q,q)
    if not isinstance(q_norm,SqrtWitness) or q_norm.radicand!=q2:
        raise ValueError('boat quaternion norm witness detached from same input')
    if q_norm.value<=F(1,10**6):
        # Shipping returns Identity for nonfinite/tiny input.  Rational inputs
        # cover the tiny branch; nonfinite deployment values remain separate.
        return Result(q,(1,0,0,0),F(1),F(0),(1,0,0,0))
    qn=tuple(x/q_norm.value for x in q)
    if M.dot(qn,qn)!=1: raise AssertionError('normalized boat quaternion not unit')
    w,x,y,z=qn
    c=w*w+x*x-y*y-z*z
    s=2*(x*y+w*z)

    if c==0 and s==0:
        # Ideal-real representative of std::atan2(0,0) -> 0.  Signed-zero
        # deployment details are intentionally left to binary32 correspondence.
        if yaw_half is not None: raise ValueError('degenerate zero heading consumes no yaw witness')
        return Result(q,qn,c,s,qn)

    if not isinstance(yaw_half,YawHalfWitness) or yaw_half.c!=c or yaw_half.s!=s:
        raise ValueError('yaw half-angle witness detached from same normalized boat quaternion')
    q_yaw_inv=(yaw_half.cos_half,F(0),F(0),-yaw_half.sin_half)
    qt=tuple(P.quat_mul(q_yaw_inv,qn))
    # Both factors are exactly unit in the ideal-real relation.  The shipping
    # post-product normalize is therefore identity here; its binary32 effect is
    # reserved for the deployment correspondence layer.
    if M.dot(qt,qt)!=1: raise AssertionError('yaw removal product lost unit norm')
    return Result(q,qn,c,s,qt)


def readiness():
    return {
      'startup_tilt_quaternion_is_not_free_operand':True,
      'boat_quaternion_normalization_relation_materialized':True,
      'heading_c_R00_and_s_R10_from_same_normalized_quaternion':True,
      'inverse_yaw_half_angle_tied_to_same_heading_direction':True,
      'left_multiply_inverse_yaw_then_normalize_relation_materialized':True,
      'atan2_AngleAxis_binary32_attached':False,
      'quaternion_normalization_binary32_attached':False,
      'input_allFinite_binary32_attached':False,
      'startup_mag_sensor_packet_attached':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
