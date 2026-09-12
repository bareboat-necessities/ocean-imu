"""Same-operand finite-real Live tilt-watchdog and preserve-yaw reset graph.

This removes the two free operands that remained in ``finite_live_tilt_prefix``:
the watchdog angle and the final preserve-yaw reset quaternion.  The watchdog
angle witness is tied to the exact post-accelerometer nominal attitude.  When
the watchdog fires, the reset witness is derived from that same predecessor
attitude and the exact guarded accelerometer sample consumed by the current IMU
prefix.

The ideal-real preserve-yaw construction is algebraic.  It uses the half-angle
identities corresponding to shipping's atan2/asin/AngleAxis calls, rather than
accepting their final quaternion as an independent value.  Binary32/libm and
Eigen normalization correspondence remain deployment obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH

R = P.rational


# Shipping compares acos(cos_tilt)*57.295779513f > 70.  The float literal is
# represented here by the exact decimal rational used by the finite-real layer.
DEG_PER_RAD_SHIPPING = F(57295779513, 1000000000)
THRESHOLD_RAD = F(70) / DEG_PER_RAD_SHIPPING


def _cos_interval_alternating(x, terms=14):
    """Rigorous rational cosine enclosure for |x|<sqrt(2)."""
    x=abs(R(x)); x2=x*x
    term=F(1); total=term
    partials=[total]
    for k in range(1,terms+1):
        term = -term*x2/F((2*k-1)*(2*k))
        total += term; partials.append(total)
    a,b=partials[-2],partials[-1]
    return (min(a,b),max(a,b))


COS70_LO, COS70_HI = _cos_interval_alternating(THRESHOLD_RAD)


def watchdog_over_limit(state:CORE.State):
    """Resolve the shipping >70deg branch from the exact current attitude.

    The only unresolved case is the tiny rational enclosure containing the
    transcendental threshold itself; it remains fail-closed for deployment
    libm correspondence rather than accepting a free tilt scalar.
    """
    if not isinstance(state,CORE.State): raise TypeError('finite core required')
    qbw=P.quat_conj(state.q_hat)
    cos_tilt=max(F(-1),min(F(1),CORE.rotation(qbw)[2][2]))
    if cos_tilt < COS70_LO: return True
    if cos_tilt >= COS70_HI: return False
    raise NotImplementedError('watchdog threshold lies inside rigorous cos(70/K) enclosure')


@dataclass(frozen=True)
class SignedHalfWitness:
    """Half angle for an angle whose sine and nonnegative cosine are known."""
    sine: F
    cosine: F
    cos_half: F
    sin_half: F
    def __post_init__(self):
        s,c,ch,sh = map(R,(self.sine,self.cosine,self.cos_half,self.sin_half))
        if c < 0 or s*s + c*c != 1:
            raise ValueError('unit sine/cosine pair with nonnegative cosine required')
        if ch < 0 or ch*ch + sh*sh != 1 or ch*ch-sh*sh != c or 2*ch*sh != s:
            raise ValueError('signed half-angle witness detached from same sine/cosine')
        object.__setattr__(self,'sine',s); object.__setattr__(self,'cosine',c)
        object.__setattr__(self,'cos_half',ch); object.__setattr__(self,'sin_half',sh)


@dataclass(frozen=True)
class AccTiltWitness:
    """Square-root witnesses for shipping ``quaternion_from_acc`` general branch."""
    acc_norm: TILT.SqrtWitness
    axis_norm: TILT.SqrtWitness | None = None
    cos_half: F | None = None
    sin_half: F | None = None


def _normalize_quaternion(q, norm:TILT.SqrtWitness):
    q=tuple(M.vec(q,4)); q2=M.dot(q,q)
    if not isinstance(norm,TILT.SqrtWitness) or norm.radicand!=q2 or norm.value<=F(1,10**8):
        raise ValueError('quaternion normalization witness detached/tiny')
    qn=tuple(x/norm.value for x in q)
    if M.dot(qn,qn)!=1: raise AssertionError('normalized quaternion lost unit norm')
    return qn


def _qref_from_acc(acc_body, witness:AccTiltWitness):
    acc=tuple(M.vec(acc_body,3)); a2=M.dot(acc,acc)
    if not isinstance(witness,AccTiltWitness) or not isinstance(witness.acc_norm,TILT.SqrtWitness):
        raise TypeError('accelerometer tilt witnesses required')
    if witness.acc_norm.radicand!=a2 or witness.acc_norm.value<=F(1,10**8):
        raise ValueError('accelerometer norm witness detached/tiny')
    an=tuple(x/witness.acc_norm.value for x in acc)
    target=tuple(-x for x in an)
    c=target[2]
    axis=(-target[1],target[0],F(0))
    aaxis2=M.dot(axis,axis)
    if aaxis2==0:
        if witness.axis_norm is not None or witness.cos_half is not None or witness.sin_half is not None:
            raise ValueError('parallel accelerometer branch consumes no axis/half-angle witnesses')
        return (F(1),F(0),F(0),F(0)) if c>0 else (F(0),F(1),F(0),F(0))
    if not isinstance(witness.axis_norm,TILT.SqrtWitness) or witness.axis_norm.radicand!=aaxis2 or witness.axis_norm.value<=0:
        raise ValueError('accelerometer rotation-axis norm witness detached')
    ch,sh=R(witness.cos_half),R(witness.sin_half)
    if ch<0 or sh<0 or ch*ch+sh*sh!=1 or ch*ch-sh*sh!=c:
        raise ValueError('accelerometer half-angle witness detached from same gravity direction')
    u=tuple(x/witness.axis_norm.value for x in axis)
    q=(ch,sh*u[0],sh*u[1],sh*u[2])
    if M.dot(q,q)!=1: raise AssertionError('accelerometer attitude quaternion lost unit norm')
    return q


def preserve_yaw_witness(state:CORE.State, sample:SENSOR.GuardedImuSample, *,
                         old_q_norm:TILT.SqrtWitness, old_yaw_half:TILT.YawHalfWitness|None,
                         acc_tilt:AccTiltWitness,
                         pitch_cos:TILT.SqrtWitness|None=None,
                         pitch_half:SignedHalfWitness|None=None,
                         roll_half:TILT.YawHalfWitness|None=None):
    """Derive the ideal-real output of shipping initialize_from_acc_preserve_yaw."""
    if not isinstance(state,CORE.State) or not isinstance(sample,SENSOR.GuardedImuSample):
        raise TypeError('finite core and same guarded sample required')
    if sample.physical.history_id!=state.reference.history_id:
        raise ValueError('preserve-yaw reset sample detached from physical history')

    old_bw=_normalize_quaternion(P.quat_conj(state.q_hat),old_q_norm)
    ow,ox,oy,oz=old_bw
    sy=2*(ow*oz+ox*oy); cy=1-2*(oy*oy+oz*oz)
    if sy==0 and cy==0:
        if old_yaw_half is not None: raise ValueError('degenerate old yaw consumes no atan2 witness')
        qyaw=(F(1),F(0),F(0),F(0))
    else:
        if not isinstance(old_yaw_half,TILT.YawHalfWitness) or (old_yaw_half.c,old_yaw_half.s)!=(cy,sy):
            raise ValueError('old-yaw witness detached from predecessor attitude')
        qyaw=(old_yaw_half.cos_half,F(0),F(0),old_yaw_half.sin_half)

    qref_tilt=_qref_from_acc(sample.conditioned_accel_body,acc_tilt)
    tilt_bw=P.quat_conj(qref_tilt)
    tw,tx,ty,tz=tilt_bw
    sinp=max(F(-1),min(F(1),2*(tw*ty-tz*tx)))
    cp2=1-sinp*sinp
    if not isinstance(pitch_cos,TILT.SqrtWitness) or pitch_cos.radicand!=cp2:
        raise ValueError('pitch cosine sqrt detached from accel-only tilt')
    cp=pitch_cos.value
    if not isinstance(pitch_half,SignedHalfWitness) or (pitch_half.sine,pitch_half.cosine)!=(sinp,cp):
        raise ValueError('pitch half-angle detached from accel-only tilt')
    qpitch=(pitch_half.cos_half,F(0),pitch_half.sin_half,F(0))

    sr=2*(tw*tx+ty*tz); cr=1-2*(tx*tx+ty*ty)
    if sr==0 and cr==0:
        if roll_half is not None: raise ValueError('degenerate roll consumes no atan2 witness')
        qroll=(F(1),F(0),F(0),F(0))
    else:
        if not isinstance(roll_half,TILT.YawHalfWitness) or (roll_half.c,roll_half.s)!=(cr,sr):
            raise ValueError('roll witness detached from accel-only tilt')
        qroll=(roll_half.cos_half,roll_half.sin_half,F(0),F(0))

    qnew_bw=tuple(P.quat_mul(P.quat_mul(qyaw,qpitch),qroll))
    if M.dot(qnew_bw,qnew_bw)!=1:
        raise AssertionError('preserve-yaw Euler composition lost unit norm')
    qnew_hat=tuple(P.quat_conj(qnew_bw))
    down=tuple(M.mv(CORE.rotation(qnew_hat),(F(0),F(0),F(1))))
    if M.dot(down,down)!=1: raise AssertionError('world-down body axis lost unit norm')
    return WATCH.PreserveYawWitness(qnew_hat,down)


def readiness():
    return {
      'watchdog_threshold_tied_to_post_accelerometer_attitude':True,
      'watchdog_threshold_rigorous_cosine_enclosure':True,
      'preserve_yaw_old_heading_tied_to_same_pre_reset_attitude':True,
      'preserve_yaw_accel_tilt_tied_to_same_guarded_accelerometer':True,
      'preserve_yaw_pitch_roll_reconstruction_algebraically_bound':True,
      'reset_down_axis_derived_from_same_final_quaternion':True,
      'free_final_preserve_yaw_quaternion_removed_from_theorem_entry':True,
      'watchdog_boundary_sliver_and_binary32_libm_closed':False,
      'sqrt_atan2_asin_angleaxis_binary32_libm_closed':False,
      'tiny_parallel_branch_binary32_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
