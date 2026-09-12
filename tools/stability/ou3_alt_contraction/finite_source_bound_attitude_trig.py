"""Bind finite attitude trig witnesses to the actual source-owned rotation angle.

``finite_attitude_runtime.TrigWitness`` deliberately carries finite values for
shipping's sin/cos/inverse-rate evaluations, but its component contract only
checks unit-circle and inverse-rate algebra.  That is not enough for a theorem:
a detached point on the unit circle could otherwise be supplied at each step.

This layer proves a necessary real-arithmetic angle relation.  Because the
component witness already enforces ``inv_omega^2 * ||w||^2 = 1``, the exact
represented angle is ``theta = t / inv_omega``.  For 0 <= theta <= 1 rad the
alternating Taylor series gives rigorous rational enclosures

  theta-theta^3/6 <= sin(theta) <= theta-theta^3/6+theta^5/120
  1-theta^2/2     <= cos(theta) <= 1-theta^2/2+theta^4/24.

The theorem-facing source word requires each supplied full/half trig result to
lie in those enclosures for the SAME angular rate and time.  Values outside the
one-radian local interval fail closed; proving that every admitted source/error
prefix stays inside it remains a separate retention obligation.  Binary32
sqrt/div/sin/cos correspondence is also still open.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT

MAX_THETA = F(1)


def trig_enclosure(theta):
    """Rigorous alternating-series enclosure for 0 <= theta <= 1 rad."""
    x=P.rational(theta)
    if x < 0 or x > MAX_THETA:
        raise ValueError('attitude trig angle outside certified <=1 rad enclosure')
    x2=x*x; x3=x2*x; x4=x2*x2; x5=x4*x
    sin_lo=x-x3/F(6)
    sin_hi=sin_lo+x5/F(120)
    cos_lo=F(1)-x2/F(2)
    cos_hi=cos_lo+x4/F(24)
    return sin_lo,sin_hi,cos_lo,cos_hi


def validate_witness(w, witness:ATT.TrigWitness):
    """Require one finite trig witness to belong to its actual-angle enclosure."""
    if not isinstance(witness,ATT.TrigWitness):
        raise TypeError('finite attitude trig witness required')
    v=tuple(P.vec(w,3)); w2=M.dot(v,v)
    if w2 <= 0:
        raise ValueError('general trig branch requires nonzero angular rate')
    # The lower component already requires inv^2*w2==1. Repeat here so this
    # theorem-facing lemma has no hidden dependency on constructor ordering.
    if witness.inv_omega*witness.inv_omega*w2 != 1:
        raise ValueError('inverse-rate witness detached from same angular rate')
    theta=witness.t/witness.inv_omega
    slo,shi,clo,chi=trig_enclosure(theta)
    if not slo <= witness.sin_theta <= shi:
        raise ValueError('sine witness detached from source-owned rotation angle')
    if not clo <= witness.cos_theta <= chi:
        raise ValueError('cosine witness detached from source-owned rotation angle')
    return theta


def validate(runtime:ATT.AngularRuntime):
    """Bind both full and half shipping trig calls to one AngularRuntime."""
    if not isinstance(runtime,ATT.AngularRuntime):
        raise TypeError('AngularRuntime required')
    if runtime.small_rate:
        if runtime.full is not None or runtime.half is not None:
            raise ValueError('small-rate branch consumes no trig witnesses')
        return runtime
    full_theta=validate_witness(runtime.w,runtime.full)
    half_theta=validate_witness(runtime.w,runtime.half)
    if full_theta != 2*half_theta:
        raise ValueError('full/half trig angles detached from same source step')
    return runtime


def readiness():
    return {
      'trig_full_half_angles_derived_from_same_angular_rate_and_step':True,
      'sin_cos_values_rigorously_enclosed_at_source_owned_angles':True,
      'detached_unit_circle_points_rejected':True,
      'one_radian_local_angle_guard_explicit':True,
      'one_radian_guard_retained_for_every_admitted_prefix':False,
      'binary32_sqrt_div_sin_cos_correspondence_closed':False,
      'ALT_LIVE_PASS':False,
    }
