"""Bind attitude trig witnesses to the actual source-owned rotation angle.

``finite_attitude_runtime.TrigWitness`` carries finite values for shipping's
sin/cos/inverse-rate evaluations. Unit-circle algebra alone is insufficient: a
detached point on the circle could otherwise be supplied at each step.

This layer now uses the global Taylor theorem rather than a <=1 rad alternating
series. For every finite rational ``theta`` it chooses a finite Taylor degree
such that the Lagrange remainder (all sin/cos derivatives have magnitude <=1)
is at most ``TRIG_TOL``. No range reduction, pi constant, estimator-error bound,
or one-radian retention assumption is needed.

The component witness already enforces ``inv_omega^2 ||w||^2 = 1``, hence the
represented source-owned angle is exactly ``theta=t/inv_omega``. Full and half
calls must lie in rigorous rational enclosures at those same angles. This
closes exact-real trig ancestry globally; binary32 sqrt/div/libm sin/cos
correspondence remains a separate deployment obligation.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import factorial

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT

TRIG_TOL = F(1,10**8)


def _remainder(abs_x:F, degree:int)->F:
    return abs_x**(degree+1)/factorial(degree+1)


def _partial_sin(x:F,m:int)->F:
    term=x; total=term
    for k in range(1,m+1):
        term *= -x*x/F((2*k)*(2*k+1))
        total += term
    return total


def _partial_cos(x:F,m:int)->F:
    term=F(1); total=term
    for k in range(1,m+1):
        term *= -x*x/F((2*k-1)*(2*k))
        total += term
    return total


def trig_enclosure(theta, *, tol=TRIG_TOL):
    """Global rigorous rational sin/cos enclosure for any finite rational angle.

    After degree n, Taylor's theorem gives |R_n(x)| <= |x|^(n+1)/(n+1)!
    because every derivative of sin/cos is bounded by one. Factorial growth
    therefore guarantees a finite degree for every finite x and positive tol.
    """
    x=P.rational(theta); eps=P.rational(tol)
    if x < 0:
        raise ValueError('attitude rotation-angle magnitude must be nonnegative')
    if eps <= 0:
        raise ValueError('positive rigorous trig enclosure tolerance required')
    if x == 0:
        return F(0),F(0),F(1),F(1)
    a=abs(x); m=0
    while True:
        # sin polynomial has degree 2m+1; cos polynomial degree 2m.
        rs=_remainder(a,2*m+1)
        rc=_remainder(a,2*m)
        if rs<=eps and rc<=eps:
            break
        m += 1
    s=_partial_sin(x,m); c=_partial_cos(x,m)
    return s-rs,s+rs,c-rc,c+rc


def validate_witness(w, witness:ATT.TrigWitness):
    """Require one finite trig witness to belong to its actual-angle enclosure."""
    if not isinstance(witness,ATT.TrigWitness):
        raise TypeError('finite attitude trig witness required')
    v=tuple(P.vec(w,3)); w2=M.dot(v,v)
    if w2 <= 0:
        raise ValueError('general trig branch requires nonzero angular rate')
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
      'global_finite_rational_angle_enclosure_available':True,
      'one_radian_local_angle_guard_required':False,
      'one_radian_guard_retained_for_every_admitted_prefix':True,
      'binary32_sqrt_div_sin_cos_correspondence_closed':False,
      'ALT_LIVE_PASS':False,
    }
