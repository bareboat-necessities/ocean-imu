"""All-time real-arithmetic guard inactivity for the exact planar MOVING record.

A source/guard component theorem, not magnetic admission or float32 transfer.
Uses the literal seeded two-stage detector, RMS EW update and weight slew.
"""
from __future__ import annotations
from fractions import Fraction as F
import math


def certificate():
    g=F(196133,20000);nu_hi=F(11,35);dt=F(1,200)
    amplitude=F(1,50);acc=amplitude*nu_hi**2;jerk=amplitude*nu_hi**3
    omega=amplitude*nu_hi
    # |d/dt Ry(-psi)(p'' e_x-g cos(theta)e_z)| <= |p'''|+|psi'|(|p''|+g).
    force_lipschitz=jerk+omega*(g+acc)
    sample_delta=dt*force_lipschitz
    # pi>3.14159, exp(x)>=sum_{k=0}^8 x^k/k!, hence exp(-pi/4)<.46.
    exponent_lo=F(314159,100000)/4
    exp_lower=sum((exponent_lo**k/F(math.factorial(k)) for k in range(9)),F(0))
    gamma=F(23,50)
    assert exp_lower>1/gamma
    hp1=gamma*sample_delta/(1-gamma)
    hp2=2*hp1
    engage_lo=F(3,100)
    assert hp2<engage_lo
    return {
        'qualification':'OU3_PLANAR_GUARD_INACTIVITY_V1',
        'result_type':'PROVED analytical theorem in real arithmetic; exact rational margins',
        'scope':'specified planar source at h=1/200 s; seeded literal guard at default detector/engagement settings',
        'source_frequency_upper_rad_s':str(nu_hi),
        'force_derivative_norm_upper':str(force_lipschitz),
        'successive_sample_force_difference_upper':str(sample_delta),
        'detector_decay_upper':str(gamma),
        'exp_lower_polynomial':str(exp_lower),
        'first_detector_output_norm_upper':str(hp1),
        'second_detector_output_norm_upper':str(hp2),
        'detector_RMS_upper_mps2':float(hp2),
        'engage_floor_mps2':str(engage_lo),
        'strict_RMS_margin':str(engage_lo-hp2),
        'proof_induction':[
            'h1_k=gamma(h1_(k-1)+a_k-a_(k-1)), with h1_0=0',
            'second detector LP is a convex average of h1, hence ||h2||<=2 sup ||h1||',
            'sum of per-axis squared-output EWs is a convex average of ||h2||^2',
            'RMS<engage_lo gives target=0, so seeded weight=0 and conditioned acceleration equals raw acceleration forever'],
        'structures_preserved':['same continuous source','seeded two-stage detector','per-axis RMS EW sum','literal target and weight slew'],
        'relaxations_introduced':['source derivative and detector norm upper bounds only; no independently selected filter state'],
        'physical_or_shipping_constants_changed':False,
        'float32_guard_transfer_verified':False,
        'all_time_magnetic_service_verified':False,'theorem_closed':False,
    }

if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
