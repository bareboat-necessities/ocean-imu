"""Exact conditional planar private-Mahony tube; not a shipping certificate.

The common quadratic tube retains pitch error and integral feedback, and uses
an effective feedback interval compatible with a non-unit raw quaternion.
Binding the stated one-step error caps to every literal float operation, and
binding the complete frontend/tuner/clock state, remain OPEN.
"""
from __future__ import annotations
from fractions import Fraction as F
from .matrix_certificates import matmul,transpose,add,ldlt,encoded
from .planar_service_frontend_binding import dyadic,round32


def certificate():
    h=dyadic(round32(F(1,200)));kp=F(1,5);ki=F(1,50)
    G=[[F(10),F(-5)],[F(-5),F(15)]]
    margin=F(27,100000);radius=F(3,500);a_ends=(F(49,100),F(51,100))
    out=[]
    for a in a_ends:
        A=[[1-a*(kp*h+ki*h*h),h/10],[-10*ki*h*a,F(1)]]
        B=[[a*(kp*h+ki*h*h)],[10*ki*h*a]]
        gap=add(add(G,matmul(matmul(transpose(A),G),A),F(-1)),G,-margin)
        _,pivots=ldlt(gap)
        assert min(pivots)>0
        bgb=matmul(matmul(transpose(B),G),B)[0][0]
        assert bgb<=F(1,500)**2
        out.append({'effective_gain':str(a),'A':encoded(A),'B':encoded(B),
            'G_minus_ATGA_minus_marginG':encoded(gap),
            'positive_LDL_pivots':[str(x) for x in pivots],'B_G_norm_upper':'1/500'})
    # ||B d||_G <= .002*.00025. Remaining caps include sample curvature,
    # arctan discretization and the stated (NOT YET source-bound) float errors.
    d=F(1,4000);e_force=F(46,10**9);i_force=F(2,10**9)
    charge=F(1,500)*d+F(19,6)*e_force+F(31,8)*i_force
    retained_margin=radius*margin/2
    assert charge<retained_margin
    return {
        'qualification':'OU3_PLANAR_MAHONY_CONDITIONAL_TUBE_V1',
        'result_type':'CONDITIONAL; common-quadratic and charge inequalities proved exactly',
        'coordinates':'z=(private normalized-direction pitch minus physical psi, 10*integralFBy)',
        'G':encoded(G),'h_binary32_exact':str(h),'squared_norm_decrement_lower':str(margin),
        'root_G_norm_radius':str(radius),'endpoint_checks':out,
        'interior_gain_argument':'A is affine in a, so A^T G A is matrix convex; endpoint inequalities cover [0.49,0.51]',
        'exogenous_angle_offset_cap':str(d),'first_coordinate_additive_cap':str(e_force),
        'scaled_integral_additive_cap':str(i_force),'G_norm_charge_upper':str(charge),
        'one_step_radial_margin_lower':str(retained_margin),
        'strict_radial_margin':str(retained_margin-charge),
        'pitch_error_abs_bound':'.0021','integralFBy_abs_bound':'.00017',
        'raw_quaternion_norm_not_assumed_one':True,
        'required_raw_squared_norm_interval':['.995','1.001'],
        'complete_literal_one_step_error_caps_verified':False,
        'literal_initialization_in_tube_verified':False,
        'all_time_frontend_tube_verified':False,
        'all_time_magnetic_service_verified':False,'theorem_closed':False,
    }

if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
