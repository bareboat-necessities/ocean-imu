#!/usr/bin/env python3
"""Exact lifted coordinate for the finite Cayley/quaternion reset in P4.

Let c be the pre-update true-minus-estimated Cayley error, d=E_theta*K*y the
actual same-Joseph attitude correction, and a the Cayley vector of the deployed
correction quaternion.  Exact composition gives

    c+ = (c-a + .5 a x c)/(1 + .25 a'c).

Against the covariance-reset tangent target G(d)(c-d), define rho by

    c+ = G(d)(c-d) + rho.

A direct rearrangement, with x=a-d, w=d x c and h=a'c, gives the exact vector
identity

 (1+.25 h) rho = -x + .5 x x c
                       -.25 h c + .25 h d -.125 h w.       (R)

This producer makes every nonlinear product in (R) an explicit lifted
coordinate.  The production augmented LDLT can therefore retain reset direction
instead of replacing rho by a free norm ball.  The lifted graph still needs
quadratic product/IQC constraints; this module supplies the exact linear
incidence relation and algebraic regression only.

No K box is created here: d must be supplied downstream as E_theta*K*y from the
same P/H/R/S/K Joseph cell.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

QUALIFICATION='OU3_P4_EXACT_FINITE_RESET_JOINT_COORDINATE_V1'


def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def add(*vs):return [sum(v[i] for v in vs) for i in range(3)]
def scale(s,v):return [s*x for x in v]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def norm(v):return math.sqrt(sum(x*x for x in v))

def exact_reset(c,d,a):
    h=dot(a,c);den=1+.25*h
    if den<=0:raise ValueError('reset Cayley denominator nonpositive')
    cp=scale(1/den,add(c,scale(-1,a),scale(.5,cross(a,c))))
    target=add(c,scale(-1,d),scale(.5,cross(d,c)))
    rho=sub(cp,target)
    return cp,rho

def lifted(c,d,a,rho):
    x=sub(a,d);xc=cross(x,c);w=cross(d,c);h=dot(a,c)
    hc=scale(h,c);hd=scale(h,d);hw=scale(h,w);hr=scale(h,rho)
    # rho + .25 h rho + x - .5 xxc + .25 hc - .25 hd + .125 hw = 0
    residual=add(rho,scale(.25,hr),x,scale(-.5,xc),scale(.25,hc),scale(-.25,hd),scale(.125,hw))
    return {'x_a_minus_d':x,'x_cross_c':xc,'w_d_cross_c':w,'h_a_dot_c':h,'h_c':hc,'h_d':hd,'h_w':hw,'h_rho':hr,'linear_incidence_residual':residual}

def build():
    # Non-axis, non-collinear regression point.  a is parallel to d as in the
    # deployed correction but has a nonlinear Cayley magnitude.
    c=[.31,-.17,.22];d=[.18,-.09,.12];r=norm(d);ratio=(2*math.tan(r/2)/r) if r else 1.;a=scale(ratio,d)
    cp,rho=exact_reset(c,d,a);z=lifted(c,d,a,rho);defect=norm(z['linear_incidence_residual'])
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','same_cell_attitude_correction_required':'d=E_theta*K*y','independent_correction_port_forbidden':True,'deployed_correction_cayley_parallel_to_d':True,'exact_reset_equation':'(1+.25*aTc)rho=-(a-d)+.5(a-d)xc-.25(aTc)c+.25(aTc)d-.125(aTc)(dxc)','lifted_coordinates':['c','d','a','rho','x=a-d','xc=x_cross_c','w=d_cross_c','h=a_dot_c','hc','hd','hw','hrho'],'exact_linear_incidence_in_lifted_coordinate_available':True,'regression_linear_incidence_residual_norm':defect,'regression_exact_reset_output':cp,'regression_reset_defect':rho,'free_reset_norm_ball_used':False,'packet_count_multiplier_used':False,'source_uniform_augmented_LDLT_closed_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('independent_correction_port_forbidden','deployed_correction_cayley_parallel_to_d','exact_linear_incidence_in_lifted_coordinate_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('free_reset_norm_ball_used','packet_count_multiplier_used','source_uniform_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if not float(d.get('regression_linear_incidence_residual_norm',math.inf))<1e-14:f.append('exact reset lifted identity regression failed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'identity_defect':d['regression_linear_incidence_residual_norm'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
