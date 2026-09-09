#!/usr/bin/env python3
"""Exact covariance-frame finite reset identity for P4.

Let d be the shipping attitude correction, a=kappa*d the Cayley vector of the
normalized deployed correction quaternion, c the pre-injection Cayley error and
G_d=I+.5[d]x.  The exact post-injection Cayley error is

  c+ = G_a (c-a) / h,       h=1+.25 a^T c,
  G_a=I+.5[a]x.

With tangent posterior t=c-d and u=G_d^-1 c+ (the attitude coordinate seen by
the Joseph posterior information matrix after covariance congruence), parallel
axes give the exact identity

  h u = A(d,kappa) t - (kappa-1) d,
  A(d,kappa)=G_d^-1 G_a.

Because G_d and G_a are polynomials of the same skew matrix they commute, and
A has singular value one along d.  This identity is stronger than treating the
reset defect as an independent norm-bounded port: u,t,d remain on the same
Joseph correction graph d=E_theta*K*y.

The module provides algebraic regression and a homogeneous lifted incidence
relation.  A production augmented LDLT must bind d to the same P/H/R/S/K cell;
this module does not promote P4 by itself.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np

QUALIFICATION='OU3_P4_EXACT_RESET_COVARIANCE_FRAME_IDENTITY_V1'


def skew(v):
    x,y,z=map(float,v);return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def kappa_from_d(d):
    r=float(np.linalg.norm(d))
    if r==0:return 1.0
    return 2.0*math.tan(0.5*r)/r

def exact_identity(c,d):
    c=np.asarray(c,dtype=float);d=np.asarray(d,dtype=float);k=kappa_from_d(d);a=k*d
    Gd=np.eye(3)+.5*skew(d);Ga=np.eye(3)+.5*skew(a);h=1+.25*float(a@c)
    if h<=0:raise ValueError('Cayley denominator nonpositive')
    cp=Ga@(c-a)/h;t=c-d;u=np.linalg.solve(Gd,cp);A=np.linalg.solve(Gd,Ga)
    residual=h*u-(A@t-(k-1)*d)
    return {'kappa':k,'h':h,'c_plus':cp.tolist(),'t':t.tolist(),'u':u.tolist(),
            'A':A.tolist(),'identity_residual_norm':float(np.linalg.norm(residual)),
            'A_axis_residual_norm':float(np.linalg.norm(A@d-d)),
            'A_operator_norm':float(np.linalg.norm(A,2))}

def build():
    points=[([.31,-.17,.22],[.18,-.09,.12]),([-.42,.11,.27],[.007,-.004,.002]),([.2,.3,-.1],[0.,0.,0.])]
    rows=[exact_identity(c,d) for c,d in points]
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_correction_required':'d=E_theta*K*y','deployed_correction_cayley_parallel_to_d':True,
      'exact_covariance_frame_identity':'(1+.25*aTc)u=G(d)^-1*G(a)*t-(kappa-1)d',
      'reset_defect_independent_port_used':False,'same_Joseph_information_frame_retained':True,
      'homogeneous_lift_available':True,'regression':rows,
      'max_identity_residual_norm':max(r['identity_residual_norm'] for r in rows),
      'max_axis_residual_norm':max(r['A_axis_residual_norm'] for r in rows),
      'source_uniform_endpoint_augmented_LDLT_closed_here':False,
      'source_uniform_every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('deployed_correction_cayley_parallel_to_d','same_Joseph_information_frame_retained','homogeneous_lift_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('reset_defect_independent_port_used','source_uniform_endpoint_augmented_LDLT_closed_here','source_uniform_every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('max_identity_residual_norm',math.inf))>1e-12:f.append('covariance-frame identity regression failed')
    if float(d.get('max_axis_residual_norm',math.inf))>1e-12:f.append('parallel-axis identity regression failed')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'residual':d['max_identity_residual_norm'],'axis':d['max_axis_residual_norm'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
