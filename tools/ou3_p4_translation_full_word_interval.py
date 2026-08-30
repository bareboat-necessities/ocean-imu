#!/usr/bin/env python3
"""Validated complete-word translation dissipation on the old limiting P3 cell.

The old P4 radius is destroyed by a translation margin obtained from one
covariance seed.  Here we accumulate a certified covariance component through a
whole one-second word instead of discarding every later process contribution.

The proof is done in the same dimensionless coordinates already certified by
P3,

 D=diag(sigma_min*h, sigma_min*h^2, sigma_min*h^3, sigma_min),
 y=D^{-1}[v,p,S,a_w].

In these coordinates every prediction injects at least rho_trans I.  This avoids
both the physical-unit condition number and the tiny h^6 S scale that made a
naive interval propagation numerically useless.

For a lower covariance component L, each correction is conservatively replaced
by the *most informative* admissible scalar observation of its translation
support.  The operator

  L -> (L^{-1}+H'R^{-1}H)^{-1}

is monotone in L, and the Joseph covariance produced by any implemented gain is
no smaller than this optimal posterior.  Therefore pretending that S=0 and the
accelerometer are accepted every IMU sample can only make the propagated lower
component smaller.  Rejected/not-due shipping branches are automatically
covered.  Magnetometer has no translation columns.

All promoted arithmetic is outward-rounded Interval arithmetic.  Tau is
adaptively subdivided if interval dependency prevents a positive Loewner lower
matrix.  Each tau subcell is certified separately and the minimum endpoint
margin is reported.  This file proves only the previously limiting translation
source cell; it deliberately does not promote complete P4 yet.
"""
from __future__ import annotations

import argparse, json, math, re
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_validated_transcendentals as VT
import ou3_p4_worst_translation_cell as WORST

REPO=Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
HORIZON_S=1.0
MAX_TAU_SPLIT_DEPTH=14


def down(x): return math.nextafter(float(x),-math.inf)
def up(x): return math.nextafter(float(x),math.inf)
def I(x): return Interval.outward_bounds(float(x),float(x))
def _point_matrix(A): return [[I(float(x)) for x in row] for row in A]


def _center_radius_lower(A):
    """Point L such that every symmetric matrix in interval box A is >= L."""
    A=matrix_symmetric_hull(A); n=len(A); C=[[0.0]*n for _ in range(n)]; rad_rows=[0.0]*n
    for i in range(n):
        for j in range(n):
            a=A[i][j]; c=0.5*a.lo+0.5*a.hi
            c=min(max(c,a.lo),a.hi); C[i][j]=c
            r=up(max(abs(a.lo-c),abs(a.hi-c))); rad_rows[i]=up(rad_rows[i]+r)
    rad=max(rad_rows); L=[row[:] for row in C]
    for i in range(n): L[i][i]=down(L[i][i]-rad)
    ok,_=symmetric_positive_definite_ldlt(_point_matrix(L))
    if not ok: raise RuntimeError(f'Loewner lower extraction lost SPD (conditioned radius={rad:.3e})')
    return L,rad


def _transition_scaled(tau:Interval,h:float):
    x=I(h)/tau
    if x.hi>=1.0e-2: raise RuntimeError('limiting cell unexpectedly left shipping small-x branch')
    alpha=VT.exp_interval(-x); em1=VT.expm1_interval(-x)
    phi_va=-(tau*em1)
    x2=x*x; x3=x2*x; x4=x3*x; x5=x4*x
    phi_pa=tau*tau*(I(0.5)*x2-I(1/6)*x3+I(1/24)*x4)
    phi_Sa=tau*tau*tau*(I(1/6)*x3-I(1/24)*x4+I(1/120)*x5)
    z=Interval.point(0.0); o=Interval.point(1.0)
    # D^-1 F D; sigma cancels exactly.
    return [[o,z,z,phi_va/I(h)],
            [o,o,z,phi_pa/I(h*h)],
            [I(0.5),o,o,phi_Sa/I(h*h*h)],
            [z,z,z,alpha]]


def _predict(L,F,rho):
    M=matrix_mul(matrix_mul(F,_point_matrix(L)),matrix_transpose(F))
    for i in range(4): M[i][i]=M[i][i]+I(rho)
    return _center_radius_lower(M)


def _measure(L,idx,rnorm):
    P=_point_matrix(L); den=P[idx][idx]+I(rnorm)
    if den.lo<=0: raise RuntimeError('nonpositive comparison innovation variance')
    M=[[P[i][j]-(P[i][idx]*P[idx][j])/den for j in range(4)] for i in range(4)]
    return _center_radius_lower(M)


def _spd_delta(L,upper,delta):
    A=_point_matrix(L); q=I(delta)
    for i in range(4): A[i][i]=A[i][i]-q*I(upper[i])
    return symmetric_positive_definite_ldlt(matrix_symmetric_hull(A))[0]


def _delta(L,upper):
    trial=1e-18; lo=0.0
    while trial<1.0 and _spd_delta(L,upper,trial): lo=trial; trial*=10.0
    if lo==0.0:
        trial=1e-36
        while trial<1e-18 and not _spd_delta(L,upper,trial): trial*=10.0
        if trial>=1e-18: return 0.0
        lo=trial; trial*=10.0
    hi=min(1.0,trial)
    for _ in range(64):
        mid=math.sqrt(lo*hi)
        if _spd_delta(L,upper,mid): lo=mid
        else: hi=mid
    return down(lo)


def _member(text,name):
    m=re.search(rf'float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;',text)
    if not m: raise RuntimeError(f'cannot source-bind {name}')
    return float(m.group(1))


def _propagate_tau(tau,h,n,rho,rS,rA,depth=0):
    try:
        F=_transition_scaled(tau,h); L=[[0.0]*4 for _ in range(4)]; maxrad=0.0
        for k in range(n):
            if k==0:
                L=[[rho if i==j else 0.0 for j in range(4)] for i in range(4)]
            else:
                L,r=_predict(L,F,rho); maxrad=max(maxrad,r)
            L,r=_measure(L,2,rS); maxrad=max(maxrad,r)
            L,r=_measure(L,3,rA); maxrad=max(maxrad,r)
        return [(tau,L,maxrad,depth)]
    except RuntimeError as e:
        if depth>=MAX_TAU_SPLIT_DEPTH: raise
        mid=math.sqrt(tau.lo*tau.hi)
        a=Interval.outward_bounds(tau.lo,mid); b=Interval.outward_bounds(mid,tau.hi)
        return _propagate_tau(a,h,n,rho,rS,rA,depth+1)+_propagate_tau(b,h,n,rho,rS,rA,depth+1)


def _mode(mode,p):
    c=WORST.build_cell(mode,p); s=WORST.serializable(c); row=c['row']; h=float(c['sched']['dt_s'])
    x=c['x']; tau=Interval.outward_bounds(h/x.hi,h/x.lo); sigma=float(c['sigma'].lo); rho=float(c['rho_translation_lower'])
    scale2=[(sigma*h)**2,(sigma*h*h)**2,(sigma*h*h*h)**2,sigma*sigma]
    u=list(map(float,row['Sigma_diagonal_upper'])); physical=[u[6],u[9],u[12],u[15]]
    upper=[up(physical[i]/scale2[i]) for i in range(4)]
    text=WRAPPER.read_text(encoding='utf-8'); rh=min(_member(text,'R_S_x_factor_'),_member(text,'R_S_y_factor_'),1.0)
    # Smallest R gives the artificial maximum-information comparison.
    rs=(I(rh)*I(float(c['rs'].lo))).lo; rs_var=I(rs).square().lo
    acc=float(c['vector']['configured_measurement_bounds']['acc_measurement_std_mps2']); acc_var=I(acc).square().lo
    # Normalize measurement variances by the corresponding D coordinate.
    rS=(I(rs_var)/I(scale2[2])).lo; rA=(I(acc_var)/I(scale2[3])).lo
    n=int(math.ceil(HORIZON_S/h)); leaves=_propagate_tau(tau,h,n,rho,rS,rA)
    certified=[]
    for t,L,rad,depth in leaves:
        d=_delta(L,upper)
        if d<=0 or not _spd_delta(L,upper,d): raise RuntimeError(f'nonpositive endpoint margin on tau leaf {t.as_list()}')
        certified.append({'tau_s':t.as_list(),'delta_lower':d,'max_conditioned_radius_removed':rad,'split_depth':depth})
    worst=min(certified,key=lambda q:q['delta_lower']); old=float(row['direct_translation_generalized_margin_lower'])
    return {'source_cell':s,'conditioned_coordinates':'D^-1[v,p,S,a_w], D=diag(sigma*h,sigma*h^2,sigma*h^3,sigma)','tau_interval_s':tau.as_list(),'tau_leaf_count':len(certified),'max_tau_split_depth_used':max(q['split_depth'] for q in certified),'steps':n,'process_injection_lower_conditioned':rho,'artificial_S_variance_conditioned':rS,'artificial_acc_aw_variance_conditioned':rA,'translation_covariance_upper_conditioned':upper,'tau_leaf_certificates':certified,'complete_word_translation_margin_lower':worst['delta_lower'],'limiting_tau_leaf':worst,'old_single_seed_translation_margin_lower':old,'margin_widening_factor_lower':down(worst['delta_lower']/old),'interval_ldlt_endpoint_recertified':True}


def build(domain_path=DEFAULT_DOMAIN):
    p=Path(domain_path).resolve(); modes={}; f=[]
    for mode in ('H','A'):
        try: modes[mode]=_mode(mode,p)
        except Exception as e: f.append(f'{mode}: {e}')
    return {'qualification':'OU3_P4_VALIDATED_WORST_CELL_COMPLETE_WORD_TRANSLATION_DISSIPATION','source_only':True,'trajectory_replay_used':False,'outward_rounded':True,'horizon_s':HORIZON_S,'modes':modes,'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS':'PASS' if not f and len(modes)==2 else 'NOT_ESTABLISHED','P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED','remaining_obligation':'extend conditioned complete-word propagation to every reachable source cell/edge and attitude-bias blocks, then validate exact nonlinear return map','failures':f}


def validate(d):
    f=list(d.get('failures',[]))
    if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False or d.get('outward_rounded') is not True: f.append('qualification flags invalid')
    for mode in ('H','A'):
        m=d.get('modes',{}).get(mode,{})
        if not float(m.get('complete_word_translation_margin_lower',0))>0: f.append(f'{mode}: no complete-word translation margin')
        if m.get('interval_ldlt_endpoint_recertified') is not True: f.append(f'{mode}: endpoint not recertified')
        if not float(m.get('margin_widening_factor_lower',0))>1: f.append(f'{mode}: complete word did not widen seed')
    if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED': f.append('partial result prematurely promoted P4')
    return f


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'translation_status':d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],'modes':{m:{'delta':d.get('modes',{}).get(m,{}).get('complete_word_translation_margin_lower'),'factor':d.get('modes',{}).get(m,{}).get('margin_widening_factor_lower'),'tau_leaves':d.get('modes',{}).get(m,{}).get('tau_leaf_count')} for m in ('H','A')},'failures':f},indent=2,sort_keys=True)); return 0 if not f else 2

if __name__=='__main__': raise SystemExit(main())
