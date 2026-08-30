#!/usr/bin/env python3
"""Validated complete-word translation covariance accumulation for P4.

This closes the first quantitative piece of the replacement P4 route on the
source cell that destroys the old certificate.  It tracks a covariance
*component* supported on one [v,p,S,a_w] axis.  The component is legitimate in
the full MEKF because translation process noise is injected in this invariant
prediction subspace; translation does not feed the attitude/bias prediction.

For every prediction the existing P3 process proof supplies

    Q >= rho_Q D_h D_h^T,
    D_h=diag(sigma*h, sigma*h^2, sigma*h^3, sigma).

We accumulate that positive contribution through the exact interval family of
the shipping integrated-OU transition instead of throwing it away after one
sample.  At every sample we then deliberately apply MORE measurement
information than shipping is guaranteed to apply: one S=0 scalar observation
and one accelerometer observation of a_w.  This remains a lower covariance
comparison because for every implemented gain K

 (I-KH)L(I-KH)' + K R K' >= (L^-1 + H'R^-1H)^-1,

and a rejected/not-due correction is less informative than the artificial
accepted correction used here.  Magnetometer has no translation columns and is
irrelevant to this covariance component.  Isotropic accelerometer R makes its
a_w information exactly axis-independent despite body/world rotation.

All theorem arithmetic below is outward-rounded Interval arithmetic.  Wide
interval matrix families are converted after each operation to a point Loewner
lower matrix C-rI, where C is an enclosed midpoint and r is an outward upper
bound on the spectral norm of the interval-radius matrix.  Every reported
endpoint generalized margin is independently re-certified by interval LDLT.

This producer certifies the previous *worst translation source cell*.  It does
not yet promote the complete H/A P4 theorem: the same construction must be
extended to all reachable source cells/edges and the attitude/bias component,
then lifted through the exact nonlinear return map.
"""
from __future__ import annotations

import argparse, json, math, re
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_validated_transcendentals as VT
import ou3_p4_frontier_margin_diagnostic as DIAG
import ou3_source_reachable_matrix_p3 as P3BASE
import ou3_vector_uco_certificate as VECTOR

REPO=Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
HORIZON_S=1.0


def down(x): return math.nextafter(float(x),-math.inf)
def up(x): return math.nextafter(float(x),math.inf)
def I(x): return Interval.outward_bounds(float(x),float(x))

def _point_matrix(A): return [[I(float(x)) for x in row] for row in A]
def _zero(n): return [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
def _diag(v):
 A=_zero(len(v))
 for i,x in enumerate(v): A[i][i]=I(x)
 return A

def _center_radius_lower(A):
 """Return point L with A >= L for every symmetric A in the interval box."""
 A=matrix_symmetric_hull(A); n=len(A)
 C=[[0.0]*n for _ in range(n)]; R=[[0.0]*n for _ in range(n)]
 for i in range(n):
  for j in range(n):
   a=A[i][j]
   c=0.5*a.lo+0.5*a.hi
   if c<a.lo: c=a.lo
   if c>a.hi: c=a.hi
   C[i][j]=c
   R[i][j]=up(max(abs(a.lo-c),abs(a.hi-c)))
 # Symmetric radius matrix => ||E||_2 <= max row sum(|E|).
 rad=0.0
 for row in R:
  s=0.0
  for x in row: s=up(s+x)
  rad=max(rad,s)
 L=[[C[i][j] for j in range(n)] for i in range(n)]
 for i in range(n): L[i][i]=down(L[i][i]-rad)
 # Enclose the chosen binary64 lower matrix and independently verify SPD.
 Li=_point_matrix(L)
 ok,_=symmetric_positive_definite_ldlt(Li)
 if not ok:
  raise RuntimeError(f'Loewner lower extraction lost SPD (radius={rad:.3e})')
 return L,rad

def _transition_interval(tau: Interval,h: float):
 hi=I(h); x=hi/tau
 if x.hi>=1.0e-2: raise RuntimeError('worst-cell backend expected shipping small-x branch')
 alpha=VT.exp_interval(-x); em1=VT.expm1_interval(-x)
 phi_va=-(tau*em1)
 x2=x*x; x3=x2*x; x4=x3*x; x5=x4*x
 phi_pa=tau*tau*(I(0.5)*x2-I(1.0/6.0)*x3+I(1.0/24.0)*x4)
 phi_Sa=tau*tau*tau*(I(1.0/6.0)*x3-I(1.0/24.0)*x4+I(1.0/120.0)*x5)
 z=Interval.point(0.0); o=Interval.point(1.0)
 return [[o,z,z,phi_va],[hi,o,z,phi_pa],[I(0.5*h*h),hi,o,phi_Sa],[z,z,z,alpha]]

def _predict_lower(L,F,Qlower):
 M=matrix_mul(matrix_mul(F,_point_matrix(L)),matrix_transpose(F))
 for i in range(4): M[i][i]=M[i][i]+I(Qlower[i])
 return _center_radius_lower(M)

def _scalar_measurement_lower(L,idx,rvar):
 P=_point_matrix(L); den=P[idx][idx]+I(rvar)
 if den.lo<=0.0: raise RuntimeError('measurement lower comparison lost positive innovation')
 M=[[P[i][j]-(P[i][idx]*P[idx][j])/den for j in range(4)] for i in range(4)]
 return _center_radius_lower(M)

def _spd_delta(Omega,upper,delta):
 A=_point_matrix(Omega)
 for i in range(4): A[i][i]=A[i][i]-I(delta*upper[i])
 ok,_=symmetric_positive_definite_ldlt(matrix_symmetric_hull(A)); return ok

def _certified_delta(Omega,upper):
 lo=0.0; hi=1.0
 # log search first because the legacy number is extraordinarily small.
 trial=1.0e-36
 while trial<1.0 and _spd_delta(Omega,upper,trial):
  lo=trial; trial*=10.0
 hi=min(1.0,trial)
 if lo==0.0:
  return 0.0
 for _ in range(64):
  mid=math.sqrt(lo*hi)
  if _spd_delta(Omega,upper,mid): lo=mid
  else: hi=mid
 return down(lo)

def _member_float(text,name):
 m=re.search(rf'float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;',text)
 if not m: raise RuntimeError(f'cannot source-bind {name}')
 return float(m.group(1))

def _mode(mode,p,diag,v):
 w=diag['modes'][mode]['p3_worst_cell']; h=float(P3BASE.source_schedule()['dt_s'])
 xlo,xhi=map(float,w['x_h_over_tau']); tau=Interval.outward_bounds(h/xhi,h/xlo)
 sigma_lo=float(w['sigma_aw'][0]); rho=float(w['rho_translation_lower'])
 # Certified P3 per-sample process lower in physical coordinates.
 scales=[sigma_lo*h,sigma_lo*h*h,sigma_lo*h*h*h,sigma_lo]
 Qlower=[down(rho*s*s) for s in scales]
 if min(Qlower)<=0.0: raise RuntimeError('process lower lost positivity')
 p3=DIAG.P3.build(p); ua=list(map(float,p3['modes'][mode]['matrix_comparison']['Sigma_diagonal_upper']))
 upper=[ua[6],ua[9],ua[12],ua[15]]
 wrapper=WRAPPER.read_text(encoding='utf-8'); rho_h=min(_member_float(wrapper,'R_S_x_factor_'),_member_float(wrapper,'R_S_y_factor_'),1.0)
 rs_std=down(rho_h*float(w['R_S'][0])); rs_var=down(rs_std*rs_std)
 racc=float(v['configured_measurement_bounds']['acc_measurement_std_mps2']); racc_var=down(racc*racc)
 F=_transition_interval(tau,h); L=[[0.0]*4 for _ in range(4)]; maxrad=0.0
 n=int(math.ceil(HORIZON_S/h))
 for k in range(n):
  if k==0:
   # First prediction starts from zero, so the certified diagonal process lower
   # is already a point lower without an interval-F multiplication.
   L=[[0.0]*4 for _ in range(4)]
   for i in range(4): L[i][i]=Qlower[i]
  else:
   L,r=_predict_lower(L,F,Qlower); maxrad=max(maxrad,r)
  L,r=_scalar_measurement_lower(L,2,rs_var); maxrad=max(maxrad,r)
  L,r=_scalar_measurement_lower(L,3,racc_var); maxrad=max(maxrad,r)
 delta=_certified_delta(L,upper)
 old=float(w['delta_translation_lower'])
 return {
  'source_cell':w,'tau_interval_s':tau.as_list(),'steps':n,'process_Q_point_lower_diagonal':Qlower,
  'artificial_measurement_policy':'S_zero_and_accelerometer_accepted_every_sample_to_overinform_lower_comparison',
  'S_zero_std_used_for_max_information':rs_std,'acc_std':racc,'max_interval_radius_norm_removed_each_operation':maxrad,
  'endpoint_Omega_translation_lower':L,'translation_covariance_upper_diagonal':upper,
  'complete_word_translation_margin_lower':delta,'old_single_seed_translation_margin_lower':old,
  'margin_widening_factor_lower':down(delta/old) if delta>0 else 0.0,
  'interval_ldlt_endpoint_recertified':bool(delta>0 and _spd_delta(L,upper,delta)),
 }

def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); diag=DIAG.build(p); v=VECTOR.build(); f=[f'diagnostic: {x}' for x in DIAG.validate(diag)]+[f'vector: {x}' for x in VECTOR.validate(v)]
 modes={}
 if not f:
  for mode in ('H','A'):
   try: modes[mode]=_mode(mode,p,diag,v)
   except Exception as e: f.append(f'{mode}: {e}')
 return {'qualification':'OU3_P4_VALIDATED_WORST_CELL_COMPLETE_WORD_TRANSLATION_DISSIPATION','source_only':True,'trajectory_replay_used':False,'outward_rounded':True,'horizon_s':HORIZON_S,'modes':modes,'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS':'PASS' if not f and len(modes)==2 else 'NOT_ESTABLISHED','P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED','remaining_obligation':'extend complete-word lower propagation to every reachable source cell/edge and attitude-bias blocks, then validate exact nonlinear return map','failures':f}
def validate(d):
 f=list(d.get('failures',[]))
 if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False or d.get('outward_rounded') is not True: f.append('qualification flags invalid')
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not float(m.get('complete_word_translation_margin_lower',0.0))>0.0: f.append(f'{mode}: no complete-word translation margin')
  if m.get('interval_ldlt_endpoint_recertified') is not True: f.append(f'{mode}: endpoint margin not LDLT recertified')
  if not float(m.get('margin_widening_factor_lower',0.0))>1.0: f.append(f'{mode}: full-word translation did not widen the seed margin')
 if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED': f.append('partial translation result prematurely promoted P4')
 return f

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'translation_status':d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],'modes':{m:{'delta':d.get('modes',{}).get(m,{}).get('complete_word_translation_margin_lower'),'factor':d.get('modes',{}).get(m,{}).get('margin_widening_factor_lower')} for m in ('H','A')},'failures':f},indent=2,sort_keys=True)); return 0 if not f else 2
if __name__=='__main__': raise SystemExit(main())
