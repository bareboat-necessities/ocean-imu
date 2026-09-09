#!/usr/bin/env python3
"""Homogeneous dense quadratic constraints for the exact finite-reset lift.

Use one projective anchor ``s`` and the lifted coordinate

 z=[s,c,d,a,rho,x,xc,w,h,hc,hd,hw,hrho].

On the physical graph s=1.  Multiplying every defining relation by s turns all
lift identities into homogeneous quadratics, so the existing augmented
S-procedure can consume them without a cubic relaxation:

  s x = s(a-d),             s xc = x cross c,
  s w = d cross c,          s h = a^T c,
  s hc = h c,               s hd = h d,
  s hw = h w,               s hrho = h rho,

and

  s[rho+.25 hrho+x-.5 xc+.25 hc-.25 hd+.125 hw] = 0.

The deployed correction Cayley vector is also parallel to d.  Compact graph
bounds on c,d,a and |a-d| are added as inequalities.  Equalities are emitted as
both +Q and -Q, compatible with the master's nonnegative multiplier interface.

Production callers MUST use ``build_for_same_cell(delta)`` so the d-ball comes
from the same-cell Joseph correction-domain certificate.  The zero-argument
``build()`` route is retained only for the historical rowwise magnitude
smoke/diagnostic and is explicitly labelled as such; it is not theorem input.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_exact_reset_joint_coordinate as JOINT

QUALIFICATION='OU3_P4_EXACT_FINITE_RESET_HOMOGENEOUS_QC_V2'
VECTORS=('c','d','a','rho','x','xc','w','hc','hd','hw','hrho')


def layout():
    out={'s':0};k=1
    for name in VECTORS:out[name]=tuple(range(k,k+3));k+=3
    out['h']=k;k+=1
    return out,k

def zero(n):return [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
def add_cross(Q,i,j,coef):
    v=Interval.point(.5*float(coef));Q[i][j]=Q[i][j]+v;Q[j][i]=Q[j][i]+v
def add_diag(Q,i,coef):Q[i][i]=Q[i][i]+Interval.point(float(coef))
def neg(Q):return [[-x for x in row] for row in Q]
def equality_pair(Q):return (Q,neg(Q))

def product_anchor_eq(n,s,out_idx,left_idx,right_idx,coef_left=1.0,coef_product=-1.0):
    Q=zero(n);add_cross(Q,s,out_idx,coef_left);add_cross(Q,left_idx,right_idx,coef_product);return Q

def cross_component_eq(n,s,out,a,b,component):
    Q=zero(n);add_cross(Q,s,out[component],1.0)
    pairs=((1,2,2,1),(2,0,0,2),(0,1,1,0))[component]
    i,j,k,l=pairs;add_cross(Q,a[i],b[j],-1.0);add_cross(Q,a[k],b[l],1.0);return Q

def dot_eq(n,s,h,a,b):
    Q=zero(n);add_cross(Q,s,h,1.0)
    for i in range(3):add_cross(Q,a[i],b[i],-1.0)
    return Q

def anchor_linear_eq(n,s,terms):
    Q=zero(n)
    for idx,coef in terms:add_cross(Q,s,idx,coef)
    return Q

def norm_ball_qc(n,s,vec,radius):
    Q=zero(n);add_diag(Q,s,float(radius)**2)
    for i in vec:add_diag(Q,i,-1.0)
    return Q

def diff_ball_qc(n,s,a,d,radius):
    Q=zero(n);add_diag(Q,s,float(radius)**2)
    for i in range(3):
        add_diag(Q,a[i],-1.0);add_diag(Q,d[i],-1.0);add_cross(Q,a[i],d[i],2.0)
    return Q

def cross_zero_eq(n,a,b,component):
    Q=zero(n);pairs=((1,2,2,1),(2,0,0,2),(0,1,1,0))[component];i,j,k,l=pairs
    add_cross(Q,a[i],b[j],1.0);add_cross(Q,a[k],b[l],-1.0);return Q

def qvalue(Q,z):
    s=0.0
    for i,row in enumerate(Q):
        for j,v in enumerate(row):s+=z[i]*v.lo*z[j]
    return s

def _build_with_delta(delta:float, *, correction_domain_source:str, same_cell:bool):
    delta=float(delta)
    if not (math.isfinite(delta) and 0.0 <= delta <= RESET.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError('correction domain must lie inside validated exact-reset utility range')
    L,n=layout();entry=ENTRY.build();bad={'entry':ENTRY.validate(entry)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset lifted-QC prerequisites failed: '+repr(bad))
    q=float(entry['coordinate_radii']['attitude_cayley_norm'])
    cb=RESET.correction_cayley_norm_bounds(delta);amax=float(cb['injected_cayley_norm_upper']);adiff=float(cb['injected_cayley_minus_delta_norm_upper'])
    eq=[];names=[]
    def emit(name,Q):
        p,m=equality_pair(Q);eq.extend((p,m));names.extend((name+'_plus',name+'_minus'))
    for i in range(3):emit('x_equals_a_minus_d_'+str(i),anchor_linear_eq(n,L['s'],[(L['x'][i],1),(L['a'][i],-1),(L['d'][i],1)]))
    for i in range(3):emit('xc_equals_x_cross_c_'+str(i),cross_component_eq(n,L['s'],L['xc'],L['x'],L['c'],i))
    for i in range(3):emit('w_equals_d_cross_c_'+str(i),cross_component_eq(n,L['s'],L['w'],L['d'],L['c'],i))
    emit('h_equals_a_dot_c',dot_eq(n,L['s'],L['h'],L['a'],L['c']))
    for outv,base in (('hc','c'),('hd','d'),('hw','w'),('hrho','rho')):
        for i in range(3):emit(outv+'_equals_h_'+base+'_'+str(i),product_anchor_eq(n,L['s'],L[outv][i],L['h'],L[base][i]))
    for i in range(3):
        terms=[(L['rho'][i],1.0),(L['hrho'][i],.25),(L['x'][i],1.0),(L['xc'][i],-.5),(L['hc'][i],.25),(L['hd'][i],-.25),(L['hw'][i],.125)]
        emit('exact_reset_incidence_'+str(i),anchor_linear_eq(n,L['s'],terms))
    for i in range(3):emit('a_parallel_d_'+str(i),cross_zero_eq(n,L['a'],L['d'],i))
    inequalities={
      'c_ball':norm_ball_qc(n,L['s'],L['c'],q),
      'd_ball':norm_ball_qc(n,L['s'],L['d'],delta),
      'a_ball':norm_ball_qc(n,L['s'],L['a'],amax),
      'a_minus_d_ball':diff_ball_qc(n,L['s'],L['a'],L['d'],adiff),
    }
    c=[.31,-.17,.22];d0=[.18,-.09,.12];r=math.sqrt(sum(v*v for v in d0));ratio=(2*math.tan(r/2)/r) if r else 1.;a=[ratio*v for v in d0];cp,rho=JOINT.exact_reset(c,d0,a);lift=JOINT.lifted(c,d0,a,rho)
    z=[0.0]*n;z[L['s']]=1.0
    for name,val in (('c',c),('d',d0),('a',a),('rho',rho),('x',lift['x_a_minus_d']),('xc',lift['x_cross_c']),('w',lift['w_d_cross_c']),('hc',lift['h_c']),('hd',lift['h_d']),('hw',lift['h_w']),('hrho',lift['h_rho'])):
        for i in range(3):z[L[name][i]]=val[i]
    z[L['h']]=lift['h_a_dot_c']
    eq_res=max((abs(qvalue(Q,z)) for Q in eq),default=0.0)
    # The fixed algebra smoke correction may lie outside a caller-supplied
    # production delta.  Only evaluate compact inequality membership when it is
    # actually a member of that domain; equality regression remains universal.
    smoke_in_domain=r <= delta
    ineq_min=min((qvalue(Q,z) for Q in inequalities.values()),default=0.0) if smoke_in_domain else None
    finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for Q in eq+list(inequalities.values()) for row in Q for x in row)
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','dimension':n,'layout':{k:(list(v) if isinstance(v,tuple) else v) for k,v in L.items()},'same_cell_d_equals_Etheta_K_y_required_downstream':True,'projective_anchor_homogenizes_all_reset_relations':True,'exact_reset_equalities_encoded_as_plus_minus_QCs':True,'deployed_a_parallel_d_encoded_exactly':True,'compact_reset_graph_inequalities':['c_ball','d_ball','a_ball','a_minus_d_ball'],'correction_norm_upper':delta,'correction_domain_source':correction_domain_source,'same_cell_correction_domain_supplied':same_cell,'rowwise_K_correction_domain_used':not same_cell,'correction_cayley_norm_upper':amax,'correction_cayley_minus_d_norm_upper':adiff,'equality_QC_names':names,'equality_QC_count':len(eq),'inequality_QC_count':len(inequalities),'dense_QCs_finite':finite,'physical_smoke_max_equality_quadratic_residual':eq_res,'physical_smoke_inside_correction_domain':smoke_in_domain,'physical_smoke_min_inequality_value':ineq_min,'independent_reset_norm_port_used':False,'packet_count_multiplier_used':False,'source_uniform_endpoint_augmented_LDLT_closed_here':False,'source_uniform_every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False}

def build_for_same_cell(correction_norm_upper:float):
    return _build_with_delta(correction_norm_upper,correction_domain_source='SAME_CELL_JOSEPH_CORRECTION_DOMAIN',same_cell=True)

def build():
    # Historical diagnostic compatibility only.  Production theorem code must
    # call build_for_same_cell().
    coeff=COEFF.build();vf=COEFF.validate(coeff)
    if vf:raise RuntimeError('legacy rowwise reset-lift diagnostic prerequisite failed: '+repr(vf))
    delta=max(float(coeff['modes'][m]['independent_rowbox_attitude_correction_norm_upper_diagnostic']) for m in ('H18','A21'))
    return _build_with_delta(delta,correction_domain_source='LEGACY_ROWWISE_MAGNITUDE_DIAGNOSTIC_ONLY',same_cell=False)

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_cell_d_equals_Etheta_K_y_required_downstream','projective_anchor_homogenizes_all_reset_relations','exact_reset_equalities_encoded_as_plus_minus_QCs','deployed_a_parallel_d_encoded_exactly','dense_QCs_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_reset_norm_port_used','packet_count_multiplier_used','source_uniform_endpoint_augmented_LDLT_closed_here','source_uniform_every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if int(d.get('equality_QC_count',0))<50:f.append('too few exact reset equality QCs')
    if int(d.get('inequality_QC_count',0))!=4:f.append('reset compact inequality count changed')
    if not float(d.get('physical_smoke_max_equality_quadratic_residual',math.inf))<1e-12:f.append('physical reset point violates equality lift')
    if d.get('physical_smoke_inside_correction_domain') is True:
        if not float(d.get('physical_smoke_min_inequality_value',-math.inf))>=-1e-12:f.append('physical reset point violates compact graph QC')
    if d.get('same_cell_correction_domain_supplied') is True and d.get('rowwise_K_correction_domain_used') is not False:f.append('same-cell lift still uses rowwise correction domain')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'dimension':d['dimension'],'equalities':d['equality_QC_count'],'delta':d['correction_norm_upper'],'source':d['correction_domain_source'],'eq_res':d['physical_smoke_max_equality_quadratic_residual'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
