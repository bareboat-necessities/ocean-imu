#!/usr/bin/env python3
"""Exact same-history BRMM-vs-OU translational prediction mismatch.

The shipping estimator predicts [v,p,S,a_w] with its integrated-OU transition
Phi(tau,h).  COMPLETE-BRMM physical truth is not required to obey that OU
model.  On one physical interval let a0,a1 be the true world accelerations at
its endpoints and J0,J1,J2 the exact same-history acceleration moments

  J0=int a, J1=int(h-s)a, J2=int (h-s)^2 a/2.

True primitives obey
  v1=v0+J0,
  p1=p0+h v0+J1,
  S1=S0+h p0+h^2 v0/2+J2,
while the shipping OU predictor uses coefficients cv,cp,cS,phi from Phi.
For error e=true-estimate, the exact additive model-mismatch source is therefore

  u_v  = J0 - cv a0,
  u_p  = J1 - cp a0,
  u_S  = J2 - cS a0,
  u_aw = a1 - phi a0.

and
  e_L,+ = Phi e_L + u_L.

This source term remains present even with zero sensor noise.  J0/J1/J2,a0,a1
must belong to one physical transition witness; they are not independent ports.
The function below performs only this exact algebra.  Source admissibility and
its coupled IQCs are separate graph constraints.
"""
from __future__ import annotations
import argparse,json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_shipping_prediction_primitives as SHIPPING
import ou3_brmm_centered_primitive_transition as PRIMITIVE
import ou3_brmm_acceleration_moment_iqc as MOM

SCHEMA=1
QUALIFICATION='OU3_P4_EXACT_BRMM_TRANSLATION_MODEL_MISMATCH_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))

def _v3(x,name):
    if len(x)!=3 or any(not isinstance(v,Interval) for v in x):raise ValueError(name+' must be Interval[3]')
    return tuple(x)

@dataclass(frozen=True)
class TranslationSourceWitness:
    source_transition_witness_id:str
    a0_world:tuple[Interval,Interval,Interval]
    a1_world:tuple[Interval,Interval,Interval]
    J0:tuple[Interval,Interval,Interval]
    J1:tuple[Interval,Interval,Interval]
    J2:tuple[Interval,Interval,Interval]

def witness(source_transition_witness_id:str,a0_world:Sequence[Interval],a1_world:Sequence[Interval],J0:Sequence[Interval],J1:Sequence[Interval],J2:Sequence[Interval])->TranslationSourceWitness:
    if not source_transition_witness_id:raise ValueError('source transition witness id required')
    return TranslationSourceWitness(source_transition_witness_id,_v3(a0_world,'a0'),_v3(a1_world,'a1'),_v3(J0,'J0'),_v3(J1,'J1'),_v3(J2,'J2'))

def mismatch(w:TranslationSourceWitness,tau:Interval,h:Interval):
    if not isinstance(w,TranslationSourceWitness):raise TypeError('TranslationSourceWitness required')
    if tau.lo<=0 or h.lo<=0:raise ValueError('positive tau/h required')
    Phi=SHIPPING.translation_axis_transition(tau,h)
    cv,cp,cS,phi=Phi[0][3],Phi[1][3],Phi[2][3],Phi[3][3]
    uv=[];up=[];uS=[];ua=[]
    for k in range(3):
        uv.append(w.J0[k]-cv*w.a0_world[k])
        up.append(w.J1[k]-cp*w.a0_world[k])
        uS.append(w.J2[k]-cS*w.a0_world[k])
        ua.append(w.a1_world[k]-phi*w.a0_world[k])
    # state order [v_xyz,p_xyz,S_xyz,a_xyz]
    return tuple(uv+up+uS+ua)

def normalized_moment_linear_coefficients(tau:Interval,h:Interval):
    """Per-axis u=[uv,up,uS,ua] map from q=[a0,a1,x0,x1,x2]."""
    Phi=SHIPPING.translation_axis_transition(tau,h);cv,cp,cS,phi=Phi[0][3],Phi[1][3],Phi[2][3],Phi[3][3]
    h2=h*h;h3=h2*h;z=I(0);o=I(1)
    return [
      [-cv,z,h,z,z],
      [-cp,z,z,h2,z],
      [-cS,z,z,z,h3],
      [-phi,o,z,z,z],
    ]
def _smoke():
    h=I(1/200);tau=I(2.0);a0=(I(.3),I(-.2),I(.1));a1=(I(.31),I(-.19),I(.08))
    J0=tuple(h*x for x in a0);J1=tuple(I(.5)*h*h*x for x in a0);J2=tuple(I(1/6)*h*h*h*x for x in a0)
    w=witness('tr0',a0,a1,J0,J1,J2);u=mismatch(w,tau,h);A=normalized_moment_linear_coefficients(tau,h)
    # Reconstruct using normalized moment coordinates.
    x=[]
    for k in range(3):x.append((a0[k],a1[k],J0[k]/h,J1[k]/(h*h),J2[k]/(h*h*h)))
    rec=[]
    for row in range(4):
        for axis in range(3):
            s=I(0)
            for c,q in zip(A[row],x[axis]):s=s+c*q
            rec.append((row,axis,s))
    target=list(u[:3])+list(u[3:6])+list(u[6:9])+list(u[9:12])
    # rec order is row-major by model component, same as target.
    defect=max(max(abs(s.lo-t.lo),abs(s.hi-t.hi)) for (_,_,s),t in zip(rec,target))
    return {'reconstruction_defect':defect,'mismatch_norm_component_upper':max(v.abs_upper() for v in u),'witness_id':w.source_transition_witness_id}
def build():
    ship=SHIPPING.build();sf=SHIPPING.validate(ship);prim=PRIMITIVE.build();pf=PRIMITIVE.validate(prim);mom=MOM.build();mf=MOM.validate(mom)
    if sf or pf or mf:raise RuntimeError(f'mismatch prerequisites failed shipping={sf} primitive={pf} moment={mf}')
    s=_smoke()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'physical_truth_not_assumed_to_follow_filter_OU':True,'exact_model_mismatch_formula_materialized':True,
      'same_transition_J0_J1_J2_a0_a1_required':True,'same_applied_tau_generates_OU_coefficients':True,
      'normalized_moment_linear_source_map_available':True,'source_term_persists_with_zero_sensor_noise':True,
      'independent_J0_J1_J2_a0_a1_ports_allowed':False,'shipping_filter_changed':False,
      'smoke':s,'production_provider_transition_witness_attached_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'attach a0/a1 and coupled J0/J1/J2 from each SAME provider transition; add this 12D mismatch to prediction state_out/source map and constrain a0/a1 plus moments on the common augmented graph'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('physical_truth_not_assumed_to_follow_filter_OU','exact_model_mismatch_formula_materialized','same_transition_J0_J1_J2_a0_a1_required','same_applied_tau_generates_OU_coefficients','normalized_moment_linear_source_map_available','source_term_persists_with_zero_sensor_noise'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_J0_J1_J2_a0_a1_ports_allowed','shipping_filter_changed','production_provider_transition_witness_attached_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('smoke',{}).get('reconstruction_defect',1))>1e-15:f.append('normalized source-map reconstruction failed')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'mismatch':d['exact_model_mismatch_formula_materialized'],'provider':d['production_provider_transition_witness_attached_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
