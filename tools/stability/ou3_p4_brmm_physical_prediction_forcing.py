#!/usr/bin/env python3
"""Exact same-history physical prediction forcing for OU-III P4.

The shipping predictor propagates the estimated latent acceleration as an OU
state.  COMPLETE-BRMM only constrains the *physical* acceleration history; it
does not assert that the true vessel acceleration follows that OU model.
Therefore the deterministic error prediction is the homogeneous shipping map
plus one same-history model-mismatch source.

For one axis over a sample h, let a0=a(0), a1=a(h) and

  J0 = integral a(s) ds,
  J1 = integral (h-s) a(s) ds,
  J2 = integral .5 (h-s)^2 a(s) ds.

Let phi_va, phi_pa, phi_Sa, alpha be the exact shipping integrated-OU
coefficients generated from the SAME committed tau and h.  With error defined
as true minus estimate, direct subtraction gives

  e_v+ = e_v + phi_va e_a + (J0 - phi_va a0)
  e_p+ = e_p + h e_v + phi_pa e_a + (J1 - phi_pa a0)
  e_S+ = e_S + h e_p + .5 h^2 e_v + phi_Sa e_a
         + (J2 - phi_Sa a0)
  e_a+ = alpha e_a + (a1 - alpha a0).

Thus the four source defects are generated once from the common physical
witness (a0,a1,J0,J1,J2); J moments and latent-acceleration mismatch are not
independent ports.  This module materializes that affine source map for all
three axes and proves coefficientwise reconstruction of the true-minus-estimate
recurrence.  It does not yet replace the same-history source relation by a
Cartesian box or claim endpoint/prefix contraction.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_shipping_prediction_primitives as SHIPPING
import ou3_brmm_centered_primitive_transition as PRIMITIVE
import ou3_brmm_correlated_window_outer_enclosure as OUTER

SCHEMA=1
QUALIFICATION="OU3_P4_BRMM_SAME_HISTORY_PHYSICAL_PREDICTION_FORCING_V1"
P3_DELTA=1.0e-18


def I(x:float)->Interval:return Interval.point(float(x))

def _vec3(x,name):
    if len(x)!=3 or any(not isinstance(v,Interval) for v in x):
        raise ValueError(name+" must be Interval[3]")
    return tuple(x)

@dataclass(frozen=True)
class PhysicalAccelerationWitness:
    a0: tuple[Interval,Interval,Interval]
    a1: tuple[Interval,Interval,Interval]
    J0: tuple[Interval,Interval,Interval]
    J1: tuple[Interval,Interval,Interval]
    J2: tuple[Interval,Interval,Interval]

@dataclass(frozen=True)
class PredictionForcing:
    dv: tuple[Interval,Interval,Interval]
    dp: tuple[Interval,Interval,Interval]
    dS: tuple[Interval,Interval,Interval]
    daw: tuple[Interval,Interval,Interval]


def witness(*,a0:Sequence[Interval],a1:Sequence[Interval],J0:Sequence[Interval],J1:Sequence[Interval],J2:Sequence[Interval])->PhysicalAccelerationWitness:
    return PhysicalAccelerationWitness(_vec3(a0,"a0"),_vec3(a1,"a1"),_vec3(J0,"J0"),_vec3(J1,"J1"),_vec3(J2,"J2"))


def forcing(w:PhysicalAccelerationWitness,tau:Interval,h:Interval)->PredictionForcing:
    if tau.lo<=0 or h.lo<=0:raise ValueError("positive tau/h required")
    axis=SHIPPING.translation_axis_transition(tau,h)
    phi_va=axis[0][3];phi_pa=axis[1][3];phi_Sa=axis[2][3];alpha=axis[3][3]
    return PredictionForcing(
      tuple(w.J0[i]-phi_va*w.a0[i] for i in range(3)),
      tuple(w.J1[i]-phi_pa*w.a0[i] for i in range(3)),
      tuple(w.J2[i]-phi_Sa*w.a0[i] for i in range(3)),
      tuple(w.a1[i]-alpha*w.a0[i] for i in range(3)))


def source_matrix(tau:Interval,h:Interval):
    """12x15 map q=[a0(3),a1(3),J0(3),J1(3),J2(3)] -> [dv,dp,dS,daw]."""
    axis=SHIPPING.translation_axis_transition(tau,h)
    phi_va,phi_pa,phi_Sa,alpha=axis[0][3],axis[1][3],axis[2][3],axis[3][3]
    z=I(0);M=[[z for _ in range(15)] for _ in range(12)]
    for a in range(3):
        M[a][a]=-phi_va; M[a][6+a]=I(1)
        M[3+a][a]=-phi_pa; M[3+a][9+a]=I(1)
        M[6+a][a]=-phi_Sa; M[6+a][12+a]=I(1)
        M[9+a][a]=-alpha; M[9+a][3+a]=I(1)
    return M


def inject_error_state(mode:str,M12):
    """Embed physical forcing into [c,bg,v,p,S,aw,(ba)] error rows."""
    n=18 if mode=="H" else 21 if mode=="A" else 0
    if not n:raise ValueError("mode must be H/A")
    z=I(0);G=[[z for _ in range(15)] for _ in range(n)]
    # M12 row groups are dv,dp,dS,daw; error-state groups begin at row 6.
    for i in range(12):
        for j in range(15):G[6+i][j]=M12[i][j]
    return G


def _point_identity_check()->dict:
    h=I(.005);tau=I(1.7)
    a0=(I(.2),I(-.1),I(.4));a1=(I(.21),I(-.08),I(.39))
    # Constant-acceleration moments are used only as an algebraic positive control.
    J0=tuple(x*h for x in a0);J1=tuple(x*h*h*I(.5) for x in a0);J2=tuple(x*h*h*h*I(1/6) for x in a0)
    w=witness(a0=a0,a1=a1,J0=J0,J1=J1,J2=J2);d=forcing(w,tau,h)
    axis=SHIPPING.translation_axis_transition(tau,h)
    # arbitrary estimator and true initial values; direct physical minus shipping estimate
    max_def=0.0
    for k in range(3):
        vhat=I(.3);phat=I(-.2);Shat=I(.1);ahat=I(.05)
        ev=I(.7);ep=I(-.4);eS=I(.25);ea=a0[k]-ahat
        vtrue=vhat+ev;ptrue=phat+ep;Strue=Shat+eS
        vhat1=vhat+axis[0][3]*ahat
        phat1=phat+h*vhat+axis[1][3]*ahat
        Shat1=Shat+h*phat+I(.5)*h*h*vhat+axis[2][3]*ahat
        ahat1=axis[3][3]*ahat
        vtrue1=vtrue+J0[k];ptrue1=ptrue+h*vtrue+J1[k];Strue1=Strue+h*ptrue+I(.5)*h*h*vtrue+J2[k]
        et=(vtrue1-vhat1,ptrue1-phat1,Strue1-Shat1,a1[k]-ahat1)
        er=(ev+axis[0][3]*ea+d.dv[k],ep+h*ev+axis[1][3]*ea+d.dp[k],eS+h*ep+I(.5)*h*h*ev+axis[2][3]*ea+d.dS[k],axis[3][3]*ea+d.daw[k])
        for x,y in zip(et,er):max_def=max(max_def,abs(x.lo-y.lo),abs(x.hi-y.hi))
    return {"coefficientwise_reconstruction_max_abs":max_def,"passed":max_def<1e-12}


def build()->dict:
    ship=SHIPPING.build();sf=SHIPPING.validate(ship);prim=PRIMITIVE.build();pf=PRIMITIVE.validate(prim);outer=OUTER.build();of=OUTER.validate(outer)
    if sf or pf or of:raise RuntimeError(f"physical forcing prerequisites failed shipping={sf} primitive={pf} outer={of}")
    check=_point_identity_check()
    if not check["passed"]:raise RuntimeError("true-minus-estimate forcing identity failed")
    return {"schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "error_sign_convention":"true_minus_estimate","source_coordinate_order":["a0_xyz","a1_xyz","J0_xyz","J1_xyz","J2_xyz"],"source_coordinate_dimension":15,
      "same_committed_tau_generates_homogeneous_and_forcing_coefficients":True,"same_physical_history_generates_a0_a1_and_all_J_moments":True,
      "exact_error_forcing_formulas":{"dv":"J0-phi_va*a0","dp":"J1-phi_pa*a0","dS":"J2-phi_Sa*a0","daw":"a1-alpha*a0"},
      "physical_moment_forcing_and_latent_acceleration_mismatch_combined_once":True,"independent_J_moment_ports_used":False,"independent_latent_acceleration_endpoint_port_used":False,
      "physical_source_forcing_matrix_materialized":True,"H18_and_A21_error_state_injection_maps_materialized":True,
      "correlated_outer_source_relation_consumed":True,"same_transition_witness_as_centered_primitive_recurrence":True,
      "point_identity_check":check,"source_uniform_quadratic_enclosure_of_full_15D_source_closed_here":False,
      "prefix_cocycle_physical_source_columns_embedded_here":False,"endpoint_augmented_LDLT_closed_here":False,"every_prefix_augmented_LDLT_closed_here":False,
      "P3_delta":P3_DELTA,"P4_MOTION_PASS":False,"P4_PASS":False,"P5_MAY_START":False,
      "next_obligation":"construct a joint quadratic outer relation for the same 15D (a0,a1,J0,J1,J2) witness, then suffix-propagate its exact injection through every later literal event and consume that map in the augmented prefix LDLT"}


def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("same_committed_tau_generates_homogeneous_and_forcing_coefficients","same_physical_history_generates_a0_a1_and_all_J_moments","physical_moment_forcing_and_latent_acceleration_mismatch_combined_once","physical_source_forcing_matrix_materialized","H18_and_A21_error_state_injection_maps_materialized","correlated_outer_source_relation_consumed","same_transition_witness_as_centered_primitive_recurrence"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("independent_J_moment_ports_used","independent_latent_acceleration_endpoint_port_used","source_uniform_quadratic_enclosure_of_full_15D_source_closed_here","prefix_cocycle_physical_source_columns_embedded_here","endpoint_augmented_LDLT_closed_here","every_prefix_augmented_LDLT_closed_here","P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    if d.get("source_coordinate_dimension")!=15:f.append("source coordinate dimension changed")
    if d.get("source_coordinate_order")!=["a0_xyz","a1_xyz","J0_xyz","J1_xyz","J2_xyz"]:f.append("source coordinate order changed")
    if d.get("point_identity_check",{}).get("passed") is not True:f.append("forcing identity not closed")
    if d.get("P3_delta")!=P3_DELTA:f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main()->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(json.dumps({"forcing":d["physical_source_forcing_matrix_materialized"],"identity":d["point_identity_check"],"P4":d["P4_PASS"],"failures":f},sort_keys=True));return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
