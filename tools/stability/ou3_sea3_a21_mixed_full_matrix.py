#!/usr/bin/env python3
"""Complete-SEA3 A21 mixed detectable/stable full-matrix bridge.

No eta9 packet condition or full A21 information inverse is introduced.  Early
bias release is closed from the source-generated Live covariance plus the
certified time-varying translation-noise memory.  Late release is closed by
appending the held b_a seed to an already-certified H18 margin and using the
first active Gauss-Markov process injection.  Promotion remains outside this
backend.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path
from ou3_interval import Interval, matrix_mul, matrix_transpose, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_information_composition as HINFO
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_riccati_tube as TUBE
import ou3_sea3_time_varying_translation_memory as MEM
import ou3_sea3_time_varying_translation_memory_certified as MEMORY

DEFAULT_DOMAIN=COMPLETE.DEFAULT_DOMAIN
SCHEMA=1
QUALIFICATION="OU3_COMPLETE_SEA3_A21_MIXED_FULL_MATRIX_BACKEND"
DELTA=1e-18
TWORD=3.0
N=21
OV,OP,OS,OA,OB=6,9,12,15,18
BASE_CELLS=32
MAX_DEPTH=8

def I(x): return Interval.outward_bounds(float(x),float(x))
def dn(x): return TUBE.down(float(x))
def up(x): return TUBE.up(float(x))
def zero(n):
    z=I(0); return [[z for _ in range(n)] for _ in range(n)]
def diag(v):
    M=zero(len(v))
    for i,x in enumerate(v): M[i][i]=I(x)
    return M

def sub(A,B):
    return matrix_symmetric_hull([[A[i][j]-B[i][j] for j in range(len(A))] for i in range(len(A))])
def interval4(b):
    if len(b)!=4 or any(len(r)!=4 for r in b): raise RuntimeError("memory matrix dimension")
    return matrix_symmetric_hull([[Interval(float(x[0]),float(x[1])) for x in r] for r in b])

def first_translation_upper(dynamic,live,word):
    inv=dynamic["dynamic_invariant"]; shi=float(inv["sigma_aw_filter_mps2"][1]); tlo=float(inv["tau_applied_s"][0])
    s=live["translation_seed"]; pv,pp,pS=map(float,(s["P_v"],s["P_p"],s["P_S"]))
    floors=int(word["imu_samples_upper"]); impulses=1+floors; s2=up(shi*shi); qc=up(2*s2/tlo); T=TWORD
    uv=up(pv+impulses*s2*T**2+qc*T**3/3)
    upp=up(pp+T**2*pv+impulses*s2*T**4/4+qc*T**5/20)
    uS=up(pS+T**2*pp+T**4*pv/4+impulses*s2*T**6/36+qc*T**7/252)
    return {"u":[uv,upp,uS,s2],"floors_assumed_every_sample":True,"floor_events_upper":floors,"qc_upper":qc}

def endpoint_upper(path,dynamic,process,hinfo,live,adet,word):
    p=H18._same_word_covariance_upper(path,dynamic,process,hinfo)
    hd=list(map(float,p["Pbar_diagonal_variance_upper"])); uth=max(hd[:3])
    qb=float(process["source_constants"]["gyro_bias_rw_variance_density"])
    ubg=up(float(live["constructor"]["P_bg_variance"])+qb*TWORD)
    tr=first_translation_upper(dynamic,live,word); uv,upp,uS,uaw=tr["u"]
    uba=float(adet["active_bias_process"]["uniform_release_and_active_variance_upper"])
    u=[uth]*3+[ubg]*3+[uv]*3+[upp]*3+[uS]*3+[uaw]*3+[uba]*3
    if len(u)!=N or any(not(math.isfinite(x) and x>0) for x in u): raise RuntimeError("endpoint covariance upper")
    return {"marginal":u,"loewner_diag":[up(N*x) for x in u],"translation":tr,
            "rule":"P<=21*diag(u) from normalized PSD trace<=21","theta_upper":uth,"bg_upper":ubg,"ba_upper":uba}

def att_bounds(process,U):
    a=process["attitude_gyro_bias"]; c=float(a["cross_norm_upper"])
    lth=dn(float(a["theta_diagonal_lower"])-c); lbg=dn(float(a["gyro_bias_diagonal_lower"])-c)
    if not(lth>0 and lbg>0): raise RuntimeError("directional Q lower")
    hhi=float(process["configured_runtime"]["imu_dt_outward_interval_s"][1]); pc=process["source_constants"]
    qg=up(max(float(x)**2 for x in pc["gyro_noise_density_rad_sqrt_s_per_axis"])); qb=float(pc["gyro_bias_rw_variance_density"])
    qth=up(qg*hhi+qb*hhi**3/3); qbg=up(qb*hhi); qx=up(qb*hhi*hhi/2)
    # F=[R,B;0,I], ||R||=1, ||B||<=h and ||B'y+z||^2<=2h^2||y||^2+2||z||^2.
    return {"lth":lth,"lbg":lbg,"Uth":up(U[0]+2*U[3]*hhi*hhi+qth+qx),"Ubg":up(2*U[3]+qbg+qx),
            "qth_upper":qth,"qbg_upper":qbg,"qcross_upper":qx,"Bnorm_upper":hhi}

def qtrans_trace(dynamic):
    h=float(dynamic["validated_rate_and_jump_bounds"]["dt_s"]); inv=dynamic["dynamic_invariant"]
    qc=up(2*float(inv["sigma_aw_filter_mps2"][1])**2/float(inv["tau_applied_s"][0]))
    return up(qc*(h+h**3/3+h**5/20+h**7/252))

def assemble(mth,mbg,mt,mba):
    M=zero(N)
    for a in range(3):
        M[a][a]=I(mth); M[3+a][3+a]=I(mbg); M[OB+a][OB+a]=I(mba)
        ix=[OV+a,OP+a,OS+a,OA+a]
        for i in range(4):
            for j in range(4): M[ix[i]][ix[j]]=mt[i][j]
    return matrix_symmetric_hull(M)

def cell(x,*,h,L,U,att,qtr,qba,qba_hi):
    try:
        F=MEM._sample_transition(x,h); Ft=matrix_transpose(F)
        Lp=matrix_symmetric_hull(matrix_mul(matrix_mul(F,L),Ft))
        U0=diag([U[OV],U[OP],U[OS],U[OA]])
        Up=matrix_symmetric_hull(matrix_mul(matrix_mul(F,U0),Ft))
        for i in range(4): Up[i][i]=Up[i][i]+I(qtr)
        Mt=sub(Lp,[[I(DELTA)*Up[i][j] for j in range(4)] for i in range(4)])
        mth=dn(att["lth"]-DELTA*att["Uth"]); mbg=dn(att["lbg"]-DELTA*att["Ubg"])
        mba=dn(qba-DELTA*up(U[OB]+qba_hi))
        if not(mth>0 and mbg>0 and mba>0): return False,-math.inf,"directional block"
        ok,p=symmetric_positive_definite_ldlt(assemble(mth,mbg,Mt,mba))
        return (True,min(x.lo for x in p),"") if ok else (False,-math.inf,"full 21x21 LDLT")
    except Exception as e: return False,-math.inf,f"{type(e).__name__}: {e}"
def split(x):
    m=math.sqrt(x.lo*x.hi); return Interval.outward_bounds(x.lo,m),Interval.outward_bounds(m,x.hi)
def certify(x,depth,kw,st):
    ok,p,e=cell(x,**kw)
    if ok: st["leaves"]+=1; st["worst"]=min(st["worst"],p); st["depth"]=max(st["depth"],depth); return True
    if depth>=MAX_DEPTH: st["fail"].append({"x":x.as_list(),"depth":depth,"error":e}); return False
    st["splits"]+=1; a,b=split(x); return certify(a,depth+1,kw,st) and certify(b,depth+1,kw,st)

def build(domain_path:Path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve(); complete=COMPLETE.build(path); dynamic=DYNAMIC.build(path); process=PROCESS.build(); h18=H18.build(path)
    hinfo=HINFO.build(path); adet=ADET.build(path); live=LIVE.build(path); event=EVENT.build(); word=WORD.build(path); mem=MEMORY.build(path)
    bad={"complete":COMPLETE.validate(complete),"dynamic":DYNAMIC.validate(dynamic),"process":PROCESS.validate(process),"H18":H18.validate(h18),
         "Hinfo":HINFO.validate(hinfo),"Adet":ADET.validate(adet),"live":LIVE.validate(live),"event":EVENT.validate(event),"word":WORD.validate(word),"memory":MEMORY.validate(mem)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError(f"A21 prerequisites: {bad}")
    if complete["canonical_P3_source"]!="COMPLETE_SEA3_NORMAL_LIVE_WORD" or float(complete["word_horizon_s"])!=TWORD: raise RuntimeError("source/horizon")
    if adet["eta9_point_packet_shortcut_used"] is not False or not h18["H18_prior_free_completion_closed"]: raise RuntimeError("route/H18")
    if mem.get("downstream_must_consume_full_interval_candidate") is not True or mem.get("entrywise_lower_endpoint_table_is_Loewner_certificate") is not False: raise RuntimeError("memory interface")

    p0=float(adet["H_to_A_release"]["bias_diagonal_floor_variance"]); qba=float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    qbad=float(process["source_constants"]["accel_bias_process_variance_density"]); h=float(dynamic["validated_rate_and_jump_bounds"]["dt_s"]); qba_hi=up(qbad*h)
    late=dn((1-DELTA)*qba-DELTA*p0); late_ok=late>0

    ep=endpoint_upper(path,dynamic,process,hinfo,live,adet,word); U=ep["loewner_diag"]; att=att_bounds(process,U)
    L=interval4(mem["word_endpoint_translation_process_measurement_noise_interval_lower"]); qtr=qtrans_trace(dynamic)
    tlo,thi=map(float,dynamic["dynamic_invariant"]["tau_applied_s"]); xlo,xhi=h/thi,h/tlo
    st={"leaves":0,"splits":0,"depth":0,"worst":math.inf,"fail":[]}; kw={"h":h,"L":L,"U":U,"att":att,"qtr":qtr,"qba":qba,"qba_hi":qba_hi}
    early=True
    for x in TUBE.interval_cells(TUBE.geom_edges(xlo,xhi,BASE_CELLS)):
        if not certify(x,0,kw,st): early=False
    early=early and not st["fail"] and st["leaves"]>0 and math.isfinite(st["worst"]) and st["worst"]>0
    pr=event["full_matrix_margin_preservation"]; suffix=all(bool(pr[k]) for k in ("covers_prediction","covers_every_due_S_update","covers_every_Normal_Live_accelerometer_update","covers_asynchronous_magnetometer_update","covers_immediate_left_error_reset","covers_aw_covariance_floor","covers_not_due_or_rejected_identity_branches"))
    closed=early and late_ok and suffix
    return {"schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":complete["canonical_P3_source"],"useful_gate":DELTA,
            "eta9_point_packet_shortcut_used":False,"full_A21_D_inverse_used":False,"constant_tau_over_word_assumed":False,"source_history_graph_consumed":False,
            "time_varying_translation_memory_full_interval_consumed":True,"translation_entrywise_lower_diagnostic_consumed":False,
            "early":{"covers_any_release_before_H18_margin":True,"active_A_measurement_bound_used_for_whole_word":True,"endpoint_upper":ep,"attitude_bias":att,
                     "qtranslation_trace_upper":qtr,"leaves":st["leaves"],"splits":st["splits"],"max_depth":st["depth"],"failures":st["fail"],
                     "worst_full_21x21_LDLT_pivot_lower":st["worst"] if early else None,"full_21x21_interval_LDLT_used":True,"closed":early},
            "late":{"H18_ba_cross_zero_at_release":True,"bias_prior_variance":p0,"Qba_lower":qba,"ba_margin_lower":late,"closed":late_ok,"another_H_word_required":False},
            "actual_applied_SpectralMSE_R_S_retained":True,"event_algebra_preserves_margin":suffix,"A21_full_matrix_closed":closed,
            "source_family_replaced":False,"trajectory_replay_used":False,"independent_tau_sigma_RS_source_created":False,"P3_promoted":False}

def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION: f.append("schema/qualification")
    if d.get("canonical_source")!="COMPLETE_SEA3_NORMAL_LIVE_WORD": f.append("source")
    for k in ("time_varying_translation_memory_full_interval_consumed","actual_applied_SpectralMSE_R_S_retained","event_algebra_preserves_margin","A21_full_matrix_closed"):
        if d.get(k) is not True: f.append(k)
    for k in ("eta9_point_packet_shortcut_used","full_A21_D_inverse_used","constant_tau_over_word_assumed","source_history_graph_consumed","translation_entrywise_lower_diagnostic_consumed","source_family_replaced","trajectory_replay_used","independent_tau_sigma_RS_source_created","P3_promoted"):
        if d.get(k) is not False: f.append(k)
    e=d.get("early",{}); l=d.get("late",{})
    if e.get("closed") is not True or e.get("full_21x21_interval_LDLT_used") is not True or e.get("failures"): f.append("early bridge")
    p=e.get("worst_full_21x21_LDLT_pivot_lower");
    if not isinstance(p,(int,float)) or not(math.isfinite(float(p)) and float(p)>0): f.append("early pivot")
    if l.get("closed") is not True or l.get("another_H_word_required") is not False or not float(l.get("ba_margin_lower",0))>0: f.append("late bridge")
    if float(d.get("useful_gate",math.nan))!=DELTA: f.append("delta")
    return list(dict.fromkeys(f))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(a.domain); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({"early":d["early"]["closed"],"pivot":d["early"]["worst_full_21x21_LDLT_pivot_lower"],"late_margin":d["late"]["ba_margin_lower"],"A21":d["A21_full_matrix_closed"],"failures":f},indent=2)); return 0 if not f else 2
if __name__=="__main__": raise SystemExit(main())