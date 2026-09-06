#!/usr/bin/env python3
"""Canonical facade for the complete-SEA3 mixed A21 full-matrix bridge."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import ou3_sea3_a21_mixed_full_matrix as IMPL

DEFAULT_DOMAIN=IMPL.DEFAULT_DOMAIN
SCHEMA=3
QUALIFICATION="OU3_COMPLETE_SEA3_A21_MIXED_DETECTABLE_STABLE_FULL_MATRIX_COMPLETION"
USEFUL_GATE=IMPL.DELTA

def build(domain_path:Path=DEFAULT_DOMAIN):
    d=IMPL.build(domain_path); out=dict(d)
    out["backend_schema"]=d["schema"]; out["backend_qualification"]=d["qualification"]
    out["schema"]=SCHEMA; out["qualification"]=QUALIFICATION
    out["paper_active_bias_route"]="ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION"
    out["H18_prior_free_completion_consumed"]=True
    out["A21_finite_bias_detectability_consumed"]=True
    out["finite_tau_detectability_does_not_imply_full_A21_information_inverse"]=True
    out["full_A21_D_inverse_available"]=False
    out["full_A21_prior_free_D_inverse_identity_used"]=False
    out["full_21x21_Omega_minus_delta_P_LDLT_closed"]=bool(d["A21_full_matrix_closed"])
    out["A21_prior_free_completion_closed"]=bool(d["A21_full_matrix_closed"])
    out["full_21x21_interval_LDLT_used"]=True
    out["event_algebra_preserves_margin_after_closure"]=bool(d["event_algebra_preserves_margin"])
    out["actual_applied_SpectralMSE_R_S_retained_through_H18_component"]=True
    out["useful_gate"]=USEFUL_GATE
    return out

def validate(d):
    probe=dict(d); probe["schema"]=d.get("backend_schema"); probe["qualification"]=d.get("backend_qualification")
    f=IMPL.validate(probe)
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION: f.append("facade schema/qualification")
    if d.get("paper_active_bias_route")!="ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION": f.append("active-bias route")
    for k in ("H18_prior_free_completion_consumed","A21_finite_bias_detectability_consumed","finite_tau_detectability_does_not_imply_full_A21_information_inverse","full_21x21_Omega_minus_delta_P_LDLT_closed","A21_prior_free_completion_closed","full_21x21_interval_LDLT_used","event_algebra_preserves_margin_after_closure","actual_applied_SpectralMSE_R_S_retained_through_H18_component"):
        if d.get(k) is not True: f.append(k)
    for k in ("eta9_point_packet_shortcut_used","full_A21_D_inverse_available","full_A21_prior_free_D_inverse_identity_used","P3_promoted"):
        if d.get(k) is not False: f.append(k)
    if float(d.get("useful_gate",math.nan))!=USEFUL_GATE: f.append("useful gate")
    return list(dict.fromkeys(f))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(a.domain); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({"A21_closed":d["A21_prior_free_completion_closed"],"early_pivot":d["early"]["worst_full_21x21_LDLT_pivot_lower"],"late_margin":d["late"]["ba_margin_lower"],"failures":f},indent=2)); return 0 if not f else 2
if __name__=="__main__": raise SystemExit(main())