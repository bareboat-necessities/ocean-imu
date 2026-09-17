#!/usr/bin/env python3
"""Non-promoting diagnostics for the corrected ALT observability formulation."""
from __future__ import annotations
import argparse, json
import numpy as np
from tools.stability.ou3_alt_contraction import magnetic_service_formulation as F

def structural_diagnostic():
 T=3.0; U=np.array([[1.0,T],[0.0,1.0]])
 yaw_only=U[1:,1:]
 return {"qualification":"OU3_ALT_OBSERVABILITY_STRUCTURAL_DIAGNOSTIC_V1","old_pair":U.tolist(),"old_pair_spectral_radius":float(max(abs(np.linalg.eigvals(U)))),"yaw_only_quotient_remaining_bg_ratio":float(yaw_only[0,0]**2),"yaw_only_quotient_strictly_contracting":False,"selected_formulation":F.selection_status(),"shipping_superword_rho_measured":False,"note":"Structural checks only. No source-qualified shipping superword is synthesized here."}

def main():
 p=argparse.ArgumentParser(); p.add_argument("--json",action="store_true"); a=p.parse_args(); r=structural_diagnostic(); print(json.dumps(r,indent=2) if a.json else r)
if __name__=="__main__": main()
