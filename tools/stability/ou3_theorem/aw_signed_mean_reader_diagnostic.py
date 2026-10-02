"""Literal carried 16-s signed-AW reader feasibility diagnostic.

This does not promote a source-uniform theorem. It reuses the committed native
AW tracking audit, whose six stress profiles include sync-locked rectification,
and compares the largest literal signed-mean error to the analytically proved
candidate full-state reader allowance.
"""
import json
from pathlib import Path
from fractions import Fraction as F
from .ag_readout_source_diagnostic import REPO
from .nonrecurring_accel_bridge import certificate as bridge

REPORT=REPO/'reports/results/ou3_stability/aw-tracking-source-feasibility.json'

def _find_numbers(obj, key_contains):
    out=[]
    if isinstance(obj,dict):
        for k,v in obj.items():
            if key_contains in k.lower() and isinstance(v,(int,float,str)):
                try: out.append((k,F(str(v))))
                except (ValueError,ZeroDivisionError): pass
            out.extend(_find_numbers(v,key_contains))
    elif isinstance(obj,list):
        for v in obj: out.extend(_find_numbers(v,key_contains))
    return out

def certificate():
    report=json.loads(REPORT.read_text())
    # The committed report names this metric by profile; tolerate naming
    # changes by searching signed-mean numeric leaves and select the maximum.
    vals=_find_numbers(report,'signed_16s_mean_error_max_mps2')
    if not vals:
        raise ValueError('committed AW audit has no dimensional signed 16-s mean error metric')
    worst_key,worst=max(vals,key=lambda kv: abs(kv[1]))
    allowance=F(bridge()['required_full_state_reader_budget_mps2'])
    return {
      'qualification':'OU3_CARRIED_SIGNED_AW_READER_16S_V1',
      'committed_report':str(REPORT.relative_to(REPO)),
      'selected_metric':worst_key,
      'worst_carried_signed_mean_error_mps2':str(worst),
      'analytic_reader_allowance_mps2':str(allowance),
      'carried_margin_mps2':str(allowance-abs(worst)),
      'carried_feasible':abs(worst)<allowance,
      'source_uniform_verified':False,
      'theorem_closed':False}

if __name__=='__main__':
 print(json.dumps(certificate(),indent=2,sort_keys=True))
