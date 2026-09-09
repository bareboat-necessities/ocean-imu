#!/usr/bin/env python3
"""Quantify why the scalar P3-margin nonlinear enclosure cannot close P4.

This is deliberately non-promoting.  It reuses the exact same constants as
``ou3_p4_uniform_closure`` but allows an arbitrarily tiny probe entry radius so
we can report the actual covariance, gain, Hessian and arithmetic scales.  It
must never set P4/P5 flags true.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_brmm_riccati_tube as TUBE
import ou3_full_process_ucc as PROCESS
import ou3_p4_uniform_closure as C


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args=ap.parse_args()
    domain=json.loads(C.DEFAULT_DOMAIN.read_text())
    closure=json.loads(C.DEFAULT_CLOSURE.read_text())
    probe=copy.deepcopy(closure)
    probe['hard_entry_search']['candidate_scale_factors']=[10.0**(-k) for k in range(1,121)]
    tube=TUBE.build(C.DEFAULT_DOMAIN)
    process=PROCESS.build()
    dynamic=DYNAMIC.build(C.DEFAULT_DOMAIN)
    modes={m:C._mode_constants(m, probe, domain, tube, process, dynamic) for m in ('H18','A21')}
    out={
      'qualification':'OU3_P4_FAILED_SCALAR_UNIFORM_ENCLOSURE_DIAGNOSTIC_V1',
      'non_promoting':True,
      'architecture':'P3 moving-Riccati delta pays whole nonlinear completed-event remainder',
      'modes':modes,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'decision':'ABANDON_THIS_SCALAR_ARCHITECTURE' if any(m['certified_hard_entry_scale'] < 1e-12 for m in modes.values()) else 'REVIEW'
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:{x:modes[k][x] for x in (
      'Pbar_lambda_max_trace_upper','completed_event_P_lambda_min_lower','source_uniform_K_norm_upper',
      'completed_event_Hessian_norm_upper','metric_condition_sqrt_upper','nonlinear_Euclidean_radius_upper',
      'certified_hard_entry_scale','finite_precision_event_abs_forcing_upper')} for k in modes},indent=2,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
