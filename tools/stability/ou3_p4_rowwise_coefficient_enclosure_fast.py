#!/usr/bin/env python3
"""Compatibility entry point for the P4 rowwise coefficient certificate.

The base certificate now consumes the canonical endpoint-referenced factored
Riccati facade directly.  That facade already performs the dependency-reduced
small-x evaluation, so installing another monkey patch is both unnecessary and
unsafe: it can detach the facade's saved endpoint-referenced BASE build.

Keep this filename for existing proof consumers/CI, but delegate unchanged.
The rowwise correction number is diagnostic only; it is not a reset-domain
certificate and cannot be promoted into the P4 storage path.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_p4_rowwise_coefficient_enclosure as BASE

def build():return BASE.build()
def validate(d):return BASE.validate(d)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d)
    d['smallx_scaled_dependency_reduction_used']=True
    d['smallx_dependency_reduction_source']='canonical_endpoint_referenced_factored_Riccati_facade'
    d['additional_raw_BASE_monkey_patch_used']=False
    d['base_certificate_unchanged']=True
    d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({
        'Joseph_gain_family':d['Joseph_gain_family_outwardly_bounded'],
        'Kalman_reset_family':d['Kalman_and_reset_coefficient_family_outwardly_bounded'],
        'H18_rowbox_correction_diagnostic':d['modes']['H18']['independent_rowbox_attitude_correction_norm_upper_diagnostic'],
        'A21_rowbox_correction_diagnostic':d['modes']['A21']['independent_rowbox_attitude_correction_norm_upper_diagnostic'],
        'reset_domain_same_graph_required':d['reset_coefficient_family_requires_same_graph_correction_domain'],
        'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
