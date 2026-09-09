#!/usr/bin/env python3
"""P4 rowwise coefficient certificate using dependency-reduced small-x tube.

This is a proof-arithmetic adapter only.  It installs the algebraically
identical pre-scaled small-x OU covariance evaluator while BASE builds the
source-uniform Riccati tube, then delegates the unchanged coefficient/reset
certificate.  It exists to avoid an interval dependency/subdivision explosion;
no deployed code, source domain, or certificate inequality is changed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_riccati_tube_smallx_scaled as FAST
import ou3_p4_rowwise_coefficient_enclosure as BASE


def build():
    with FAST.installed(max_depth=14):
        return BASE.build()

def validate(d):return BASE.validate(d)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['smallx_scaled_dependency_reduction_used']=True;d['base_certificate_unchanged']=True;d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'kalman_reset_family':d['Kalman_and_reset_coefficient_family_outwardly_bounded'],'H18_correction':d['modes']['H18']['attitude_correction_norm_upper'],'A21_correction':d['modes']['A21']['attitude_correction_norm_upper'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
