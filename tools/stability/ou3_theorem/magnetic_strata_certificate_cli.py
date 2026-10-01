"""Emit current magnetic 1-s exhaustive certificate status.

Until a theorem-domain seed cover is constructed, this intentionally runs
with no seeds and emits the exact unresolved coverage reason.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .magnetic_strata_certificate import exhaustive_certificate
from .magnetic_seed_cover import theorem_seed_cover
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()
 seeds,audit=theorem_seed_cover()
 z=exhaustive_certificate(seeds)
 z["seed_audit"]=audit
 Path(a.output).write_text(json.dumps(z,indent=2,sort_keys=True)+"\n")
 print(json.dumps(z,sort_keys=True))
if __name__=="__main__":main()
