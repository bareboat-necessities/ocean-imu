#!/usr/bin/env python3
"""Validate committed theorem status and shipping-source provenance."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from theorem_status import status_report

REPO=Path(__file__).resolve().parents[3]
PROVENANCE=REPO/"reports/results/ou3_stability/provenance.json"
STATUS=REPO/"reports/results/ou3_stability/theorem-status.json"
RICCATI_STATUS=REPO/"reports/results/ou3_stability/interval-riccati-status.json"

def git_blob_sha(path: Path) -> str:
    data=path.read_bytes(); h=hashlib.sha1(); h.update(f"blob {len(data)}\0".encode("ascii")); h.update(data); return h.hexdigest()

def validate() -> dict:
    provenance=json.loads(PROVENANCE.read_text(encoding="utf-8"))
    committed=json.loads(STATUS.read_text(encoding="utf-8"))
    failures=[]
    sources = provenance["authoritative_shipping_sources"] + provenance.get("operation_lemma_sources", [])
    for row in sources:
        path=REPO/row["path"]
        if not path.is_file(): failures.append(f"missing bound source: {row['path']}"); continue
        actual=git_blob_sha(path)
        if actual!=row["git_blob_sha"]: failures.append(f"source provenance changed: {row['path']} {row['git_blob_sha']} -> {actual}")
    expected=status_report()
    if committed!=expected: failures.append("committed theorem-status.json differs from theorem_status.status_report()")
    riccati=json.loads(RICCATI_STATUS.read_text(encoding="utf-8"))
    if riccati.get("qualification")!="OU3_A21_INTERVAL_RICCATI_V1":
        failures.append("interval Riccati status has wrong qualification")
    if riccati.get("certificate_complete") is True:
        from interval_riccati import load_certificate
        verified=load_certificate(RICCATI_STATUS).verify()
        if not verified["certificate_complete"]:
            failures.append("interval Riccati status claims completion without verified inclusion")
    return {"validation_pass":not failures,"failures":failures,"base_main_commit":provenance["base_main_commit"],
            "shipping_behavior_authority":"source implementation","theorem_closed":expected["theorem_closed"]}

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--output",type=Path); args=ap.parse_args()
    report=validate()
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(report,sort_keys=True)); return 0 if report["validation_pass"] else 1

if __name__=="__main__": raise SystemExit(main())
