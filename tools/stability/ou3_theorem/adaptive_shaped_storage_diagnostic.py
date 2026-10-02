"""Non-promoting carried shaped-storage diagnostic on generated tuner history."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .adaptive_shaped_storage import scaling,tuner_coboundary

def analyze(path):
    tr=json.loads(Path(path).read_text())
    a=[e for e in tr["events"] if e.get("kind")=="adaptive_state"]
    if len(a)<2: raise ValueError("adaptive trace missing")
    max_step=0.; product=np.eye(4)
    for x,y in zip(a,a[1:]):
        C=tuner_coboundary(x["tau"],x["sigma_aw"],y["tau"],y["sigma_aw"])
        max_step=max(max_step,float(np.linalg.norm(C-np.eye(4),2)))
        product=C@product
    direct=scaling(a[-1]["tau"],a[-1]["sigma_aw"])@np.linalg.inv(scaling(a[0]["tau"],a[0]["sigma_aw"]))
    return {
      "profile":tr["profile"],"horizon_s":tr.get("horizon_s"),
      "adaptive_samples":len(a),
      "tau_range":[min(x["tau"] for x in a),max(x["tau"] for x in a)],
      "sigma_range":[min(x["sigma_aw"] for x in a),max(x["sigma_aw"] for x in a)],
      "RS_range":[min(x["R_S"] for x in a),max(x["R_S"] for x in a)],
      "TS_range":[min(x["T_S"] for x in a),max(x["T_S"] for x in a)],
      "maximum_single_commit_coordinate_change_2norm":max_step,
      "tuner_coboundary_telescoping_error_inf":float(np.linalg.norm(product-direct,np.inf)),
      "structures_preserved":["one generated tuner history","joint tau/sigma/R_S/T_S tuples","persistent Mahony-driven frequency/variance history"],
      "relaxations_introduced":["finite carried replay only","no source-uniform claim"],
      "theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--directory",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--horizons",default="16,30,60,100");a=p.parse_args()
 out={"kind":"NON_PROMOTING_ADAPTIVE_SHAPED_STORAGE","words":[analyze(a.directory/f"word-{prof}-{h}s.json") for h in map(int,a.horizons.split(",")) for prof in ("0","wave")]}
 a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()
