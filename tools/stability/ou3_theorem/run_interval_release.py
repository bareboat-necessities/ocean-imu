"""Execute the source-uniform closed-loop release interval calculation."""
from __future__ import annotations
import json,traceback
from .source_uniform_root import root_cells
from .literal_history_leaf_propagator import propagate_history_cell
from .closed_loop_release_stream import generate_interval_stream

def main():
 root=root_cells(60.)[0]
 try:
  payload=propagate_history_cell(root,60.,.005)
 except Exception as e:
  traceback.print_exc()
  print(json.dumps({"verified":False,"stage":"history_leaf","first_failure_sample":0,
                    "reason":str(e),"classification":"D_ENCLOSURE_FAILURE"},indent=2))
  return 2
 out=generate_interval_stream(payload,.005)
 slim={k:v for k,v in out.items() if k not in ("P","events","mean_mid","mean_rad")}
 slim["event_count"]=len(out.get("events",()))
 if "mean_rad" in out:
  r=out["mean_rad"];slim["terminal_radii"]={
   "AW":float(max(r[15:18])),"BA":float(max(r[18:21])),
   "S":float(max(r[12:15]))}
 print(json.dumps(slim,indent=2))
 return 0 if out.get("verified") else 2

if __name__=="__main__":raise SystemExit(main())
