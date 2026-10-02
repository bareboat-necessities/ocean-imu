"""Diagnose physical HistoryCell subdivision for causal guard/release failures."""
from __future__ import annotations
import json,re
from .source_uniform_root import root_cells
from .literal_history_leaf_propagator import propagate_history_cell
from .physical_history_subdivision import split_guard_failure,assert_no_generated_coordinates
from .enclosure_failure import EnclosureFailure

def failure_time(msg):
 m=re.search(r"t=([0-9.]+)",msg);return float(m.group(1)) if m else 0.

def main(depth=6):
 cells=list(root_cells(60.));trace=[]
 for level in range(depth):
  nxt=[]
  for cell in cells:
   assert_no_generated_coordinates(cell)
   try:
    propagate_history_cell(cell,60.,.005)
    trace.append({"level":level,"status":"history_leaf_survived","token":cell.prefix_token})
   except (EnclosureFailure,ArithmeticError,ValueError,OverflowError) as e:
    t=e.time if isinstance(e,EnclosureFailure) and e.time is not None else failure_time(str(e));z=split_guard_failure(cell,t)
    trace.append({"level":level,"status":"split","failure_time":t,
                  "coordinate":z["coordinate"],"reason":str(e)})
    nxt.extend((z["left"],z["right"]))
  if not nxt:break
  cells=nxt
 print(json.dumps({"trace":trace,"remaining_cells":len(cells),
                   "physical_only":True,"depth":depth},indent=2))
 return 0
if __name__=="__main__":raise SystemExit(main())
