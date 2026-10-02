# ruff: noqa: F401, F811
"""Source-uniform causal history-cover generator.

This is the sole admissible cover driver for the constructive shaped certificate.
It subdivides shared theorem-input coordinates only and delegates ALL generated
state to one literal causal propagator.
"""
from __future__ import annotations
from dataclasses import dataclass
import heapq,math
from .reachable_history_enclosure import HistoryCell,FORBIDDEN_SPLIT

@dataclass(frozen=True)
class CoverResult:
 verified:bool;leaves:tuple;unresolved:tuple;split_count:int;reason:str

def normalized_width(cell,name,scales):
 x=dict(cell.coordinates)[name];s=max(float(scales.get(name,1.0)),1e-30);return x.width/s

def choose_split(cell,diagnostic,scales):
 """Split only shared inputs; prefer causal sensitivity supplied by propagator."""
 names=[k for k,_ in cell.coordinates if k not in FORBIDDEN_SPLIT]
 if not names:raise ArithmeticError("no admissible history split coordinate")
 scores=diagnostic.get("input_sensitivity",{}) if isinstance(diagnostic,dict) else {}
 return max(names,key=lambda k:(float(scores.get(k,0.0)),normalized_width(cell,k,scales)))

def generate(root_cells,propagate_leaf,certify_leaf,*,scales=None,max_cells=100000,max_depth=40):
 scales=scales or {};heap=[];counter=0;verified=[];unresolved=[];splits=0
 def push(cell,depth):
  nonlocal counter;counter+=1
  # widest normalized root cells first
  w=max((normalized_width(cell,k,scales) for k,_ in cell.coordinates),default=0.)
  heapq.heappush(heap,(-w,counter,depth,cell))
 for c in root_cells:push(c,0)
 while heap:
  _,_,depth,cell=heapq.heappop(heap)
  try:leaf=propagate_leaf(cell)
  except ArithmeticError as err:
   diag={"reason":str(err),"input_sensitivity":getattr(err,"input_sensitivity",{})}
   leaf=None
  if leaf is not None:
   cert=certify_leaf(leaf)
   if cert.get("verified",False):
    verified.append((cell,leaf,cert));continue
   diag=cert
   if cert.get("outside_admissible_domain",False):continue
  if depth>=max_depth or len(verified)+len(unresolved)+len(heap)+2>max_cells:
   unresolved.append((cell,diag));continue
  name=choose_split(cell,diag,scales);a,b=cell.split(name);splits+=1;push(a,depth+1);push(b,depth+1)
 return CoverResult(not unresolved,tuple(verified),tuple(unresolved),splits,
                    "complete" if not unresolved else "unresolved history cells")

def max_ratio(result:CoverResult):
 if not result.verified:return {"verified":False,"reason":result.reason,"unresolved":len(result.unresolved)}
 vals=[x[2]["ratio"].ratio_upper for x in result.leaves]
 if not vals:return {"verified":False,"reason":"empty admissible cover"}
 i=max(range(len(vals)),key=vals.__getitem__)
 return {"verified":True,"max_ratio_upper":vals[i],"leaf_count":len(vals),
         "limiting_token":result.leaves[i][0].prefix_token,"split_count":result.split_count,
         "generated_outputs_independently_boxed":False}
