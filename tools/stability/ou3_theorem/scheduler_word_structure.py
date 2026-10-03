"""Exact scheduler-word structure and neighboring-cell differences."""
from __future__ import annotations
from .scheduler_branch_enumeration import representatives,due_word

def adjacent_word_differences(period=.1363605111837387,dt=.005,steps=4000):
 cells=representatives(period,dt,steps)
 words=[due_word(m,period,dt,steps)[0] for _,_,m in cells]
 hist={}
 maxdiff=0
 for a,b in zip(words,words[1:]):
  idx=tuple(i for i,(x,y) in enumerate(zip(a,b)) if x!=y)
  hist[len(idx)]=hist.get(len(idx),0)+1;maxdiff=max(maxdiff,len(idx))
 return {"cells":len(cells),"adjacent_hamming_histogram":hist,"max_adjacent_hamming":maxdiff}
if __name__=="__main__":
 import json;print(json.dumps(adjacent_word_differences(),indent=2,sort_keys=True))
