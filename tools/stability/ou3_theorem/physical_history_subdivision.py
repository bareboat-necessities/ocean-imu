"""Physical-input-only subdivision driver for the causal release leaf."""
from __future__ import annotations
import re

PHYSICAL_PREFIXES=("physical_p_","physical_v_","slow_accel_","slow_gyro_",
                   "fast_accel_primitive_","fast_gyro_primitive_",
                   "gravity_dir_","mag_field_")
GENERATED_TOKENS=("tau","sigma","R_S","T_S","mahony","guard","weight","covariance","gain","scheduler")

def is_physical_coordinate(name):
 return name.startswith(PHYSICAL_PREFIXES) and not any(x in name for x in GENERATED_TOKENS)

def physical_names(cell):
 return [k for k,_ in cell.coordinates if is_physical_coordinate(k)]

def guard_floor_sensitivity(cell,failure_time):
 """Prioritize coordinates causally capable of affecting the failed prefix."""
 scores={}
 for name,iv in cell.coordinates:
  if not is_physical_coordinate(name):continue
  m=re.search(r"_([0-9]+(?:\.[0-9]+)?)$",name)
  t=float(m.group(1)) if m else 0.
  if t>failure_time+30.:continue
  base=iv.width
  if name.startswith(("gravity_dir_","slow_accel_","fast_accel_primitive_")):base*=4.
  elif name.startswith(("physical_v_","physical_p_")):base*=2.
  scores[name]=base/(1.+abs(t-failure_time))
 return scores

def split_guard_failure(cell,failure_time):
 names=physical_names(cell)
 if not names:raise ArithmeticError("no physical HistoryCell coordinate available")
 scores=guard_floor_sensitivity(cell,failure_time)
 name=max(names,key=lambda k:(scores.get(k,0.),dict(cell.coordinates)[k].width))
 a,b=cell.split(name)
 return {"coordinate":name,"left":a,"right":b,"generated_split":False,
         "failure_time":failure_time}

def assert_no_generated_coordinates(cell):
 bad=[k for k,_ in cell.coordinates if not is_physical_coordinate(k)]
 if bad:raise AssertionError("nonphysical/generated root coordinates: "+",".join(bad))
 return True
