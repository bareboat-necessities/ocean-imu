"""Scheduler branch enumeration for a fixed one-second/20-second word.

For fixed dt and T_S, periodic_update_due depends on initial elapsed only
through threshold crossings. Exact branch endpoints are residues of n*dt modulo
T_S. This module enumerates open cells and representatives; interval promotion
then needs only endpoint/one-sided treatment of the finite threshold set.
"""
from __future__ import annotations
import math

def branch_points(period,dt,steps):
    pts={0.0,float(period)}
    for n in range(1,steps+1):
        # initial e that makes a due event exactly at step n modulo period
        x=(period-(n*dt)%period)%period
        if 0<x<period:pts.add(x)
    return sorted(pts)

def representatives(period,dt,steps):
    p=branch_points(period,dt,steps)
    return [(p[i],p[i+1],(p[i]+p[i+1])/2) for i in range(len(p)-1)]

def due_word(elapsed,period,dt,steps):
    out=[]
    for _ in range(steps):
        total=elapsed+dt
        if total < period:
            elapsed=total;out.append(False)
        else:
            elapsed=math.fmod(total,period);out.append(True)
    return tuple(out),elapsed

def certificate(period=.1363605111837387,dt=.005,steps=4000):
    cells=representatives(period,dt,steps)
    words={due_word(m,period,dt,steps)[0] for _,_,m in cells}
    return {"qualification":"OU3_SCHEDULER_BRANCH_ENUMERATION_V1",
            "period":period,"dt":dt,"steps":steps,
            "threshold_points":len(branch_points(period,dt,steps)),
            "open_cells":len(cells),"distinct_due_words":len(words),
            "finite_branch_reduction":True,
            "endpoint_float32_tolerance_handling_verified":False,
            "theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
