"""MAGNETIC SERVICE information accounting."""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence

Row2=tuple[float,float]

def _row2(x: Sequence[float]) -> Row2:
    if len(x)!=2: raise ValueError("information row must have two columns")
    y=(float(x[0]),float(x[1]))
    if not all(math.isfinite(v) for v in y): raise ValueError("non-finite information row")
    return y

@dataclass(frozen=True)
class MagneticEvent:
    history_id: str
    t_s: float
    applied: bool
    gauged: bool
    finite: bool
    saturated: bool
    transported_whitened_rows: tuple[Row2,...]
    def __post_init__(self) -> None:
        if not self.history_id or not math.isfinite(self.t_s): raise ValueError("valid history and time required")
        for row in self.transported_whitened_rows: _row2(row)
    @property
    def informative_candidate(self) -> bool:
        return self.applied and self.gauged and self.finite and not self.saturated and bool(self.transported_whitened_rows)

def min_eigenvalue_2x2(a: float,b: float,c: float) -> float:
    tr=a+c; disc=math.sqrt(max(0.0,(a-c)*(a-c)+4.0*b*b)); return 0.5*(tr-disc)

def audit_window(events: Sequence[MagneticEvent],start_s: float,T_M_s: float,mu_M: float) -> dict:
    if not (math.isfinite(start_s) and math.isfinite(T_M_s) and T_M_s>0.0): raise ValueError("finite positive service window required")
    if not (math.isfinite(mu_M) and mu_M>0.0): raise ValueError("positive information floor required")
    end_s=start_s+T_M_s; chosen=[]; history=None; previous=-math.inf; attempted=0
    for e in events:
        if e.t_s<previous: raise ValueError("events must be chronological")
        previous=e.t_s
        if e.t_s<start_s or e.t_s>end_s: continue
        attempted+=1
        if history is None: history=e.history_id
        elif e.history_id!=history: raise ValueError("detached magnetic history")
        if e.informative_candidate: chosen.append(e)
    a=b=c=0.0
    for e in chosen:
        for row in e.transported_whitened_rows:
            x,y=_row2(row); a+=x*x; b+=x*y; c+=y*y
    lam=min_eigenvalue_2x2(a,b,c)
    times=[e.t_s for e in chosen]
    max_gap=None
    if times:
        bounds=[start_s,*times,end_s]; max_gap=max(v-u for u,v in zip(bounds,bounds[1:]))
    return {"start_s":start_s,"end_s":end_s,"attempted_events":attempted,
            "actually_applied_informative_events":len(chosen),
            "information_gramian":[[a,b],[b,c]],"information_min_eigenvalue":lam,
            "information_pass":lam>=mu_M,"max_usable_event_gap_s_diagnostic_only":max_gap,
            "gap_alone_implies_service":False,"finite_window_audit_only":True}
