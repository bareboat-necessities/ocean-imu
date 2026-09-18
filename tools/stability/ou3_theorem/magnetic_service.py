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


def _matrix(x: Sequence[Sequence[float]], rows: int | None = None,
            cols: int | None = None) -> list[list[float]]:
    a=[[float(v) for v in row] for row in x]
    if not a or not a[0]: raise ValueError("nonempty matrix required")
    width=len(a[0])
    if any(len(row)!=width for row in a): raise ValueError("ragged matrix")
    if rows is not None and len(a)!=rows: raise ValueError("matrix row mismatch")
    if cols is not None and width!=cols: raise ValueError("matrix column mismatch")
    if not all(math.isfinite(v) for row in a for v in row): raise ValueError("non-finite matrix")
    return a

def _cholesky3(s: Sequence[Sequence[float]]) -> list[list[float]]:
    a=_matrix(s,3,3)
    # Symmetrize exactly as the shipping diagnostic interpretation does before
    # whitening; a non-SPD actual innovation covariance is not informative service.
    q=[[0.5*(a[i][j]+a[j][i]) for j in range(3)] for i in range(3)]
    L=[[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(i+1):
            v=q[i][j]-sum(L[i][k]*L[j][k] for k in range(j))
            if i==j:
                if not (v>0.0 and math.isfinite(v)): raise ValueError("innovation covariance must be SPD")
                L[i][j]=math.sqrt(v)
            else:
                L[i][j]=v/L[j][j]
    return L

def transported_whitened_heading_bias_rows(
    H_m: Sequence[Sequence[float]],
    S_m: Sequence[Sequence[float]],
    Phi_E_hb: Sequence[Sequence[float]],
) -> tuple[Row2,...]:
    """Construct G=L_S^{-1} H_m Phi E_hb from actual shipping quantities.

    H_m is the 3-by-n magnetometer sensitivity used at the applied update,
    S_m is that update's actual 3-by-3 innovation covariance, and Phi_E_hb is
    the complete preceding same-history differential from normalized root
    heading/axial-gyro-bias coordinates to the pre-update n-state coordinate.
    With S_m=L_S L_S^T, G^T G=(H Phi E)^T S_m^{-1}(H Phi E).
    """
    H=_matrix(H_m)
    if len(H)!=3: raise ValueError("shipping magnetometer sensitivity must have three rows")
    n=len(H[0])
    PE=_matrix(Phi_E_hb,n,2)
    A=[[sum(H[i][k]*PE[k][j] for k in range(n)) for j in range(2)] for i in range(3)]
    L=_cholesky3(S_m)
    G=[[0.0,0.0] for _ in range(3)]
    for i in range(3):
        for j in range(2):
            rhs=A[i][j]-sum(L[i][k]*G[k][j] for k in range(i))
            G[i][j]=rhs/L[i][i]
    return tuple((row[0],row[1]) for row in G)

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
