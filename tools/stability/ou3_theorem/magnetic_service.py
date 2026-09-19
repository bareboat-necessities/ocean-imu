"""MAGNETIC SERVICE information accounting.

The service Gramian is built from the same innovation covariance and magnetic
Jacobian used by an actually applied shipping correction.  Callers may record
the resulting transported/whitened rows in MagneticEvent; callback cadence or
an attempted packet is never a substitute for this information.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence

Row2=tuple[float,float]
Matrix=tuple[tuple[float,...],...]

def _matrix(name: str,x: Sequence[Sequence[float]]) -> Matrix:
    rows=tuple(tuple(float(v) for v in row) for row in x)
    if not rows or not rows[0]: raise ValueError(f"{name} must be nonempty")
    n=len(rows[0])
    if any(len(row)!=n for row in rows): raise ValueError(f"{name} must be rectangular")
    if not all(math.isfinite(v) for row in rows for v in row):
        raise ValueError(f"{name} must be finite")
    return rows

def _matmul(a: Matrix,b: Matrix) -> Matrix:
    if len(a[0])!=len(b): raise ValueError("matrix dimension mismatch")
    bt=tuple(zip(*b))
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in bt) for row in a)

def _cholesky_spd(s: Matrix) -> Matrix:
    n=len(s)
    if any(len(row)!=n for row in s): raise ValueError("innovation covariance must be square")
    scale=max(1.0,max(abs(v) for row in s for v in row))
    tol=64.0*math.ulp(1.0)*scale
    for i in range(n):
        for j in range(n):
            if abs(s[i][j]-s[j][i])>tol:
                raise ValueError("innovation covariance must be symmetric")
    l=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            r=s[i][j]-sum(l[i][k]*l[j][k] for k in range(j))
            if i==j:
                if not r>0.0: raise ValueError("innovation covariance must be SPD")
                l[i][j]=math.sqrt(r)
            else:
                l[i][j]=r/l[j][j]
    return tuple(tuple(row) for row in l)

def _solve_lower(l: Matrix,b: Matrix) -> Matrix:
    n=len(l)
    if len(b)!=n: raise ValueError("lower-solve dimension mismatch")
    cols=len(b[0])
    x=[[0.0]*cols for _ in range(n)]
    for i in range(n):
        if l[i][i]==0.0: raise ValueError("singular lower factor")
        for j in range(cols):
            x[i][j]=(b[i][j]-sum(l[i][k]*x[k][j] for k in range(i)))/l[i][i]
    return tuple(tuple(row) for row in x)

def transported_whitened_rows(
    H_m: Sequence[Sequence[float]],
    S_m_actual: Sequence[Sequence[float]],
    Phi_from_window_root: Sequence[Sequence[float]],
    E_hb: Sequence[Sequence[float]],
) -> tuple[Row2,...]:
    """Return G=L^{-1} H_m Phi E_hb for S_m_actual=L L^T.

    S_m_actual is the literal innovation covariance presented to the shipping
    magnetic update's factorization.  H_m is that update's literal magnetic
    sensitivity, Phi is the complete preceding same-history differential, and
    E_hb injects the normalized heading/axial-gyro-bias root coordinates.
    Phi may be rectangular across a dimension-changing error-coordinate map;
    H uses destination coordinates and E uses the original root coordinates.
    Accepting dimensions does not certify any particular release map.
    Cholesky is used here only as an algebraically equivalent square root of the
    same SPD S; the shipping estimator remains authoritative and unchanged.
    """
    h=_matrix("H_m",H_m)
    s=_matrix("S_m_actual",S_m_actual)
    phi=_matrix("Phi_from_window_root",Phi_from_window_root)
    e=_matrix("E_hb",E_hb)
    if len(s)!=len(h) or len(s[0])!=len(h):
        raise ValueError("innovation covariance/measurement dimension mismatch")
    if len(h[0])!=len(phi):
        raise ValueError("H/Phi dimension mismatch")
    if len(phi[0])!=len(e):
        raise ValueError("Phi/E dimension mismatch")
    if len(e[0])!=2:
        raise ValueError("E_hb must inject exactly two information coordinates")
    raw=_matmul(_matmul(h,phi),e)
    g=_solve_lower(_cholesky_spd(s),raw)
    return tuple((row[0],row[1]) for row in g)

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


@dataclass(frozen=True)
class MagneticServiceContinuationCertificate:
    """All-time source-uniform MAGNETIC SERVICE certificate.

    A finite event replay cannot establish this object. It is a theorem-side
    continuation certificate for one persistent physical execution, analogous
    to MARINE MOTION and IMU BIAS continuation certificates.
    """
    history_id: str
    service_window_s: float
    information_floor: float
    every_window_same_history_certified: bool
    applied_event_semantics_certified: bool
    transported_actual_innovation_covariance_certified: bool
    all_time_continuation_certified: bool

    def __post_init__(self) -> None:
        if not self.history_id:
            raise ValueError("nonempty history id required")
        if not all(math.isfinite(x) and x>0 for x in
                   (self.service_window_s,self.information_floor)):
            raise ValueError("positive finite service constants required")


def continuation_admitted(cert: MagneticServiceContinuationCertificate,
                          *, required_window_s: float,
                          required_information_floor: float) -> bool:
    """Source-uniform service is an all-time contract, never a finite replay."""
    if not all(math.isfinite(x) and x>0 for x in
               (required_window_s,required_information_floor)):
        raise ValueError("positive required service constants")
    return (cert.service_window_s <= required_window_s
            and cert.information_floor >= required_information_floor
            and cert.every_window_same_history_certified
            and cert.applied_event_semantics_certified
            and cert.transported_actual_innovation_covariance_certified
            and cert.all_time_continuation_certified)
