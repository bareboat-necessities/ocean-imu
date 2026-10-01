"""Proof-side literal magnetic event-box exporter.

This is intentionally separate from shipping estimator behavior.  It carries
the complete root differential Phi through a supplied literal operation trace
and captures pre-correction H_m, S_m_actual, Phi for actually-applied magnetic
events.  Every accepted/rejected/gate stratum must be supplied separately;
this module never hulls different event decisions.
"""
from __future__ import annotations
from dataclasses import dataclass
from .interval_riccati_21 import IMat,matmul
from .rank_loss_interval_factor import eye


@dataclass(frozen=True)
class PredictionBox:
    F: IMat

@dataclass(frozen=True)
class CorrectionBox:
    A: IMat                 # literal homogeneous mean map I-KH
    sensor: str
    applied: bool
    H: IMat|None=None
    S_actual: IMat|None=None
    t_lo: float=0.0
    t_hi: float=0.0

@dataclass(frozen=True)
class ResetBox:
    G: IMat

@dataclass(frozen=True)
class SyncBox:
    """Covariance-only shipping sync; homogeneous mean differential is identity."""
    pass

@dataclass(frozen=True)
class MagneticLiteralEventBox:
    t_lo: float
    t_hi: float
    H_m: IMat
    S_m_actual: IMat
    Phi_from_window_root: IMat


def export_magnetic_event_boxes(root_dim:int,ops) -> tuple[MagneticLiteralEventBox,...]:
    """Carry Phi and capture literal PRE-correction magnetic event boxes."""
    Phi=eye(root_dim);out=[]
    for op in ops:
        if isinstance(op,PredictionBox):
            Phi=matmul(op.F,Phi)
        elif isinstance(op,ResetBox):
            Phi=matmul(op.G,Phi)
        elif isinstance(op,SyncBox):
            continue
        elif isinstance(op,CorrectionBox):
            if op.sensor=="mag" and op.applied:
                if op.H is None or op.S_actual is None:
                    raise ValueError("applied magnetic event needs H and actual S")
                if op.H.shape[1]!=Phi.shape[0]:
                    raise ValueError("magnetic H/Phi dimension mismatch")
                out.append(MagneticLiteralEventBox(op.t_lo,op.t_hi,op.H,op.S_actual,Phi))
            # Correction differential applies AFTER the captured pre-correction row.
            Phi=matmul(op.A,Phi)
        else:
            raise TypeError(f"unknown operation {type(op)!r}")
    return tuple(out)


def event_tuples(events):
    """Adapter consumed by magnetic_nuisance_interval.service_schur_information."""
    return tuple((e.H_m,e.S_m_actual,e.Phi_from_window_root) for e in events)


def validate_one_second_stratum(events,window_lo=0.0,window_hi=1.0):
    if not events:
        return {"verified":False,"reason":"no actually-applied magnetic event boxes"}
    if any(e.t_lo<window_lo or e.t_hi>window_hi or e.t_lo>e.t_hi for e in events):
        return {"verified":False,"reason":"event time interval leaves service window"}
    # Event decisions are fixed by the caller's closed combinatorial stratum.
    return {"verified":True,"event_count":len(events),
            "accepted_event_stratum_fixed":True,
            "window":[window_lo,window_hi]}
