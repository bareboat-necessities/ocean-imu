"""Dependency-preserving reachable-history enclosure for OU-III shaped supply.

Only shared physical/sensor history coordinates may be subdivided.  Generated
Mahony/tuner/covariance/gain/scheduler quantities are outputs of one propagator.
This module is proof-side and changes no shipping behavior.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Any
import math

FORBIDDEN_SPLIT={"tau","sigma_aw","R_S","T_S","covariance","gain","mahony",
                 "frequency","variance","scheduler_phase","mag_reference","ba_gate"}

@dataclass(frozen=True)
class HistoryInterval:
    lo: float
    hi: float
    def __post_init__(self):
        if not (math.isfinite(self.lo) and math.isfinite(self.hi) and self.lo<=self.hi):
            raise ValueError("finite ordered history interval required")
    @property
    def width(self): return self.hi-self.lo

@dataclass(frozen=True)
class HistoryCell:
    """One shared causal input cell; generated outputs are deliberately absent."""
    coordinates: tuple[tuple[str,HistoryInterval], ...]
    prefix_token: str
    def __post_init__(self):
        names=[k for k,_ in self.coordinates]
        if len(names)!=len(set(names)): raise ValueError("duplicate history coordinate")
        bad=FORBIDDEN_SPLIT.intersection(names)
        if bad: raise ValueError("generated outputs cannot be history coordinates: "+",".join(sorted(bad)))
    def split(self,name:str):
        d=dict(self.coordinates)
        if name in FORBIDDEN_SPLIT: raise ValueError("cannot split generated output")
        if name not in d: raise KeyError(name)
        x=d[name]; m=(x.lo+x.hi)/2
        def repl(v):
            return HistoryCell(tuple((k,v if k==name else q) for k,q in self.coordinates),self.prefix_token)
        return repl(HistoryInterval(x.lo,m)),repl(HistoryInterval(m,x.hi))

@dataclass(frozen=True)
class JointPairBound:
    """Joint enclosure of dissipation and supply from ONE propagated cell."""
    dissipation_lower: float
    dissipation_upper: float
    supply_lower: float
    supply_upper: float
    dependency_token: str
    zero_set_excluded: bool
    def quotient_upper(self):
        if not self.zero_set_excluded or not self.dissipation_lower>0: return math.inf
        return max(0.0,self.supply_upper)/self.dissipation_lower

@dataclass(frozen=True)
class PropagatedCell:
    root: HistoryCell
    generated_state: Any
    pair: JointPairBound
    prefix_pairs: tuple[JointPairBound,...]

def propagate(cell:HistoryCell, shipping_propagator:Callable[[HistoryCell],PropagatedCell]):
    out=shipping_propagator(cell)
    if out.root!=cell: raise ArithmeticError("propagator changed shared history root")
    tokens={out.pair.dependency_token,*[p.dependency_token for p in out.prefix_pairs]}
    if tokens!={cell.prefix_token}: raise ArithmeticError("source/loss dependency token separated")
    return out

def certify(leaves:tuple[PropagatedCell,...], entry_budget:float):
    if not leaves: raise ValueError("nonempty reachable cover required")
    q=max(x.pair.quotient_upper() for x in leaves)
    prefix_ok=all(p.zero_set_excluded and p.dissipation_lower>0
                  for x in leaves for p in x.prefix_pairs)
    return {"source_uniform_verified":math.isfinite(q),
            "linked_supply_to_dissipation_upper":q,
            "entry_budget":entry_budget,
            "entry_budget_passed":math.isfinite(q) and q<entry_budget,
            "every_prefix_positive_dissipation":prefix_ok,
            "history_cells":len(leaves),
            "generated_outputs_independently_boxed":False}

def rigorous_cell_from_interval_joint(root:HistoryCell,iq,zero_set_excluded:bool):
    """Fail-closed scalar enclosure from one dependency-preserving interval joint form.

    This deliberately does not divide by a global unrelated eigenvalue.  A
    later kernel-aware generalized solver must provide the positive homogeneous
    restriction; until then the cell cannot promote.
    """
    from tools.stability.ou3_theorem.interval_riccati import symmetric_interval_gershgorin
    alo,ahi=symmetric_interval_gershgorin(iq.A.mid,iq.A.rad)
    # Supply includes cross/source blocks and cannot be bounded independently
    # without a source-domain metric. Fail closed until the SAME cell supplies it.
    return {"dependency_token":root.prefix_token,"homogeneous_action_lower":alo,
            "homogeneous_action_upper":ahi,"zero_set_excluded":zero_set_excluded,
            "source_metric_attached":False,"constructive_quotient_verified":False,
            "reason":"same-history source-domain metric not yet attached"}
