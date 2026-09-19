"""Analytic lemmas for the single OU-III stability theorem.

H18 is a finite bridge to the implemented accelerometer-bias release.
Asymptotic contraction is required on the recurring magnetically informed A21
tail, not on the held-bias bridge.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

def finite_bridge_bound(v0: float, steps: int, gain: float, supply: float) -> float:
    """Iterate V+ <= gain*V+supply for a finite H18 bridge."""
    if not all(math.isfinite(x) for x in (v0,gain,supply)): raise ValueError("finite arguments required")
    if v0 < 0 or steps < 0 or gain < 0 or supply < 0: raise ValueError("nonnegative bridge data required")
    if steps == 0: return v0
    if gain == 1.0: return v0 + steps*supply
    return gain**steps*v0 + supply*(gain**steps-1.0)/(gain-1.0)

def release_bound_from_service(remaining_applied_updates: int, service_window_s: float,
                               guard_s: float, reference_release_s: float) -> float:
    """Finite gate-clear bound after finite magnetic-reference refinement."""
    if remaining_applied_updates < 0: raise ValueError("remaining update count must be nonnegative")
    if not all(math.isfinite(x) for x in (service_window_s,guard_s,reference_release_s)): raise ValueError("finite release data required")
    if service_window_s <= 0 or guard_s < 0 or reference_release_s < 0: raise ValueError("positive window and nonnegative times required")
    return reference_release_s + remaining_applied_updates*service_window_s + guard_s

@dataclass(frozen=True)
class A21TailPremises:
    linear_storage_ratio: float
    nonlinear_storage_lipschitz: float
    supply_gain: float
    coercivity_lower: float
    prefix_gain: float
    def __post_init__(self) -> None:
        v=(self.linear_storage_ratio,self.nonlinear_storage_lipschitz,self.supply_gain,self.coercivity_lower,self.prefix_gain)
        if not all(math.isfinite(x) for x in v): raise ValueError("finite premises required")
        if self.linear_storage_ratio < 0 or self.nonlinear_storage_lipschitz < 0: raise ValueError("nonnegative contraction data required")
        if self.supply_gain < 0 or self.coercivity_lower <= 0 or self.prefix_gain < 0: raise ValueError("valid supply/coercivity/prefix data required")

def nonlinear_tail_ratio(p: A21TailPremises) -> float:
    """Small-gain lift: rho=(sqrt(rho0)+eta)^2."""
    return (math.sqrt(p.linear_storage_ratio)+p.nonlinear_storage_lipschitz)**2

def practical_radius_bound(p: A21TailPremises, disturbance_bound: float) -> float:
    if not math.isfinite(disturbance_bound) or disturbance_bound < 0: raise ValueError("finite nonnegative disturbance bound required")
    rho=nonlinear_tail_ratio(p)
    if not rho < 1.0: raise ValueError("A21 tail is not strictly contractive")
    return math.sqrt(p.supply_gain/(p.coercivity_lower*(1.0-rho)))*disturbance_bound

def proof_route_status(*, reference_refinement_finite: bool, h18_bridge_retained: bool,
                       a21_linear_uniform: bool, a21_nonlinear_bound: bool,
                       a21_prefix_retained: bool) -> dict:
    flags=(reference_refinement_finite,h18_bridge_retained,a21_linear_uniform,a21_nonlinear_bound,a21_prefix_retained)
    if any(type(x) is not bool for x in flags): raise ValueError("literal booleans required")
    return {"h18_role":"finite bridge only","a21_role":"recurring asymptotic tail",
            "release_finite":reference_refinement_finite,"h18_bridge_closed":h18_bridge_retained,
            "a21_tail_closed":all(flags[2:]),"route_closed":all(flags)}
