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


@dataclass(frozen=True)
class RefinementPremises:
    """Sufficient conditions for MagAutoTuner refinement to finish."""
    min_samples: int
    min_window_s: float
    usable_sample_gap_s: float
    field_norm_lower: float
    field_norm_upper: float
    max_norm_ratio_from_mean: float
    horizontal_mean_lower: float
    min_horizontal_fraction: float
    def __post_init__(self) -> None:
        if self.min_samples < 1: raise ValueError("positive sample count required")
        vals=(self.min_window_s,self.usable_sample_gap_s,self.field_norm_lower,
              self.field_norm_upper,self.max_norm_ratio_from_mean,
              self.horizontal_mean_lower,self.min_horizontal_fraction)
        if not all(math.isfinite(x) for x in vals): raise ValueError("finite refinement premises required")
        if self.min_window_s < 0 or self.usable_sample_gap_s <= 0 or self.field_norm_lower <= 0:
            raise ValueError("positive refinement timing/field bounds required")
        if self.field_norm_upper < self.field_norm_lower or self.max_norm_ratio_from_mean < 0:
            raise ValueError("ordered field bounds required")
        if self.horizontal_mean_lower <= 0 or not 0 < self.min_horizontal_fraction < 1:
            raise ValueError("positive horizontal field requirement")

def refinement_sample_gate_uniform(p: RefinementPremises) -> bool:
    """Sufficient all-sample conditions for the unweighted MagAutoTuner gate.

    The running accepted mean norm remains in [field_norm_lower,field_norm_upper].
    Hence every next norm differs from that mean by at most
    (upper-lower)/lower.  The horizontal mean is also nondegenerate.
    """
    norm_ratio=(p.field_norm_upper-p.field_norm_lower)/p.field_norm_lower
    horizontal_fraction=p.horizontal_mean_lower/p.field_norm_upper
    return (norm_ratio <= p.max_norm_ratio_from_mean and
            horizontal_fraction >= p.min_horizontal_fraction)

def refinement_completion_bound(start_s: float, p: RefinementPremises) -> float:
    """Finite refinement time once recurring usable samples satisfy the gate."""
    if not math.isfinite(start_s) or start_s < 0: raise ValueError("finite nonnegative start required")
    if not refinement_sample_gate_uniform(p):
        raise ValueError("premises do not guarantee MagAutoTuner sample acceptance")
    # One usable sample every gap. The accepted-window clock advances by the
    # actual intersample time, so min_window and min_samples are both covered.
    return start_s + max(p.min_window_s, p.min_samples*p.usable_sample_gap_s)
