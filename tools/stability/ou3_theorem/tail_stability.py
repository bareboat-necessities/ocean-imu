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
    """Sufficient conditions for the deployed unweighted MagAutoTuner."""
    min_samples: int
    min_window_s: float
    service_window_s: float
    true_field_norm_lower: float
    true_field_norm_upper: float
    measurement_residual_norm: float
    true_horizontal_lower: float
    tilt_error_upper_rad: float
    max_norm_ratio_from_mean: float
    min_horizontal_fraction: float
    def __post_init__(self) -> None:
        if self.min_samples < 1: raise ValueError("positive sample count required")
        vals=(self.min_window_s,self.service_window_s,
              self.true_field_norm_lower,self.true_field_norm_upper,
              self.measurement_residual_norm,self.true_horizontal_lower,
              self.tilt_error_upper_rad,self.max_norm_ratio_from_mean,
              self.min_horizontal_fraction)
        if not all(math.isfinite(x) for x in vals): raise ValueError("finite refinement premises required")
        if self.min_window_s < 0 or self.service_window_s <= 0:
            raise ValueError("positive service-window timing required")
        if self.true_field_norm_lower <= self.measurement_residual_norm or self.true_field_norm_upper < self.true_field_norm_lower:
            raise ValueError("field must dominate residual")
        if self.measurement_residual_norm < 0 or self.true_horizontal_lower <= 0 or self.tilt_error_upper_rad < 0:
            raise ValueError("valid field/residual/tilt bounds required")
        if self.max_norm_ratio_from_mean < 0 or not 0 < self.min_horizontal_fraction < 1:
            raise ValueError("valid tuner gates required")

def refinement_gate_margins(p: RefinementPremises) -> dict:
    """Derive literal MagAutoTuner gates from physical field/residual bounds.

    Rotation preserves the true field norm. With ||r||<=R, every corrected
    sample norm lies in [B_min-R,B_max+R]. For one physical field magnitude the
    sharper running-norm variation is 2R/(B_min-R), independent of attitude.
    A tilt-frame error eps changes a vector by at most 2 B_max sin(eps/2);
    adding the measurement residual gives a conservative horizontal-mean loss.
    """
    norm_ratio=2.0*p.measurement_residual_norm/(p.true_field_norm_lower-p.measurement_residual_norm)
    rotation_loss=2.0*p.true_field_norm_upper*math.sin(min(math.pi,p.tilt_error_upper_rad)/2.0)
    horizontal_lower=p.true_horizontal_lower-rotation_loss-p.measurement_residual_norm
    sample_norm_upper=p.true_field_norm_upper+p.measurement_residual_norm
    horizontal_fraction=horizontal_lower/sample_norm_upper
    return {"norm_ratio_upper":norm_ratio,
            "norm_ratio_margin":p.max_norm_ratio_from_mean-norm_ratio,
            "horizontal_mean_lower":horizontal_lower,
            "horizontal_fraction_lower":horizontal_fraction,
            "horizontal_fraction_margin":horizontal_fraction-p.min_horizontal_fraction}

def refinement_tilt_limit_rad(p: RefinementPremises) -> float:
    """Largest tilt-frame error allowed by the literal horizontal gate."""
    required=p.min_horizontal_fraction*(p.true_field_norm_upper+p.measurement_residual_norm)
    budget=p.true_horizontal_lower-p.measurement_residual_norm-required
    if budget <= 0.0: return 0.0
    x=min(1.0,budget/(2.0*p.true_field_norm_upper))
    return 2.0*math.asin(x)

def refinement_sample_gate_uniform(p: RefinementPremises) -> bool:
    m=refinement_gate_margins(p)
    return m["norm_ratio_margin"] >= 0.0 and m["horizontal_fraction_margin"] >= 0.0

def refinement_completion_bound(start_s: float, p: RefinementPremises) -> float:
    """Completion bound from recurring MAGNETIC SERVICE.

    Once the physical margins make every service event acceptable to the
    unweighted tuner, at least one such event occurs per T_M. The count reaches
    N in at most N*T_M. Because accepted-window time telescopes between accepted
    callbacks, the elapsed-window gate reaches W in at most W+T_M.
    """
    if not math.isfinite(start_s) or start_s < 0: raise ValueError("finite nonnegative start required")
    if not refinement_sample_gate_uniform(p):
        raise ValueError("physical bounds do not guarantee MagAutoTuner acceptance")
    wait=max(p.min_samples*p.service_window_s,p.min_window_s+p.service_window_s)
    return start_s+wait
