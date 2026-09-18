"""MARINE MOTION contract helpers.

Finite sampled traces can falsify pointwise or same-history constraints, but a
finite replay cannot prove the all-time bounded primitive.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Iterable, Sequence

Vec3 = tuple[float, float, float]

def _vec3(x: Sequence[float]) -> Vec3:
    if len(x) != 3:
        raise ValueError("expected three-vector")
    y = tuple(float(v) for v in x)
    if not all(math.isfinite(v) for v in y):
        raise ValueError("non-finite vector")
    return y  # type: ignore[return-value]

def norm(x: Sequence[float]) -> float:
    a = _vec3(x)
    return math.sqrt(sum(v * v for v in a))

def sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    x, y = _vec3(a), _vec3(b)
    return tuple(u - v for u, v in zip(x, y))  # type: ignore[return-value]

def add(a: Sequence[float], b: Sequence[float]) -> Vec3:
    x, y = _vec3(a), _vec3(b)
    return tuple(u + v for u, v in zip(x, y))  # type: ignore[return-value]

def scale(c: float, a: Sequence[float]) -> Vec3:
    x = _vec3(a)
    return tuple(c * v for v in x)  # type: ignore[return-value]

@dataclass(frozen=True)
class MarineLimits:
    p_max_m: float
    v_max_mps: float
    a_max_mps2: float
    omega_max_rad_s: float
    p_ac_max_m_s: float
    def __post_init__(self) -> None:
        vals=(self.p_max_m,self.v_max_mps,self.a_max_mps2,self.omega_max_rad_s,self.p_ac_max_m_s)
        if not all(math.isfinite(v) and v > 0.0 for v in vals):
            raise ValueError("marine limits must be positive and finite")

@dataclass(frozen=True)
class MarineContinuationCertificate:
    """All-time theorem-side MARINE MOTION certificate.

    These fields are proof obligations for one continuous physical history;
    finite replay data may falsify them but cannot establish all-time
    continuation by itself.
    """
    history_id: str
    p_norm_upper_m: float
    v_norm_upper_mps: float
    a_norm_upper_mps2: float
    omega_norm_upper_rad_s: float
    primitive_span_upper_m_s: float
    translation_kinematics_certified: bool
    attitude_rate_kinematics_certified: bool
    persistent_primitive_state_certified: bool
    reference_acceleration_accounted: bool
    all_time_continuation_certified: bool

    def __post_init__(self) -> None:
        if not self.history_id:
            raise ValueError("nonempty history_id required")
        vals=(
            self.p_norm_upper_m,self.v_norm_upper_mps,self.a_norm_upper_mps2,
            self.omega_norm_upper_rad_s,self.primitive_span_upper_m_s,
        )
        if not all(math.isfinite(v) and v>=0.0 for v in vals):
            raise ValueError("certificate bounds must be finite and nonnegative")


def continuation_admitted(cert: MarineContinuationCertificate,
                          limits: MarineLimits) -> bool:
    """Require one persistent motion/attitude continuation and bounded q."""
    return (
        cert.p_norm_upper_m <= limits.p_max_m
        and cert.v_norm_upper_mps <= limits.v_max_mps
        and cert.a_norm_upper_mps2 <= limits.a_max_mps2
        and cert.omega_norm_upper_rad_s <= limits.omega_max_rad_s
        and cert.primitive_span_upper_m_s <= limits.p_ac_max_m_s
        and cert.translation_kinematics_certified
        and cert.attitude_rate_kinematics_certified
        and cert.persistent_primitive_state_certified
        and cert.reference_acceleration_accounted
        and cert.all_time_continuation_certified
    )


def quiet_water_continuation_certificate(history_id: str="quiet") -> MarineContinuationCertificate:
    return MarineContinuationCertificate(
        history_id=history_id,
        p_norm_upper_m=0.0,
        v_norm_upper_mps=0.0,
        a_norm_upper_mps2=0.0,
        omega_norm_upper_rad_s=0.0,
        primitive_span_upper_m_s=0.0,
        translation_kinematics_certified=True,
        attitude_rate_kinematics_certified=True,
        persistent_primitive_state_certified=True,
        reference_acceleration_accounted=True,
        all_time_continuation_certified=True,
    )


@dataclass(frozen=True)
class MarineSample:
    history_id: str
    t_s: float
    p_m: Vec3
    v_mps: Vec3
    a_mps2: Vec3
    omega_rad_s: Vec3
    q_m_s: Vec3 | None = None
    def __post_init__(self) -> None:
        if not self.history_id or not math.isfinite(self.t_s):
            raise ValueError("valid history and time required")
        for value in (self.p_m,self.v_mps,self.a_mps2,self.omega_rad_s):
            _vec3(value)
        if self.q_m_s is not None:
            _vec3(self.q_m_s)

def constant_displacement_continuation_admitted(p0_m: Sequence[float]) -> bool:
    return norm(p0_m) == 0.0

def quiet_water_admitted() -> bool:
    return constant_displacement_continuation_admitted((0.0,0.0,0.0))

def primitive_span(q_values: Iterable[Sequence[float]]) -> float:
    q=[_vec3(v) for v in q_values]
    if not q:
        raise ValueError("at least one primitive value required")
    return max((norm(sub(a,b)) for i,a in enumerate(q) for b in q[i+1:]),default=0.0)

def audit_sampled_trace(samples: Sequence[MarineSample], limits: MarineLimits, *,
                        kinematic_tolerance: float=0.0,
                        primitive_tolerance: float=0.0) -> dict:
    if not samples:
        raise ValueError("empty trace")
    if kinematic_tolerance < 0.0 or primitive_tolerance < 0.0:
        raise ValueError("tolerances must be nonnegative")
    hid=samples[0].history_id
    failures=[]
    q_values=[]
    for i,s in enumerate(samples):
        if s.history_id != hid: failures.append(f"sample {i}: detached history")
        if norm(s.p_m) > limits.p_max_m: failures.append(f"sample {i}: displacement bound")
        if norm(s.v_mps) > limits.v_max_mps: failures.append(f"sample {i}: velocity bound")
        if norm(s.a_mps2) > limits.a_max_mps2: failures.append(f"sample {i}: acceleration bound")
        if norm(s.omega_rad_s) > limits.omega_max_rad_s: failures.append(f"sample {i}: angular-rate bound")
        if s.q_m_s is not None: q_values.append(_vec3(s.q_m_s))
        if i == 0: continue
        p=samples[i-1]
        dt=s.t_s-p.t_s
        if not (math.isfinite(dt) and dt > 0.0):
            failures.append(f"sample {i}: non-increasing time")
            continue
        p_pred=add(p.p_m,scale(0.5*dt,add(p.v_mps,s.v_mps)))
        v_pred=add(p.v_mps,scale(0.5*dt,add(p.a_mps2,s.a_mps2)))
        if norm(sub(s.p_m,p_pred)) > kinematic_tolerance: failures.append(f"sample {i}: p/v same-history enclosure")
        if norm(sub(s.v_mps,v_pred)) > kinematic_tolerance: failures.append(f"sample {i}: v/a same-history enclosure")
        if p.q_m_s is not None and s.q_m_s is not None:
            q_pred=add(p.q_m_s,scale(0.5*dt,add(p.p_m,s.p_m)))
            if norm(sub(s.q_m_s,q_pred)) > primitive_tolerance: failures.append(f"sample {i}: q/p same-history enclosure")
    if q_values and primitive_span(q_values) > limits.p_ac_max_m_s + primitive_tolerance:
        failures.append("sampled primitive span")
    return {
        "finite_prefix_pass":not failures,
        "failures":failures,
        "one_history_id":hid,
        "all_time_primitive_required":True,
        "all_time_membership_certified_by_finite_trace":False,
        "constant_nonzero_displacement_continuation_admitted":False,
        "quiet_water_admitted":True,
    }
