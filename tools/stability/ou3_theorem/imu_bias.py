"""IMU BIAS same-history contract and literal estimator operations."""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence

Vec3=tuple[float,float,float]

def vec(x: Sequence[float]) -> Vec3:
    if len(x)!=3: raise ValueError("expected three-vector")
    y=tuple(float(v) for v in x)
    if not all(math.isfinite(v) for v in y): raise ValueError("non-finite vector")
    return y  # type: ignore[return-value]

def norm(x: Sequence[float]) -> float:
    y=vec(x); return math.sqrt(sum(v*v for v in y))

def add(a: Sequence[float],b: Sequence[float]) -> Vec3:
    x,y=vec(a),vec(b); return tuple(u+v for u,v in zip(x,y))  # type: ignore[return-value]

def sub(a: Sequence[float],b: Sequence[float]) -> Vec3:
    x,y=vec(a),vec(b); return tuple(u-v for u,v in zip(x,y))  # type: ignore[return-value]

def scale(c: float,x: Sequence[float]) -> Vec3:
    y=vec(x); return tuple(c*v for v in y)  # type: ignore[return-value]

@dataclass(frozen=True)
class BiasLimits:
    B_a_mps2: float
    D_a_mps3: float
    B_g_rad_s: float
    D_g_rad_s2: float
    def __post_init__(self) -> None:
        vals=(self.B_a_mps2,self.D_a_mps3,self.B_g_rad_s,self.D_g_rad_s2)
        if not all(math.isfinite(v) and v>=0.0 for v in vals):
            raise ValueError("bias limits must be finite and nonnegative")

@dataclass(frozen=True)
class BiasSample:
    history_id: str
    t_s: float
    b_a_mps2: Vec3
    b_g_rad_s: Vec3
    def __post_init__(self) -> None:
        if not self.history_id or not math.isfinite(self.t_s): raise ValueError("valid history and time required")
        vec(self.b_a_mps2); vec(self.b_g_rad_s)

def successor_allowed(current: Sequence[float],successor: Sequence[float],dt_s: float,
                      amplitude_limit: float,rate_limit: float,*,tolerance: float=0.0) -> bool:
    if not (math.isfinite(dt_s) and dt_s>0.0): return False
    if min(amplitude_limit,rate_limit,tolerance)<0.0: return False
    return (norm(current)<=amplitude_limit+tolerance and
            norm(successor)<=amplitude_limit+tolerance and
            norm(sub(successor,current))<=rate_limit*dt_s+tolerance)

def audit_bias_trace(samples: Sequence[BiasSample],limits: BiasLimits,*,tolerance: float=0.0) -> dict:
    if not samples: raise ValueError("empty bias history")
    failures=[]; hid=samples[0].history_id
    for i,s in enumerate(samples):
        if s.history_id!=hid: failures.append(f"sample {i}: detached history")
        if norm(s.b_a_mps2)>limits.B_a_mps2+tolerance: failures.append(f"sample {i}: accelerometer amplitude")
        if norm(s.b_g_rad_s)>limits.B_g_rad_s+tolerance: failures.append(f"sample {i}: gyro amplitude")
        if i:
            p=samples[i-1]; dt=s.t_s-p.t_s
            if not successor_allowed(p.b_a_mps2,s.b_a_mps2,dt,limits.B_a_mps2,limits.D_a_mps3,tolerance=tolerance):
                failures.append(f"sample {i}: accelerometer predecessor/rate")
            if not successor_allowed(p.b_g_rad_s,s.b_g_rad_s,dt,limits.B_g_rad_s,limits.D_g_rad_s2,tolerance=tolerance):
                failures.append(f"sample {i}: gyro predecessor/rate")
    return {"finite_prefix_pass":not failures,"failures":failures,"history_id":hid,
            "independent_successor_boxes_allowed":False,"physical_bias_reset_at_estimator_release":False}

def estimator_prediction_factor(mode: str,phi_ou: float) -> float:
    """Literal accelerometer-bias estimate prediction factor for the shipping mode."""
    if not math.isfinite(phi_ou) or not (0.0 < phi_ou <= 1.0):
        raise ValueError("phi_ou must be finite in (0,1]")
    if mode=="H18": return 1.0
    if mode=="A21": return phi_ou
    raise ValueError("mode must be H18 or A21")

@dataclass(frozen=True)
class BiasPredictionRelation:
    """One physical bias law with a mode-dependent estimator prediction."""
    mode: str
    phi_e: float
    e_b_plus: Vec3
    b_true_k: Vec3
    w_k: Vec3
    model_mismatch: Vec3
    e_b_minus_next: Vec3

def accel_prediction_relation(mode: str,phi_ou: float,e_b_plus: Sequence[float],
                              b_true_k: Sequence[float],w_k: Sequence[float]) -> BiasPredictionRelation:
    """First-class relation e^- = phi_e e^+ + (1-phi_e)b + w.

    H18 uses phi_e=1. A21 uses the literal shipping OU factor. Physical truth
    is the same bounded/rate-bounded bias history in both modes.
    """
    phi_e=estimator_prediction_factor(mode,phi_ou)
    e_plus=vec(e_b_plus); b_true=vec(b_true_k); w=vec(w_k)
    mismatch=scale(1.0-phi_e,b_true)
    e_minus=add(add(scale(phi_e,e_plus),mismatch),w)
    return BiasPredictionRelation(mode,phi_e,e_plus,b_true,w,mismatch,e_minus)

def accel_prediction_error(mode: str,phi_ou: float,e_b_plus: Sequence[float],
                           b_true_k: Sequence[float],w_k: Sequence[float]) -> Vec3:
    return accel_prediction_relation(mode,phi_ou,e_b_plus,b_true_k,w_k).e_b_minus_next

def gyro_prediction_error(e_b_plus: Sequence[float],w_k: Sequence[float]) -> Vec3:
    """Shipping gyro-bias mean predictor is identity."""
    return add(e_b_plus,w_k)

def correction_error(e_minus: Sequence[float],estimated_bias_increment: Sequence[float]) -> Vec3:
    """Physical truth is unchanged by a Kalman correction: e_corr=e_minus-delta_bhat."""
    return sub(e_minus,estimated_bias_increment)

def release_error(e_before: Sequence[float]) -> Vec3:
    """H18-to-A21 release changes update permission/covariance, not bias truth/state."""
    return vec(e_before)

def project_ball(x: Sequence[float],radius: float) -> Vec3:
    if not (math.isfinite(radius) and radius>=0.0): raise ValueError("finite nonnegative radius required")
    y=vec(x); n=norm(y)
    return y if n<=radius or n==0.0 else scale(radius/n,y)

def projection_defect(estimated_bias_corrected: Sequence[float],radius: float) -> Vec3:
    """r_proj=bhat_corr-Proj_R(bhat_corr), matching the shipping estimate projection."""
    y=vec(estimated_bias_corrected)
    return sub(y,project_ball(y,radius))

def projected_error(b_true: Sequence[float],estimated_bias_corrected: Sequence[float],
                    radius: float) -> Vec3:
    """Post-projection error b_true-Proj_R(bhat_corr)."""
    return sub(b_true,project_ball(estimated_bias_corrected,radius))
