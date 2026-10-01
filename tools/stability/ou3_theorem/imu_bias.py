"""IMU BIAS: two-timescale same-history contract and literal estimator algebra.

Every b_true/b_true_k in the estimator algebra below denotes the SLOW component.
The shipping estimated-bias state is unchanged; FAST measurement error remains
an actual sensor forcing on that SAME execution, never a second slow state.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence

from .imu_temporal import FastWindow, slow_change, audit_fast_prefix

Vec3=tuple[float,float,float]

def vec(x: Sequence[float]) -> Vec3:
    if len(x)!=3: raise ValueError("expected three-vector")
    y=tuple(float(v) for v in x)
    if not all(math.isfinite(v) for v in y): raise ValueError("non-finite vector")
    return y  # type: ignore[return-value]

def norm(x: Sequence[float]) -> float:
    y=vec(x); return math.hypot(*y)

def add(a: Sequence[float],b: Sequence[float]) -> Vec3:
    x,y=vec(a),vec(b); return tuple(u+v for u,v in zip(x,y))  # type: ignore[return-value]

def sub(a: Sequence[float],b: Sequence[float]) -> Vec3:
    x,y=vec(a),vec(b); return tuple(u-v for u,v in zip(x,y))  # type: ignore[return-value]

def scale(c: float,x: Sequence[float]) -> Vec3:
    y=vec(x); return tuple(c*v for v in y)  # type: ignore[return-value]

@dataclass(frozen=True)
class BiasLimits:
    """One SLOW + FAST profile; missing temporal parameters remain OPEN."""
    B_a_s_mps2: float
    D_a_s_mps3: float
    B_g_s_rad_s: float
    D_g_s_rad_s2: float
    B_a_f_mps2: float
    B_g_f_rad_s: float
    accel_fast_window: FastWindow | None = None
    gyro_fast_window: FastWindow | None = None

    def __post_init__(self) -> None:
        vals=(self.B_a_s_mps2,self.D_a_s_mps3,self.B_g_s_rad_s,
              self.D_g_s_rad_s2,self.B_a_f_mps2,self.B_g_f_rad_s)
        if not all(math.isfinite(v) and v>=0.0 for v in vals):
            raise ValueError("bias limits must be finite and nonnegative")
        for w,b in ((self.accel_fast_window,self.B_a_f_mps2),
                    (self.gyro_fast_window,self.B_g_f_rad_s)):
            if w is not None:
                w.validate_amplitude(b)

    @classmethod
    def from_constants(cls):
        import json
        from pathlib import Path
        c=json.loads(Path(__file__).with_name("constants.json").read_text())["imu_bias"]
        def window(h,c):
            if (h is None)!=(c is None):
                raise ValueError("fast horizon and accumulation cap must be qualified together")
            return None if h is None else FastWindow(h,c)
        return cls(c["B_a_s_mps2"],c["D_a_s_mps3"],c["B_g_s_rad_s"],c["D_g_s_rad_s2"],
                   c["B_a_f_mps2"],c["B_g_f_rad_s"],
                   window(c["fast_accel_horizon_s"],c["fast_accel_accumulation_cap_mps"]),
                   window(c["fast_gyro_horizon_s"],c["fast_gyro_accumulation_cap_rad"]))

    @property
    def temporal_parameters_present(self) -> bool:
        return self.accel_fast_window is not None and self.gyro_fast_window is not None


@dataclass(frozen=True)
class BiasContinuationCertificate:
    """External all-time evidence for ONE carried decomposition, BOTH sensors.

    Calibration gates and a finite capture cannot manufacture this certificate.
    The window fields certify every placement, including cross-word windows.
    """
    history_id: str
    accel_slow_norm_upper_mps2: float
    accel_slow_rate_norm_upper_mps3: float
    gyro_slow_norm_upper_rad_s: float
    gyro_slow_rate_norm_upper_rad_s2: float
    predecessor_continuity_certified: bool
    calibration_scope_qualified: bool
    all_time_continuation_certified: bool
    accel_fast_norm_upper_mps2: float
    gyro_fast_norm_upper_rad_s: float
    accel_fast_window: FastWindow | None = None
    gyro_fast_window: FastWindow | None = None
    decomposition_id: str = ""
    temporal_evidence_id: str = ""
    delivered_sample_scope_certified: bool = False
    all_placed_windows_certified: bool = False

    def __post_init__(self) -> None:
        if not self.history_id:
            raise ValueError("nonempty history_id required")
        vals=(self.accel_slow_norm_upper_mps2,self.accel_slow_rate_norm_upper_mps3,
              self.gyro_slow_norm_upper_rad_s,self.gyro_slow_rate_norm_upper_rad_s2,
              self.accel_fast_norm_upper_mps2,self.gyro_fast_norm_upper_rad_s)
        if not all(math.isfinite(v) and v>=0.0 for v in vals):
            raise ValueError("bias certificate bounds must be finite and nonnegative")


def continuation_admitted(cert: BiasContinuationCertificate,
                          limits: BiasLimits) -> bool:
    """Fail closed: neither a noise box nor finite-prefix statistics is FAST."""
    def covered(given,required):
        return (given is not None and required is not None
                and given.horizon_s == required.horizon_s and given.cap <= required.cap)
    return bool(
        cert.accel_slow_norm_upper_mps2 <= limits.B_a_s_mps2
        and cert.accel_slow_rate_norm_upper_mps3 <= limits.D_a_s_mps3
        and cert.gyro_slow_norm_upper_rad_s <= limits.B_g_s_rad_s
        and cert.gyro_slow_rate_norm_upper_rad_s2 <= limits.D_g_s_rad_s2
        and cert.accel_fast_norm_upper_mps2 <= limits.B_a_f_mps2
        and cert.gyro_fast_norm_upper_rad_s <= limits.B_g_f_rad_s
        and covered(cert.accel_fast_window,limits.accel_fast_window)
        and covered(cert.gyro_fast_window,limits.gyro_fast_window)
        and cert.decomposition_id and cert.temporal_evidence_id
        and cert.delivered_sample_scope_certified and cert.all_placed_windows_certified
        and cert.predecessor_continuity_certified
        and cert.calibration_scope_qualified and cert.all_time_continuation_certified
    )


@dataclass(frozen=True)
class BiasSample:
    history_id: str
    t_s: float
    b_a_s_mps2: Vec3
    b_g_s_rad_s: Vec3
    b_a_f_mps2: Vec3 = (0.0,0.0,0.0)
    b_g_f_rad_s: Vec3 = (0.0,0.0,0.0)
    decomposition_id: str = "single-carried-split"

    def __post_init__(self) -> None:
        if not self.history_id or not self.decomposition_id or not math.isfinite(self.t_s):
            raise ValueError("valid history, decomposition and time required")
        for x in (self.b_a_s_mps2,self.b_g_s_rad_s,self.b_a_f_mps2,self.b_g_f_rad_s):
            vec(x)


def successor_allowed(current: Sequence[float],successor: Sequence[float],dt_s: float,
                      amplitude_limit: float,rate_limit: float,*,tolerance: float=0.0) -> bool:
    """SLOW predecessor test only; FAST has a separate carried-window test."""
    if not (math.isfinite(dt_s) and dt_s>0.0): return False
    if not all(math.isfinite(v) and v >= 0.0 for v in (amplitude_limit,rate_limit,tolerance)): return False
    return (norm(current)<=amplitude_limit+tolerance and
            norm(successor)<=amplitude_limit+tolerance and
            norm(sub(successor,current))<=slow_change(amplitude_limit,rate_limit,dt_s)+tolerance)


def audit_bias_trace(samples: Sequence[BiasSample],limits: BiasLimits,*,tolerance: float=0.0) -> dict:
    """Finite-prefix audit. Last timestamp closes the preceding fast hold cell.

    No proof boundary resets a bias or its accumulation. The last fast value
    has no observed hold duration and its temporal extension is not certified.
    Passing this audit never certifies an all-time physical continuation.
    """
    if not samples: raise ValueError("empty bias history")
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("finite nonnegative tolerance required")
    failures=[]; hid=samples[0].history_id; split=samples[0].decomposition_id
    for i,s in enumerate(samples):
        if s.history_id!=hid or s.decomposition_id!=split:
            failures.append(f"sample {i}: detached history/decomposition")
        for name,value,bound in (("accelerometer slow",s.b_a_s_mps2,limits.B_a_s_mps2),
                                 ("gyro slow",s.b_g_s_rad_s,limits.B_g_s_rad_s),
                                 ("accelerometer fast",s.b_a_f_mps2,limits.B_a_f_mps2),
                                 ("gyro fast",s.b_g_f_rad_s,limits.B_g_f_rad_s)):
            if norm(value)>bound+tolerance: failures.append(f"sample {i}: {name} amplitude")
        if i:
            p=samples[i-1]; dt=s.t_s-p.t_s
            if not successor_allowed(p.b_a_s_mps2,s.b_a_s_mps2,dt,limits.B_a_s_mps2,limits.D_a_s_mps3,tolerance=tolerance):
                failures.append(f"sample {i}: accelerometer slow predecessor/rate")
            if not successor_allowed(p.b_g_s_rad_s,s.b_g_s_rad_s,dt,limits.B_g_s_rad_s,limits.D_g_s_rad_s2,tolerance=tolerance):
                failures.append(f"sample {i}: gyro slow predecessor/rate")
    slow_amplitude_pass=not failures
    temporal={}
    if len(samples)>1 and all(b.t_s>a.t_s for a,b in zip(samples,samples[1:])):
        for name,field,bound,window in (
            ("accel","b_a_f_mps2",limits.B_a_f_mps2,limits.accel_fast_window),
            ("gyro","b_g_f_rad_s",limits.B_g_f_rad_s,limits.gyro_fast_window)):
            temporal[name]=audit_fast_prefix([s.t_s for s in samples],
                [getattr(s,field) for s in samples[:-1]],bound,window,tolerance=tolerance)
            if temporal[name]["finite_prefix_pass"] is False:
                failures.append(f"{name}: fast all-placed-window prefix")
    known=len(temporal)==2 and all(v["finite_prefix_pass"] is not None for v in temporal.values())
    return {"finite_prefix_pass":False if failures else (True if known else None),
            "slow_and_amplitude_prefix_pass":slow_amplitude_pass,
            "failures":failures,"history_id":hid,"decomposition_id":split,
            "fast_prefixes":temporal,"temporal_qualification":"CONDITIONAL_PREFIX" if known else "OPEN",
            "all_time_certified":False,"last_fast_cell_extension_certified":False,
            "independent_successor_boxes_allowed":False,"physical_bias_reset_at_estimator_release":False}

def estimator_prediction_factor(mode: str,phi_ou: float) -> float:
    """Literal accelerometer-bias estimate prediction factor for the shipping mode."""
    if not math.isfinite(phi_ou) or not (0.0 <= phi_ou <= 1.0):
        raise ValueError("phi_ou must be finite in [0,1]")
    if mode=="H18": return 1.0
    if mode=="A21": return phi_ou
    raise ValueError("mode must be H18 or A21")

@dataclass(frozen=True)
class BiasPredictionRelation:
    """One physical SLOW bias law with a mode-dependent estimator prediction."""
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
    is the same bounded/rate-bounded SLOW history in both modes. FAST error
    enters the actual sensor operation, not this physical slow-bias recurrence.
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

def project_estimate(x: Sequence[float], radius: float) -> Vec3:
    """Shipping estimate-projection branches, interpreted in real arithmetic.

    A nonpositive radius disables projection; a nonfinite estimate is reset
    only when projection is enabled. Hardware norm overflow and rounding are
    separate arithmetic obligations, not silently absorbed into physical bias.
    """
    if not math.isfinite(radius):
        raise ValueError("finite radius required for the proof domain")
    if len(x) != 3:
        raise ValueError("expected three-vector")
    y = tuple(float(v) for v in x)
    if radius <= 0.0:
        return y
    if not all(math.isfinite(v) for v in y):
        return (0.0, 0.0, 0.0)
    n = math.hypot(*y)
    return y if n <= radius else scale(radius / n, y)


def projection_defect(estimated_bias_corrected: Sequence[float],radius: float) -> Vec3:
    """r_proj=bhat_corr-Proj_R(bhat_corr), matching the shipping estimate projection."""
    y=vec(estimated_bias_corrected)
    return sub(y,project_estimate(y,radius))

def projected_error(b_true: Sequence[float],estimated_bias_corrected: Sequence[float],
                    radius: float) -> Vec3:
    """Post-projection error b_true-Proj_R(bhat_corr)."""
    return sub(b_true,project_estimate(estimated_bias_corrected,radius))


def prediction_joint_blocks(mode: str, phi_ou: float) -> tuple:
    """Coefficients for z=[e_b;b_true], each scalar multiplying I_3.

    z_next^- = (A tensor I_3) z^+ + (B tensor I_3) w.
    B has ONE shared physical-increment column, not independent error/truth
    disturbances. Correction, projection and frame changes are NOT prediction.
    """
    phi_e = estimator_prediction_factor(mode, phi_ou)
    return ((phi_e, 1.0 - phi_e), (0.0, 1.0)), ((1.0,), (1.0,))


@dataclass(frozen=True)
class BiasCorrectionRelation:
    """Separate correction and projection on one unchanged physical bias."""
    b_true: Vec3
    estimate_minus: Vec3
    applied_increment: Vec3
    estimate_corrected: Vec3
    estimate_plus: Vec3
    e_minus: Vec3
    e_corrected: Vec3
    projection_defect: Vec3
    e_plus: Vec3


def correction_projection_relation(b_true: Sequence[float],
                                   estimate_minus: Sequence[float],
                                   applied_increment: Sequence[float],
                                   radius: float) -> BiasCorrectionRelation:
    """Retain the actual Kalman increment, then the estimate projection.

    These are finite real-arithmetic identities, not a full-state contraction
    claim. The increment must be sourced from the shipping update, with zero
    increment on a held/rejected correction. A floating-point trace must carry
    its arithmetic residual separately.
    """
    b, estimate, delta = vec(b_true), vec(estimate_minus), vec(applied_increment)
    corrected = add(estimate, delta)
    plus = project_estimate(corrected, radius)
    return BiasCorrectionRelation(
        b, estimate, delta, corrected, plus, sub(b, estimate),
        sub(b, corrected), sub(corrected, plus), sub(b, plus),
    )
