"""Shipping-faithful shaped-storage algebra for OU-III.

No independent tuner boxes.  This module only defines exact coordinate maps
used on one generated shipping history.
"""
from __future__ import annotations
import math
import numpy as np

def scaling(tau: float, sigma: float) -> np.ndarray:
    if not (math.isfinite(tau) and tau > 0 and math.isfinite(sigma) and sigma > 0):
        raise ValueError("reachable positive applied tau/sigma required")
    return np.diag([1/(sigma*tau),1/(sigma*tau*tau),
                    1/(sigma*tau**3),1/sigma])

def normalized_prediction_exact(h: float, tau: float) -> np.ndarray:
    x=h/tau
    a=math.exp(-x)
    em1=math.expm1(-x)
    return np.array([
        [1.0,0.0,0.0,1.0-a],
        [x,1.0,0.0,x+em1],
        [0.5*x*x,x,1.0,0.5*x*x-x-em1],
        [0.0,0.0,0.0,a]],dtype=float)

def physical_prediction(h: float, tau: float) -> np.ndarray:
    x=h/tau; a=math.exp(-x); em1=math.expm1(-x)
    return np.array([
        [1.0,0.0,0.0,-tau*em1],
        [h,1.0,0.0,tau*tau*(x+em1)],
        [0.5*h*h,h,1.0,tau**3*(0.5*x*x-x-em1)],
        [0.0,0.0,0.0,a]],dtype=float)

def normalized_factor(A: np.ndarray, tau_before: float, sigma_before: float,
                      tau_after: float, sigma_after: float) -> np.ndarray:
    return scaling(tau_after,sigma_after) @ A @ np.linalg.inv(scaling(tau_before,sigma_before))

def tuner_coboundary(tau_before: float, sigma_before: float,
                      tau_after: float, sigma_after: float) -> np.ndarray:
    return scaling(tau_after,sigma_after) @ np.linalg.inv(scaling(tau_before,sigma_before))

def prediction_parity(h: float,tau: float,sigma:float) -> float:
    got=normalized_factor(physical_prediction(h,tau),tau,sigma,tau,sigma)
    return float(np.linalg.norm(got-normalized_prediction_exact(h,tau),ord=np.inf))
