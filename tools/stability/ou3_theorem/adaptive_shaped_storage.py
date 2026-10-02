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

def lin_scaling21(tau: float, sigma: float) -> np.ndarray:
    """Identity outside (v,p,S,aw); shipping normalization inside LIN."""
    T=np.eye(21); T[6:18,6:18]=np.kron(scaling(tau,sigma),np.eye(3))
    return T

def normalized_correction21(K: np.ndarray,H: np.ndarray,tau:float,sigma:float) -> np.ndarray:
    T=lin_scaling21(tau,sigma)
    A=np.eye(21)-np.asarray(K,float)@np.asarray(H,float)
    return T@A@np.linalg.inv(T)

def normalized_S_correction21(K: np.ndarray,H: np.ndarray,tau:float,sigma:float) -> np.ndarray:
    H=np.asarray(H,float)
    if not (np.linalg.norm(H[:,12:15]-np.eye(3),ord=np.inf)<1e-10 and
            np.linalg.norm(np.c_[H[:,:12],H[:,15:]],ord=np.inf)<1e-10):
        raise ValueError("not literal S=0 observation row")
    return normalized_correction21(K,H,tau,sigma)
