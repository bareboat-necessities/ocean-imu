"""Exact two-line rank-one covariance recurrence."""
from __future__ import annotations
import math

def adjacent_line_update(a:float,rho:float,j:float,c_next:float):
    if a<0 or j<0 or c_next<=0 or abs(rho)>1+1e-15:
        raise ValueError("invalid covariance/information parameters")
    x=min(1.0,max(0.0,rho*rho));s2=1.0-x;m=1.0/c_next
    den=1.0+a*(m*x+j*s2)
    post=a/den
    return {"overlap_sq":x,"denominator":den,"total_rank_one_variance":post,
            "next_kernel_variance":post*x,
            "next_quotient_variance":post*s2,
            "next_kernel_ceiling":c_next,
            "kernel_invariance":post*x<=c_next*(1+1e-14)}

def kernel_carry_derivative(a:float,x:float,j:float,c_next:float):
    m=1.0/c_next
    den=1+a*(j+(m-j)*x)
    return a*(1+a*j)/(den*den)

def total_variance_derivative(a:float,x:float,j:float,c_next:float):
    m=1.0/c_next
    den=1+a*(j+(m-j)*x)
    return -a*a*(m-j)/(den*den)
