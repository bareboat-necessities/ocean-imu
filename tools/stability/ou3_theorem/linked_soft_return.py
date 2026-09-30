"""Same-word linked soft recurrence; avoids independent dbar*Hbar suprema."""
from __future__ import annotations
import math

def linked_cauchy(d:float,h:float,direct:float):
    if d<0 or h<0: raise ValueError("nonnegative")
    return {"linked_product":d*h,"direct_sq":direct*direct,
            "cauchy_holds":direct*direct<=d*h*(1+1e-14)}

def rank_one_soft_return(a:float,g2:float,rho:float,lambda2:float,c_next:float):
    if a<0 or g2<0 or lambda2<=0 or c_next<=0 or abs(rho)>1+1e-15:
        raise ValueError("invalid")
    x=min(1.0,max(0.0,rho*rho));m=1.0/c_next
    prior=a*g2
    info=m*x+lambda2*(1.0-x)
    post=prior/(1.0+prior*info)
    return {"propagated_prior_variance":prior,"effective_information":info,
            "posterior_rank_one_variance":post,
            "next_kernel_variance":post*x,
            "next_quotient_variance":post*(1.0-x)}

def linked_product_status():
    return {"independent_dbar_Hbar_required":False,
            "same_word_Pi_cancels_generally":False,
            "replacement":"one-dimensional propagated soft gain and overlap recurrence"}
