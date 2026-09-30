"""Homogeneous compatibility chart for the reduced two-word Gram.

Avoids normalizing a moving rank-one kernel.  On a rank-one stratum the
generalized eigenproblem is scalar; a two-column homogeneous chart handles
rank transitions without division by a vanishing singular vector.
"""
from __future__ import annotations
from .interval_riccati_21 import IMat,add,matmul,scale,transpose
from .rank_loss_interval_factor import generalized_ratio_lower

def reduced_pair(R:IMat,Q_nonmag:IMat,N_terminal:IMat,kernel_gram:IMat,c:float):
    if c<=0: raise ValueError("positive kernel ceiling")
    Q=add(Q_nonmag,scale(kernel_gram,1.0/c))
    return (matmul(matmul(transpose(R),Q),R),
            matmul(matmul(transpose(R),N_terminal),R))

def scalar_ratio(q:float,n:float):
    if q<=0:
        return {"verified":False,"upper":float("inf"),"reason":"nonpositive denominator"}
    return {"verified":True,"upper":n/q,"denominator":q,"numerator":n}

def transition_ratio_lower_action(R:IMat,Q_nonmag:IMat,N_terminal:IMat,
                                  kernel_gram:IMat,c:float):
    q,n=reduced_pair(R,Q_nonmag,N_terminal,kernel_gram,c)
    return generalized_ratio_lower(q,n)

def kernel_prior_floor_for_unit_direction(c:float):
    if c<=0: raise ValueError("positive c")
    return 1.0/c

def terminal_bound_implies_ratio(c:float,B_terminal:float):
    if c<=0 or B_terminal<0: raise ValueError("bounds")
    return c*B_terminal
