"""FAST contribution to guard LP-minus-raw state from signed primitive.

For d_k=(L-I)u at a fixed k, collect the causal FIR/IIR coefficients q_j over
the qualified history. With U_j=h sum_{i<=j}u_i, Abel summation gives a bound
(C/h)*variation(q), plus the pre-horizon IIR tail charged by B. This is applied
ONCE to d_k, not recursively as 2B at every sample.
"""
from __future__ import annotations
import math

def one_pole_q(alpha,n):
 # lags 0..n-1 for L-I
 return [-alpha]+[(1-alpha)*(alpha**j) for j in range(1,n)]

def cascade_impulse(alpha,poles,n):
 # repeated convolution of positive one-pole impulse, truncated to n
 h=[1.]+[0.]*(n-1)
 base=[(1-alpha)*(alpha**j) for j in range(n)]
 for _ in range(poles):
  z=[0.]*n
  for i in range(n):
   for j in range(n-i):z[i+j]+=h[i]*base[j]
  h=z
 return h

def lp_minus_i_q(alpha,poles,n):
 h=cascade_impulse(alpha,poles,n);h[0]-=1.;return h

def abel_charge(q):
 if not q:return 0.
 # U before the first retained sample is a boundary primitive; both endpoints.
 return abs(q[0])+abs(q[-1])+sum(abs(q[j+1]-q[j]) for j in range(len(q)-1))

def fast_lp_minus_raw_bound(alpha,poles,dt,B,C,H,elapsed):
 if elapsed<=0:return 0.
 n=min(max(1,int(math.ceil(elapsed/dt))+1),max(1,int(math.floor(H/dt))))
 q=lp_minus_i_q(alpha,poles,n)
 inside=(C/dt)*abel_charge(q)
 # State inherited from input older than H: positive LP cascade tail <= B times
 # remaining unit mass. Raw -I term has no old-input tail.
 h=cascade_impulse(alpha,poles,n)
 tail_mass=max(0.,1.-sum(h))
 tail=B*tail_mass
 # Pointwise 2B is always valid; select the tighter theorem consequence.
 return min(2*B,inside+tail)

def audit(alpha,poles,dt,B,C,H,elapsed):
 z=fast_lp_minus_raw_bound(alpha,poles,dt,B,C,H,elapsed)
 return {"bound":z,"pointwise":2*B,"uses_signed_primitive":True,
         "paid_once_per_lp_state":True,"elapsed":elapsed}
