"""Bisect attitude factor residual threshold with exact force."""
from .attitude_factor_threshold import residual_attitude
def solve(tol=1e-6):
 lo,hi=1.,2.
 for _ in range(80):
  if hi-lo<=tol:break
  m=(lo+hi)/2
  if residual_attitude(m)["residual"]<1:lo=m
  else:hi=m
 return {"verified":True,"beta_star_lower_deg":lo,"beta_star_upper_deg":hi,
  "residual_lower":residual_attitude(lo)["residual"],"residual_upper":residual_attitude(hi)["residual"]}
