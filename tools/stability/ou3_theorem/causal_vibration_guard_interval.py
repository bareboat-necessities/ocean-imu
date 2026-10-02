# ruff: noqa: F401, F811
"""Interval replica of the literal AccelVibrationGuard state.

Uses scalar interval arithmetic componentwise. Branch crossings at target clamp
or weight rail are represented by hulls, never midpoint choices.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from .causal_tuner_interval import I,expi,sqrti

def hull(a,b):return I(min(a.lo,b.lo),max(a.hi,b.hi))
def clamp01(x):
 if x.hi<=0:return I(0.,0.)
 if x.lo>=1:return I(1.,1.)
 return I(max(0.,x.lo),min(1.,x.hi))
def square_outer(x):
 if x.lo<=0<=x.hi:return I(0.,max(x.lo*x.lo,x.hi*x.hi))
 return I(min(x.lo*x.lo,x.hi*x.hi),max(x.lo*x.lo,x.hi*x.hi))

@dataclass
class GuardBox:
 stages:list
 detect:list
 removed_ms:list
 weight:I
 lp_raw_diff_norm:list
 smooth_lp_raw_diff_norm:list
 elapsed:float=0.0
 initialized:bool=False
 cutoff_hz:float=3.0
 poles:int=2
 detect_hz:float=25.
 engage_lo:float=.03
 engage_hi:float=.08
 slew_tau:float=5.
 def step(self,acc,dt,raw_norm_lower=None,raw_increment_norm_upper=None,fast_state_bound=None):
  acc=tuple(acc)
  if not self.initialized:
   self.stages=[list(acc) for _ in range(4)]
   self.detect=[list(acc),[I(0,0) for _ in range(3)]]
   self.removed_ms=[I(0,0) for _ in range(3)];self.weight=I(0,0);self.lp_raw_diff_norm=[0.0]*4;self.smooth_lp_raw_diff_norm=[0.0]*4;self.elapsed=0.0
   self.initialized=True
   return {"conditioned":acc,"rms":I(0,0),"excess":I(0,0),"weight":self.weight,
           "conditioning_delta_norm_upper":0.0,
           "conditioned_norm_lower":raw_norm_lower}
  alpha=expi(I(-2*math.pi*self.cutoff_hz*dt,-2*math.pi*self.cutoff_hz*dt))
  # Coupled recurrence for d_p = stage_p - raw_k.  Let
  # delta_k = ||raw_k-raw_{k-1}||. For p=0,
  # d_0,k = alpha(d_0,k-1 - (raw_k-raw_{k-1})).
  # For p>0, stage_p,k=(1-alpha)stage_{p-1,k}+alpha stage_p,k-1,
  # so ||d_p,k|| <= (1-alpha)||d_{p-1,k}||
  #                  + alpha(||d_p,k-1||+delta_k).
  # raw_increment_norm_upper is supplied from the same physical/sensor history.
  delta_raw=raw_increment_norm_upper
  low=list(acc)
  old_diff=list(self.lp_raw_diff_norm)
  old_smooth=list(self.smooth_lp_raw_diff_norm)
  new_diff=list(old_diff);new_smooth=list(old_smooth)
  for p in range(self.poles):
   nxt=[]
   for a in range(3):
    self.stages[p][a]=(1-alpha)*low[a]+alpha*self.stages[p][a]
    nxt.append(self.stages[p][a])
   low=nxt
   if delta_raw is not None:
    ah=alpha.hi
    if p==0:new_smooth[p]=ah*(old_smooth[p]+delta_raw)
    else:new_smooth[p]=(1-alpha.lo)*new_smooth[p-1]+ah*(old_smooth[p]+delta_raw)
  self.smooth_lp_raw_diff_norm=new_smooth
  self.elapsed+=dt
  # FAST is not recursively charged in delta_raw. Its complete LP-I state is
  # supplied by one signed-primitive Abel bound at the current elapsed time.
  fb=0.0 if fast_state_bound is None else float(fast_state_bound)
  new_diff=[x+fb for x in new_smooth]
  self.lp_raw_diff_norm=new_diff
  gamma=expi(I(-2*math.pi*self.detect_hz*dt,-2*math.pi*self.detect_hz*dt))
  hp=list(acc)
  for p in range(2):
   nxt=[]
   for a in range(3):
    self.detect[p][a]=(1-gamma)*hp[a]+gamma*self.detect[p][a]
    nxt.append(hp[a]-self.detect[p][a])
   hp=nxt
  beta=1-math.exp(-2*math.pi*.05*dt)
  for a in range(3):
   self.removed_ms[a]=(1-beta)*self.removed_ms[a]+beta*square_outer(hp[a])
  total=self.removed_ms[0]+self.removed_ms[1]+self.removed_ms[2]
  rms=I(math.sqrt(max(0,total.lo)),math.sqrt(max(0,total.hi)))
  excess=I(max(0,rms.lo-self.engage_lo),max(0,rms.hi-self.engage_lo))
  if self.engage_hi>self.engage_lo:
   target=clamp01((rms-self.engage_lo)*(1/(self.engage_hi-self.engage_lo)))
  else:
   if rms.hi<self.engage_lo:target=I(0,0)
   elif rms.lo>=self.engage_lo:target=I(1,1)
   else:target=I(0,1)
  slew=1-math.exp(-dt/self.slew_tau)
  w=self.weight+slew*(target-self.weight)
  # Rail parking is discontinuous only inside epsilon neighborhoods; hull both.
  eps=1e-4
  if w.hi<eps:w=I(0,0)
  elif w.lo>1-eps:w=I(1,1)
  else:
   lo=0. if w.lo<eps else w.lo;hi=1. if w.hi>1-eps else w.hi;w=I(lo,hi)
  self.weight=w
  conditioned=tuple(acc[a]+w*(low[a]-acc[a]) for a in range(3))
  delta=tuple(conditioned[a]-acc[a] for a in range(3))
  dn=math.sqrt(sum(max(abs(x.lo),abs(x.hi))**2 for x in delta))
  # Exact output relation is conditioned=raw+w*(LP-raw).  Use the coupled
  # LP/raw difference recurrence, not the component-box delta.  Since w>=0,
  # ||conditioned|| >= ||raw||-w_hi*||LP-raw||.
  coupled=None
  lp_diff=new_diff[self.poles-1] if self.poles else 0.0
  if raw_norm_lower is not None:
   coupled=max(0.0,float(raw_norm_lower)-w.hi*lp_diff)
  return {"conditioned":conditioned,"rms":rms,"excess":excess,"weight":w,
          "conditioning_delta_norm_upper":dn,"lp_raw_diff_norm_upper":lp_diff,"conditioned_norm_lower":coupled,
          "branch_midpoint_used":False}

def initial_guard():
 z=I(0,0)
 return GuardBox([[z,z,z] for _ in range(4)],[[z,z,z] for _ in range(2)],[z,z,z],z,[0.0]*4,[0.0]*4,0.0)
