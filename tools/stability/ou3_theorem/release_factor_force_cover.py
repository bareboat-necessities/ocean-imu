"""Adaptive force-sector cover for PSD-factor innovation certificate."""
from __future__ import annotations
import math,heapq,numpy as np
from dataclasses import dataclass
from .golive_release_seed import covariance_factor_interval
from .release_interval_propagation import one_sample_boxes
from .accel_geometry_cell import accel_h_from_force_rotation,rotation_cell_from_ball
from .psd_factor_innovation import factor_innovation_probe

RHO=18.60665
@dataclass(frozen=True)
class FCell:
 lo:np.ndarray;hi:np.ndarray;depth:int=0
 def split(self,j):
  m=(self.lo[j]+self.hi[j])/2;a=self.hi.copy();b=self.lo.copy();a[j]=m;b[j]=m
  return FCell(self.lo.copy(),a,self.depth+1),FCell(b,self.hi.copy(),self.depth+1)
 def disjoint_ball(self):
  near=np.where((self.lo<=0)&(self.hi>=0),0,np.minimum(abs(self.lo),abs(self.hi)))
  return np.linalg.norm(near)>RHO
 def contained_ball(self):
  far=np.maximum(abs(self.lo),abs(self.hi));return np.linalg.norm(far)<=RHO

def attempt(c):
 _,_,_,Racc,*_=one_sample_boxes();Rm,Rr=rotation_cell_from_ball(math.radians(6.9))
 fm=(c.lo+c.hi)/2;fr=(c.hi-c.lo)/2
 H=accel_h_from_force_rotation(fm,fr,Rm,Rr)
 return factor_innovation_probe(H,covariance_factor_interval(),Racc)

def cover(max_depth=18,max_cells=200000):
 root=FCell(np.full(3,-RHO),np.full(3,RHO));heap=[(-2*RHO,0,root)];n=0;good=[];bad=[]
 while heap:
  _,_,c=heapq.heappop(heap)
  if c.disjoint_ball():continue
  z=attempt(c)
  if z["verified"]:good.append((c,z));continue
  if c.depth>=max_depth or len(heap)+len(good)+len(bad)>=max_cells:bad.append((c,z));continue
  j=int(np.argmax(c.hi-c.lo));a,b=c.split(j)
  for x in (a,b):n+=1;heapq.heappush(heap,(-float(np.max(x.hi-x.lo)),n,x))
 return {"verified":not bad,"leaf_count":len(good),"unresolved_count":len(bad),
  "max_depth":max((c.depth for c,_ in good),default=0),
  "worst_verified_residual":max((z["residual"] for _,z in good),default=None),
  "first_unresolved":bad[0][1] if bad else None}
