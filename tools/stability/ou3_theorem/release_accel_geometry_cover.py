"""Adaptive one-step release innovation cover in shared attitude/force coordinates."""
from __future__ import annotations
from dataclasses import dataclass
import heapq,math,numpy as np
from .golive_release_seed import covariance_interval
from .release_interval_propagation import one_sample_boxes
from .interval_riccati_21 import predict_covariance,verified_joseph_update_psd,diagonal_interval
from .accel_geometry_cell import accel_h_from_force_rotation,rotation_cell_from_ball

@dataclass(frozen=True)
class AFCell:
 theta_lo:np.ndarray;theta_hi:np.ndarray;f_lo:np.ndarray;f_hi:np.ndarray;depth:int=0
 def widths(self):return np.r_[self.theta_hi-self.theta_lo,self.f_hi-self.f_lo]
 def split(self,j):
  if j<3:lo=self.theta_lo;hi=self.theta_hi;fl=self.f_lo;fh=self.f_hi
  else:lo=self.f_lo;hi=self.f_hi;fl=self.theta_lo;fh=self.theta_hi
  if j<3:
   m=(self.theta_lo[j]+self.theta_hi[j])/2
   a,b=self.theta_hi.copy(),self.theta_lo.copy();a[j]=m;b[j]=m
   return AFCell(self.theta_lo.copy(),a,self.f_lo.copy(),self.f_hi.copy(),self.depth+1),AFCell(b,self.theta_hi.copy(),self.f_lo.copy(),self.f_hi.copy(),self.depth+1)
  k=j-3;m=(self.f_lo[k]+self.f_hi[k])/2
  a,b=self.f_hi.copy(),self.f_lo.copy();a[k]=m;b[k]=m
  return AFCell(self.theta_lo.copy(),self.theta_hi.copy(),self.f_lo.copy(),a,self.depth+1),AFCell(self.theta_lo.copy(),self.theta_hi.copy(),b,self.f_hi.copy(),self.depth+1)

def root_cell():
 # Corrected captured domain <=6.9 deg. Force vector = gravity + physical a,
 # so each component lies in +/- (g+Amax); norm coupling is retained elsewhere.
 th=math.radians(6.9);fm=9.80665+8.8
 return AFCell(np.full(3,-th),np.full(3,th),np.full(3,-fm),np.full(3,fm))

def attempt(cell):
 F,Q,_,Racc,*_=one_sample_boxes();P=predict_covariance(covariance_interval(),F,Q)
 fm=(cell.f_lo+cell.f_hi)/2;fr=(cell.f_hi-cell.f_lo)/2
 try:
  Rm,Rr=rotation_cell_from_ball(math.radians(6.9))
  H=accel_h_from_force_rotation(fm,fr,Rm,Rr)
  P1,c=verified_joseph_update_psd(P,H,Racc)
  return {"verified":True,"P":P1,"innovation":c}
 except Exception as e:return {"verified":False,"reason":str(e)}

def adaptive_cover(max_depth=10,max_cells=4096):
 root=root_cell();heap=[(-float(np.max(root.widths())),0,root)];leaves=[];bad=[];counter=0
 while heap:
  _,_,c=heapq.heappop(heap);z=attempt(c)
  if z["verified"]:leaves.append((c,z));continue
  if c.depth>=max_depth or len(heap)+len(leaves)+len(bad)>=max_cells:
   bad.append((c,z));continue
  # Normalize theta by captured radius and force by force radius.
  # Exclude force boxes wholly outside the physical sphere before splitting.
  fnear=np.where((c.f_lo<=0)&(c.f_hi>=0),0,np.minimum(abs(c.f_lo),abs(c.f_hi)))
  if np.linalg.norm(fnear)>18.60665: continue
  w=(c.f_hi-c.f_lo)/18.60665
  j=3+int(np.argmax(w));a,b=c.split(j)
  for x in (a,b):
   counter+=1;heapq.heappush(heap,(-float(np.max(x.widths())),counter,x))
 return {"verified":not bad,"leaf_count":len(leaves),"unresolved_count":len(bad),
         "max_depth":max((x.depth for x,_ in leaves),default=0),
         "first_unresolved":bad[0][1] if bad else None}
