# ruff: noqa: F401, F811
"""Verified branch-and-bound QCQP for linked shaped supply on one history cell.

The domain is the SAME temporal source polytope.  Bounds use interval quadratic
evaluation on source boxes plus exact feasibility pruning by linear facets.
Nonlinear MARINE/MAGNETIC side constraints are callbacks and must certify each
leaf; UNKNOWN fails closed.
"""
from __future__ import annotations
from dataclasses import dataclass
import heapq,math
import numpy as np

@dataclass(frozen=True)
class Box:
 lo:np.ndarray;hi:np.ndarray
 def split(self,j):
  m=(self.lo[j]+self.hi[j])/2;a=self.hi.copy();a[j]=m;b=self.lo.copy();b[j]=m
  return Box(self.lo.copy(),a),Box(b,self.hi.copy())

def linear_feasible(box,A,b):
 # min A_i u over box <= b_i is necessary/sufficient for box intersection not excluded
 mn=np.where(A>=0,A*box.lo,A*box.hi).sum(axis=1)
 return bool(np.all(mn<=b+1e-15))

def quadratic_interval(M,box):
 """Outward-safe elementary interval bound for u'Mu on a box."""
 M=(np.asarray(M,float)+np.asarray(M,float).T)/2;n=len(box.lo);lo=hi=0.
 for i in range(n):
  for j in range(n):
   c=M[i,j]
   vals=(box.lo[i]*box.lo[j],box.lo[i]*box.hi[j],box.hi[i]*box.lo[j],box.hi[i]*box.hi[j])
   if i==j and box.lo[i]<=0<=box.hi[i]: vals=vals+(0.,)
   a=min(vals);z=max(vals)
   lo+=c*(a if c>=0 else z);hi+=c*(z if c>=0 else a)
 return math.nextafter(lo,-math.inf),math.nextafter(hi,math.inf)

def linear_interval(v,box):
 v=np.asarray(v,float);lo=np.where(v>=0,v*box.lo,v*box.hi).sum();hi=np.where(v>=0,v*box.hi,v*box.lo).sum()
 return math.nextafter(float(lo),-math.inf),math.nextafter(float(hi),math.inf)

def source_supply_upper(B,C,e_box,u_box):
 # upper of -2 e'B u-u'C u on independent e/u BOX coordinates, while both
 # remain from one history cell. This primitive is used only inside a leaf.
 # Stack z=[e,u] and evaluate exact quadratic interval.
 ne=len(e_box.lo);nu=len(u_box.lo);Z=Box(np.r_[e_box.lo,u_box.lo],np.r_[e_box.hi,u_box.hi])
 M=np.block([[np.zeros((ne,ne)),-B],[-B.T,-C]])
 return quadratic_interval(M,Z)[1]

@dataclass
class QCQPCertificate:
 upper:float;leaves:int;verified:bool;reason:str;max_width:float

def verified_supply_bnb(B,C,e_box,u_box,A,b,side_callback,*,tol=1e-5,max_leaves=200000):
 """Rigorous upper bound of linked supply on one shared-history cell.

 side_callback(box)-> True/False/None: True whole box satisfies nonlinear side
 constraints, False disjoint, None unresolved and therefore subdivision needed.
 """
 heap=[];counter=0
 def push(ub,box):
  nonlocal counter;counter+=1;heapq.heappush(heap,(-ub,counter,box))
 if not linear_feasible(u_box,A,b):return QCQPCertificate(-math.inf,0,True,"empty linear domain",0)
 ub=source_supply_upper(B,C,e_box,u_box);push(ub,u_box);best=-math.inf;leaves=0;maxw=0.
 while heap:
  neg,_,box=heapq.heappop(heap);ub=-neg
  if ub<=best:continue
  side=side_callback(box)
  if side is False:continue
  widths=box.hi-box.lo;maxw=max(maxw,float(widths.max(initial=0)))
  if side is True and widths.max(initial=0)<=tol:
   best=max(best,ub);leaves+=1;continue
  if leaves+len(heap)>=max_leaves:return QCQPCertificate(max(best,ub),leaves,False,"leaf limit",maxw)
  j=int(np.argmax(widths))
  if widths[j]<=tol:
   if side is None:return QCQPCertificate(max(best,ub),leaves,False,"nonlinear side unresolved",maxw)
   best=max(best,ub);leaves+=1;continue
  for child in box.split(j):
   if linear_feasible(child,A,b):push(source_supply_upper(B,C,e_box,child),child)
 return QCQPCertificate(best,leaves,True,"complete",maxw)

def shipping_side_callback(witness_builder,constants):
 """Adapt a shared-history box to the nonlinear shipping premise callback."""
 from .marine_magnetic_qcqp import callback
 def f(box):
  w=witness_builder(box)
  if w is None:return None
  return callback(w,constants)
 return f

def interval_matrix_quadratic_upper(Mmid,Mrad,box):
 """Rigorous upper of z'Mz for interval M on box."""
 mid=quadratic_interval(Mmid,box)[1]
 R=np.asarray(Mrad,float);absmax=np.maximum(np.abs(box.lo),np.abs(box.hi))
 err=float(absmax@R@absmax)
 return math.nextafter(mid+err,math.inf)

def interval_source_supply_upper(Bmid,Brad,Cmid,Crad,e_box,u_box):
 ne=len(e_box.lo);nu=len(u_box.lo);Z=Box(np.r_[e_box.lo,u_box.lo],np.r_[e_box.hi,u_box.hi])
 M=np.block([[np.zeros((ne,ne)),-np.asarray(Bmid,float)],[-np.asarray(Bmid,float).T,-np.asarray(Cmid,float)]])
 R=np.block([[np.zeros((ne,ne)),np.asarray(Brad,float)],[np.asarray(Brad,float).T,np.asarray(Crad,float)]])
 return interval_matrix_quadratic_upper(M,R,Z)

def verified_interval_supply_bnb(B,C,e_box,u_box,A,b,side_callback,*,tol=1e-5,max_leaves=200000):
 heap=[];counter=0
 def ub(box):return interval_source_supply_upper(B.mid,B.rad,C.mid,C.rad,e_box,box)
 def push(v,box):
  nonlocal counter;counter+=1;heapq.heappush(heap,(-v,counter,box))
 if not linear_feasible(u_box,A,b):return QCQPCertificate(-math.inf,0,True,"empty linear domain",0)
 push(ub(u_box),u_box);best=-math.inf;leaves=0;maxw=0.
 while heap:
  neg,_,box=heapq.heappop(heap);v=-neg
  if v<=best:continue
  side=side_callback(box)
  if side is False:continue
  widths=box.hi-box.lo;maxw=max(maxw,float(widths.max(initial=0)))
  if side is True and widths.max(initial=0)<=tol:best=max(best,v);leaves+=1;continue
  if leaves+len(heap)>=max_leaves:return QCQPCertificate(max(best,v),leaves,False,"leaf limit",maxw)
  j=int(np.argmax(widths))
  if widths[j]<=tol:
   if side is None:return QCQPCertificate(max(best,v),leaves,False,"nonlinear side unresolved",maxw)
   best=max(best,v);leaves+=1;continue
  for child in box.split(j):
   if linear_feasible(child,A,b):push(ub(child),child)
 return QCQPCertificate(best,leaves,True,"complete",maxw)

def interval_source_supply_upper(Bi,Ci,e_box,u_box):
 """Rigorous supply upper for interval B,C over one leaf box."""
 Bm=np.asarray(Bi.mid,float);Br=np.asarray(Bi.rad,float);Cm=np.asarray(Ci.mid,float);Cr=np.asarray(Ci.rad,float)
 ne,nu=Bm.shape;Z=Box(np.r_[e_box.lo,u_box.lo],np.r_[e_box.hi,u_box.hi])
 M=np.block([[np.zeros((ne,ne)),-Bm],[-Bm.T,-Cm]])
 ub=quadratic_interval(M,Z)[1]
 # coefficient uncertainty: 2 sum Br_ij |e_i u_j| + sum Cr_ij |u_i u_j|
 ea=np.maximum(abs(e_box.lo),abs(e_box.hi));ua=np.maximum(abs(u_box.lo),abs(u_box.hi))
 ub += 2*float(ea@Br@ua)+float(ua@Cr@ua)
 return math.nextafter(ub,math.inf)

def verified_interval_supply_bnb(Bi,Ci,e_box,u_box,A,b,side_callback,*,tol=1e-5,max_leaves=200000):
 heap=[];counter=0
 def push(ub,box):
  nonlocal counter;counter+=1;heapq.heappush(heap,(-ub,counter,box))
 if not linear_feasible(u_box,A,b):return QCQPCertificate(-math.inf,0,True,"empty linear domain",0)
 push(interval_source_supply_upper(Bi,Ci,e_box,u_box),u_box);best=-math.inf;leaves=0;maxw=0.
 while heap:
  neg,_,box=heapq.heappop(heap);ub=-neg
  if ub<=best:continue
  side=side_callback(box)
  if side is False:continue
  widths=box.hi-box.lo;maxw=max(maxw,float(widths.max(initial=0)))
  if side is True and widths.max(initial=0)<=tol:best=max(best,ub);leaves+=1;continue
  if leaves+len(heap)>=max_leaves:return QCQPCertificate(max(best,ub),leaves,False,"leaf limit",maxw)
  j=int(np.argmax(widths))
  if widths[j]<=tol:
   if side is None:return QCQPCertificate(max(best,ub),leaves,False,"nonlinear side unresolved",maxw)
   best=max(best,ub);leaves+=1;continue
  for child in box.split(j):
   if linear_feasible(child,A,b):push(interval_source_supply_upper(Bi,Ci,e_box,child),child)
 return QCQPCertificate(best,leaves,True,"complete",maxw)
