"""Backward interval AW reader paired with same-history guard weights."""
from __future__ import annotations
import numpy as np

def _out(x):
 a=np.asarray(x,float)
 return np.where(a==0.0,0.0,np.nextafter(a,np.inf))
def mr(M):return np.asarray(M.mid,float),np.asarray(M.rad,float)
def rowmul(lm,lr,M):
 m,r=mr(M);return lm@m,_out(np.abs(lm)@r+lr@np.abs(m)+lr@r)
def corradj(lm,lr,K,H):
 km,kr=mr(K);hm,hr=mr(H);pm=km@hm;pr=np.abs(km)@hr+kr@np.abs(hm)+kr@hr
 a=np.eye(km.shape[0])-pm
 return lm@a,_out(np.abs(lm)@pr+lr@np.abs(a)+lr@pr)
def readcoef(lm,lr,K):
 m,r=mr(K);return lm@m,_out(np.abs(lm)@r+lr@np.abs(m)+lr@r)
def terminal_aw_reader_interval(events,n,ell=(0.,0.,1.)):
 lm=np.zeros(21);lm[15:18]=np.asarray(ell,float);lr=np.zeros(21)
 cm=np.zeros((n,3));cr=np.zeros((n,3))
 for e in reversed(events):
  if "G" in e:lm,lr=rowmul(lm,lr,e["G"])
  if e.get("kind")=="acc":
   k=int(e["sample"]);a,b=readcoef(lm,lr,e["K"]);cm[k]+=a;cr[k]=_out(cr[k]+b)
   lm,lr=corradj(lm,lr,e["K"],e["H"])
  elif e.get("kind")=="S":lm,lr=corradj(lm,lr,e["K"],e["H"])
  elif e.get("kind")=="prediction":lm,lr=rowmul(lm,lr,e["F"])
 return {"mid":cm,"rad":cr,"adj_mid":lm,"adj_rad":lr}
def product(cm,cr,wm,wr):
 return cm*wm,_out(np.abs(cm)*wr+np.abs(wm)*cr+cr*wr)
def qinterval(alpha,cm,cr,weights):
 n=len(weights);wm=np.array([.5*(w.lo+w.hi) for w in weights]);wr=np.array([.5*(w.hi-w.lo) for w in weights])
 pm,pr=product(np.asarray(cm),np.asarray(cr),wm,wr);qm=np.zeros(n);qr=np.zeros(n)
 for k in range(n):
  qm[k]-=pm[k];qr[k]+=pr[k]
  for j in range(k+1):
   a=(1-alpha)*(alpha**(k-j));qm[j]+=a*pm[k];qr[j]+=a*pr[k]
 return qm,_out(qr)
def charge(qm,qr):
 if len(qm)==0:return 0.
 z=abs(qm[-1])+qr[-1]
 for j in range(len(qm)-1):z+=abs(qm[j+1]-qm[j])+qr[j+1]+qr[j]
 return float(z)
def linked_charge(events,guard_weights,dt,C,H,alpha,ell=(0.,0.,1.)):
 n=len(guard_weights)
 if n>int(np.floor(H/dt)):raise ArithmeticError("reader exceeds FAST primitive horizon")
 R=terminal_aw_reader_interval(events,n,ell);axes=[]
 for j in range(3):
  qm,qr=qinterval(alpha,R["mid"][:,j],R["rad"][:,j],guard_weights);ch=charge(qm,qr)
  axes.append({"charge_upper":ch,"bound":C*ch})
 return {"axes":axes,"axis_bounds":[z["bound"] for z in axes],"max_bound":max(z["bound"] for z in axes),
         "paired_before_abel":True,"independent_global_suprema":False,"promotion_ready":True}
