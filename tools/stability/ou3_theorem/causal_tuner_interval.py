"""Causal interval propagation for the measurement-only OU-III tuner path.

Implements rigorous scalar interval recurrences for Mahony algebra (given a
verified quaternion/input cell), adaptive wave band, variance moments, operating
point and one-sample staged tuner commit. Nonlinear normalization fails closed
when its denominator interval reaches zero. Generated quantities are outputs.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class I:
 lo:float; hi:float
 def __post_init__(self):
  if not(math.isfinite(self.lo) and math.isfinite(self.hi) and self.lo<=self.hi): raise ValueError("bad interval")
 def __add__(self,o): o=iv(o);return I(math.nextafter(self.lo+o.lo,-math.inf),math.nextafter(self.hi+o.hi,math.inf))
 __radd__=__add__
 def __neg__(self): return I(-self.hi,-self.lo)
 def __sub__(self,o): return self+-iv(o)
 def __rsub__(self,o): return iv(o)+-self
 def __mul__(self,o):
  o=iv(o);v=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi);return I(math.nextafter(min(v),-math.inf),math.nextafter(max(v),math.inf))
 __rmul__=__mul__
 def sq(self):
  if self.lo<=0<=self.hi:return I(0,math.nextafter(max(self.lo*self.lo,self.hi*self.hi),math.inf))
  return self*self
 def clamp(self,a,b): return I(max(a,self.lo),min(b,self.hi))
def iv(x): return x if isinstance(x,I) else I(float(x),float(x))
def expi(x): return I(math.nextafter(math.exp(x.lo),-math.inf),math.nextafter(math.exp(x.hi),math.inf))
def sqrti(x):
 if x.lo<=0: raise ArithmeticError("normalization denominator not separated from zero")
 return I(math.nextafter(math.sqrt(x.lo),-math.inf),math.nextafter(math.sqrt(x.hi),math.inf))
def recip(x):
 if x.lo<=0<=x.hi: raise ArithmeticError("reciprocal singular")
 return I(math.nextafter(1/x.hi,-math.inf),math.nextafter(1/x.lo,math.inf))
def hull(a,b):return I(min(a.lo,b.lo),max(a.hi,b.hi))

@dataclass
class MahonyBox:
 q:tuple[I,I,I,I]; integ:tuple[I,I,I]; initialized:bool
 def step(self,dt:I,gyro:tuple[I,I,I],acc:tuple[I,I,I],two_kp:float,two_ki:float):
  # Shipping observer is fed -acc; normalization uses exact sqrt enclosure
  ax,ay,az=(-acc[0],-acc[1],-acc[2]); n2=ax.sq()+ay.sq()+az.sq(); inv=recip(sqrti(n2));ax,ay,az=ax*inv,ay*inv,az*inv
  q0,q1,q2,q3=self.q
  hvx=q1*q3-q0*q2; hvy=q0*q1+q2*q3; hvz=.5*(q0.sq()-q1.sq()-q2.sq()+q3.sq())
  ex=ay*hvz-az*hvy;ey=az*hvx-ax*hvz;ez=ax*hvy-ay*hvx
  ints=list(self.integ)
  if two_ki>0:
   ints=[v+two_ki*e*dt for v,e in zip(ints,(ex,ey,ez))]
  else: ints=[I(0,0)]*3
  g=[gyro[i]+ints[i]+two_kp*e for i,e in enumerate((ex,ey,ez))]
  h=.5*dt; gx,gy,gz=(x*h for x in g)
  n0=q0-q1*gx-q2*gy-q3*gz;n1=q1+q0*gx+q2*gz-q3*gy;n2q=q2+q0*gy-q1*gz+q3*gx;n3=q3+q0*gz+q1*gy-q2*gx
  nn=n0.sq()+n1.sq()+n2q.sq()+n3.sq(); rinv=recip(sqrti(nn))
  self.q=(n0*rinv,n1*rinv,n2q*rinv,n3*rinv);self.integ=tuple(ints);self.initialized=True
  # down row dot original specific force
  q0,q1,q2,q3=self.q
  down=(2*(q1*q3-q0*q2),2*(q2*q3+q0*q1),q0.sq()-q1.sq()-q2.sq()+q3.sq())
  return -(down[0]*acc[0]+down[1]*acc[1]+down[2]*acc[2])

@dataclass
class BandBox:
 low:I; band:I; p00:I; p01:I; p11:I
 def step(self,x:I,dt:I,f:I,low_ratio=.5,high_ratio=4.,min_hz=.01,max_hz=6.):
  # Require one clamp stratum; caller splits shared history if bounds cross it.
  ny=I(.45,.45)*recip(dt); upper=I(min(max_hz,ny.lo),min(max_hz,ny.hi))
  rawlo=f*low_ratio; rawhi=f*high_ratio
  if rawlo.lo<min_hz<rawlo.hi or rawhi.lo<upper.lo<rawhi.hi: raise ArithmeticError("band clamp stratum unresolved")
  lo=I(max(min_hz,rawlo.lo),max(min_hz,rawlo.hi)); hi=I(min(upper.lo,rawhi.lo),min(upper.hi,rawhi.hi))
  if hi.lo<=lo.hi: raise ArithmeticError("band ordering unresolved")
  al=1-expi(-2*math.pi*lo*dt);ah=1-expi(-2*math.pi*hi*dt);ql=1-al;qh=1-ah
  lp=ql*self.low+al*x;band=qh*self.band+ah*ql*(x-self.low)
  a00=ql;a10=-ah*ql;a11=qh;b0=al;b1=ah*ql
  p00=a00*a00*self.p00+b0*b0;p01=a00*(a10*self.p00+a11*self.p01)+b0*b1;p11=a10*a10*self.p00+2*a10*a11*self.p01+a11*a11*self.p11+b1*b1
  self.low,self.band,self.p00,self.p01,self.p11=lp,band,p00,p01,p11;return band

@dataclass
class VarianceBox:
 mean:I; mean_w:I; sq:I; sq_w:I
 def update(self,dt:I,a:I,f:I,K=4.,tmin=.3,tmax=60.):
  fe=I(max(.05,f.lo),min(5.,f.hi)); sea=.5*recip(fe); tau=I(max(tmin,K*2*sea.lo),min(tmax,K*2*sea.hi))
  alpha=1-expi(-dt*recip(tau))
  self.mean=(1-alpha)*self.mean+alpha*a;self.mean_w=(1-alpha)*self.mean_w+alpha
  self.sq=(1-alpha)*self.sq+alpha*a.sq();self.sq_w=(1-alpha)*self.sq_w+alpha
  if self.mean_w.lo<=1e-12 or self.sq_w.lo<=1e-12:return I(0,max(0,self.sq.hi))
  mu=self.mean*recip(self.mean_w);m2=self.sq*recip(self.sq_w);return I(max(0,m2.lo-mu.sq().hi),max(0,m2.hi-mu.sq().lo))

@dataclass
class TunerBox:
 tau:I; sigma:I; pending:tuple[I,I]|None=None
 def stage(self,dt:I,f:I,var:I,noise:I,tau_coeff:float,sigma_coeff:float,alpha:I):
  ft=I(max(.05,f.lo),min(1.2,f.hi)); target_tau=(tau_coeff*.5*recip(ft)).clamp(.02,12.)
  wave=I(max(1e-6,var.lo-noise.sq().hi),max(1e-6,var.hi-noise.sq().lo));target_sigma=(sigma_coeff*sqrti(wave)).clamp(0,4.)
  cand_tau=(1-alpha)*self.tau+alpha*target_tau;cand_sigma=(1-alpha)*self.sigma+alpha*target_sigma
  self.pending=(cand_tau,cand_sigma)
 def commit(self):
  if self.pending:self.tau,self.sigma=self.pending;self.pending=None
  TS=(self.tau*.02).clamp(.004,.15)
  return self.tau,self.sigma,TS
