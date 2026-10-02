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

def pow_pos(x:I,p:float)->I:
 if x.lo<=0: raise ArithmeticError("positive power base not separated from zero")
 return I(math.nextafter(x.lo**p,-math.inf),math.nextafter(x.hi**p,math.inf))

def spectral_mse_RS(tau:I,sigma:I,TS:I,c_sigma=0.90,cj=.0538,
                    r_a=.0148*.0148*.005)->I:
 # Literal deployed law: C_J*(2 r_a)^(1/14)*(sigma/c_sigma*tau^4)^(6/7)/sqrt(TS)
 sigma_ab=I(max(1e-6,sigma.lo/c_sigma),max(1e-6,sigma.hi/c_sigma))
 u=sigma_ab*tau*tau*tau*tau
 return cj*((2*r_a)**(1/14))*pow_pos(u,6/7)*recip(sqrti(TS))

@dataclass
class ShippingTunerBox:
 """Applied/staged tuple including deployed SpectralMSE R_S channel."""
 tau:I; sigma:I; RS:I; pending:tuple[I,I,I]|None=None
 def stage(self,dt:I,f:I,var:I,noise:I,tau_coeff:float,sigma_coeff:float,
           alpha:I,alpha_rs:I):
  ft=I(max(.05,f.lo),min(1.2,f.hi));tt=(tau_coeff*.5*recip(ft)).clamp(.02,12.)
  wave=I(max(1e-6,var.lo-noise.sq().hi),max(1e-6,var.hi-noise.sq().lo));ss=(sigma_coeff*sqrti(wave)).clamp(0,4.)
  TS=(tt*.02).clamp(.004,.15);rr=spectral_mse_RS(tt,ss,TS,sigma_coeff).clamp(.15,100.)
  ct=(1-alpha)*self.tau+alpha*tt;cs=(1-alpha)*self.sigma+alpha*ss;cr=(1-alpha_rs)*self.RS+alpha_rs*rr
  self.pending=(ct,cs,cr)
 def commit(self):
  if self.pending:self.tau,self.sigma,self.RS=self.pending;self.pending=None
  TS=(self.tau*.02).clamp(.004,.15)
  return self.tau,self.sigma,self.RS,TS

@dataclass
class CausalAdaptationBox:
 mahony:MahonyBox; band:BandBox; variance:VarianceBox; tuner:ShippingTunerBox
 def step(self,dt:I,gyro,acc,wave_frequency:I,noise_sigma:I,two_kp=.2,two_ki=.02,
          tau_coeff=1.38,sigma_coeff=.90,adapt_periods=.40,rs_mult=1.5):
  # Commit y_{k-1}'s staged candidate before y_k, matching shipping timing.
  applied=self.tuner.commit()
  vertical=self.mahony.step(dt,gyro,acc,two_kp,two_ki)
  banded=self.band.step(vertical,dt,wave_frequency)
  var=self.variance.update(dt,banded,wave_frequency)
  sea=.5*recip(wave_frequency.clamp(.05,1.2))
  adapt_sec=adapt_periods*sea
  alpha=1-expi(-dt*recip(adapt_sec))
  # Default slew-log is zero: RS horizon = mult*tau_target, clamped by shipping helper.
  tau_target=(tau_coeff*.5*recip(wave_frequency.clamp(.05,1.2))).clamp(.02,12.)
  rs_sec=rs_mult*tau_target
  alpha_rs=1-expi(-dt*recip(rs_sec))
  self.tuner.stage(dt,wave_frequency,var,noise_sigma,tau_coeff,sigma_coeff,alpha,alpha_rs)
  return {"vertical":vertical,"band":banded,"variance":var,
          "tau":applied[0],"sigma_aw":applied[1],"R_S":applied[2],"T_S":applied[3],
          "pending":self.tuner.pending}
