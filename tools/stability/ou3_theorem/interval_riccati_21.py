"""Midpoint-radius matrix kernels for the literal 21-state A21 covariance map.

These kernels preserve every covariance entry.  They are deliberately generic:
shipping-specific F/Q/H/R interval constructors feed them, while innovation
inverses are accepted only with a separate verified residual certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

from tools.stability.ou3_theorem.interval_riccati import (
    symmetric_interval_gershgorin, verified_interval_innovation_inverse,
)

N=21

@dataclass(frozen=True)
class IMat:
    mid: tuple[tuple[float,...],...]
    rad: tuple[tuple[float,...],...]

    def __post_init__(self) -> None:
        n=len(self.mid)
        if n==0 or len(self.rad)!=n:
            raise ValueError("nonempty matrix required")
        m=len(self.mid[0])
        if m==0 or any(len(r)!=m for r in self.mid) or any(len(r)!=m for r in self.rad):
            raise ValueError("rectangular midpoint/radius matrices required")
        for a,b in zip(self.mid,self.rad):
            for x,r in zip(a,b):
                if not math.isfinite(x) or not math.isfinite(r) or r<0:
                    raise ValueError("finite midpoint and nonnegative radius required")

    @property
    def shape(self) -> tuple[int,int]:
        return len(self.mid),len(self.mid[0])


def _out(x: float) -> float:
    return math.nextafter(x,math.inf)


def add(a: IMat,b: IMat) -> IMat:
    if a.shape!=b.shape: raise ValueError("shape mismatch")
    m=[];r=[]
    for ar,br,aa,bb in zip(a.mid,b.mid,a.rad,b.rad):
        m.append(tuple(x+y for x,y in zip(ar,br)))
        r.append(tuple(_out(x+y) for x,y in zip(aa,bb)))
    return IMat(tuple(m),tuple(r))


def transpose(a: IMat) -> IMat:
    n,m=a.shape
    return IMat(tuple(tuple(a.mid[i][j] for i in range(n)) for j in range(m)),
                tuple(tuple(a.rad[i][j] for i in range(n)) for j in range(m)))


def matmul(a: IMat,b: IMat) -> IMat:
    n,k=a.shape; k2,m=b.shape
    if k!=k2: raise ValueError("shape mismatch")
    cm=[];cr=[]
    for i in range(n):
        mr=[];rr=[]
        for j in range(m):
            mid=0.0;rad=0.0
            for t in range(k):
                am,ar=a.mid[i][t],a.rad[i][t]
                bm,br=b.mid[t][j],b.rad[t][j]
                mid += am*bm
                rad += abs(am)*br+abs(bm)*ar+ar*br
            mr.append(mid);rr.append(_out(rad))
        cm.append(tuple(mr));cr.append(tuple(rr))
    return IMat(tuple(cm),tuple(cr))


def scale(a: IMat,s: float) -> IMat:
    if not math.isfinite(s): raise ValueError("finite scale required")
    return IMat(tuple(tuple(s*x for x in row) for row in a.mid),
                tuple(tuple(abs(s)*x for x in row) for row in a.rad))


def symmetrize(a: IMat) -> IMat:
    if a.shape[0]!=a.shape[1]: raise ValueError("square matrix required")
    return scale(add(a,transpose(a)),.5)


def predict_covariance(p: IMat,f: IMat,q: IMat) -> IMat:
    """Literal P-=F P F'+Q with all cross-covariances retained."""
    if p.shape!=(N,N) or f.shape!=(N,N) or q.shape!=(N,N):
        raise ValueError("A21 prediction is 21-state")
    return symmetrize(add(matmul(matmul(f,p),transpose(f)),q))


def innovation_covariance(p: IMat,h: IMat,r: IMat) -> IMat:
    """S=H P H'+R for the shipping 3-D acc/S/mag updates."""
    if p.shape!=(N,N) or h.shape!=(3,N) or r.shape!=(3,3):
        raise ValueError("shipping update dimensions required")
    return symmetrize(add(matmul(matmul(h,p),transpose(h)),r))


def innovation_inverse_spectral_certificate(s: IMat) -> dict:
    """Verify an innovation inverse using interval Gershgorin residual radius."""
    if s.shape!=(3,3): raise ValueError("3x3 innovation required")
    lo,hi=symmetric_interval_gershgorin(s.mid,s.rad)
    # ||E||_2 <= max Gershgorin row sum of the radius matrix.
    r=max(sum(row) for row in s.rad)
    a=min(sum(1 for _ in [0])*s.mid[i][i]-sum(abs(s.mid[i][j]) for j in range(3) if j!=i)
          for i in range(3))
    out=verified_interval_innovation_inverse(
        midpoint_min_eigenvalue=a,spectral_radius_bound=r)
    return {"interval_eigen_lower":lo,"interval_eigen_upper":hi,**out}


def joseph_covariance(p: IMat,k: IMat,h: IMat,r: IMat) -> IMat:
    """Literal Joseph form (I-KH)P(I-KH)'+K R K'."""
    if p.shape!=(N,N) or k.shape!=(N,3) or h.shape!=(3,N) or r.shape!=(3,3):
        raise ValueError("shipping Joseph dimensions required")
    eye=IMat(tuple(tuple(1.0 if i==j else 0.0 for j in range(N)) for i in range(N)),
             tuple(tuple(0.0 for _ in range(N)) for _ in range(N)))
    kh=matmul(k,h)
    a=add(eye,scale(kh,-1.0))
    return symmetrize(add(matmul(matmul(a,p),transpose(a)),
                          matmul(matmul(k,r),transpose(k))))


def aw_covariance_floor_event(p: IMat, target_aw_variance_upper: float,
                              off_aw: int=15) -> IMat:
    """Interval enclosure of shipping PSD a_w floor: only AW block can increase."""
    if p.shape!=(N,N) or target_aw_variance_upper<0:
        raise ValueError("valid AW floor required")
    m=[list(x) for x in p.mid]; r=[list(x) for x in p.rad]
    for i in range(off_aw,off_aw+3):
        # Enclose max(current,target) without assuming which branch fires.
        lo=m[i][i]-r[i][i]; hi=max(m[i][i]+r[i][i],target_aw_variance_upper)
        m[i][i]=.5*(lo+hi);r[i][i]=_out(.5*(hi-lo))
    return IMat(tuple(tuple(x) for x in m),tuple(tuple(x) for x in r))


def accel_bias_release_event(p: IMat, target_variance: float,
                             off_ba: int=18) -> IMat:
    """Literal release floor max(P_ba, sigma_bacc0^2), cross terms unchanged."""
    if p.shape!=(N,N) or target_variance<0:
        raise ValueError("valid bias release floor required")
    m=[list(x) for x in p.mid]; r=[list(x) for x in p.rad]
    for i in range(off_ba,off_ba+3):
        lo=max(0.0,m[i][i]-r[i][i],target_variance)
        hi=max(m[i][i]+r[i][i],target_variance)
        m[i][i]=.5*(lo+hi);r[i][i]=_out(.5*(hi-lo))
    return IMat(tuple(tuple(x) for x in m),tuple(tuple(x) for x in r))


def spectral_box(p: IMat) -> tuple[float,float]:
    if p.shape!=(N,N): raise ValueError("21-state covariance required")
    return symmetric_interval_gershgorin(p.mid,p.rad)


def diagonal_interval(values_lo: tuple[float,...], values_hi: tuple[float,...]) -> IMat:
    if len(values_lo)!=len(values_hi): raise ValueError("equal diagonal lengths required")
    n=len(values_lo); mid=[];rad=[]
    for i,(lo,hi) in enumerate(zip(values_lo,values_hi)):
        if not all(math.isfinite(x) for x in (lo,hi)) or hi<lo:
            raise ValueError("ordered finite diagonal intervals required")
        mr=[0.0]*n;rr=[0.0]*n
        mr[i]=.5*(lo+hi);rr[i]=_out(.5*(hi-lo))
        mid.append(tuple(mr));rad.append(tuple(rr))
    return IMat(tuple(mid),tuple(rad))


def selector_h(offset: int) -> IMat:
    """Exact 3x21 selector used by the integral S update."""
    if offset<0 or offset+3>N: raise ValueError("selector outside state")
    rows=[]
    for a in range(3):
        row=[0.0]*N;row[offset+a]=1.0;rows.append(row)
    return IMat(tuple(tuple(x) for x in rows),
                tuple(tuple(0.0 for _ in row) for row in rows))


def shipping_integral_update_intervals(r_s_std_min: float,
                                       r_s_std_max: float) -> tuple[IMat,IMat]:
    """Literal H/R interval constructor for applyIntegralZeroPseudoMeas."""
    if not all(math.isfinite(x) for x in (r_s_std_min,r_s_std_max)) or not 0<r_s_std_min<=r_s_std_max:
        raise ValueError("ordered positive S-noise std required")
    h=selector_h(12)
    r=diagonal_interval((r_s_std_min**2,)*3,(r_s_std_max**2,)*3)
    return h,r


def shipping_mag_update_intervals(field_norm_max: float,
                                  mag_noise_std_min: float,
                                  mag_noise_std_max: float) -> tuple[IMat,IMat]:
    """Conservative literal H/R box for J_att=-[R_wb B]x.

    Every skew entry is in [-|B|,|B|], with exact zeros on the diagonal of the
    3x3 attitude block. The remaining 18 columns are exactly zero.
    """
    vals=(field_norm_max,mag_noise_std_min,mag_noise_std_max)
    if not all(math.isfinite(x) for x in vals) or field_norm_max<=0 or not 0<mag_noise_std_min<=mag_noise_std_max:
        raise ValueError("valid magnetic interval data required")
    hm=[];hr=[]
    for i in range(3):
        mr=[0.0]*N;rr=[0.0]*N
        for j in range(3):
            if i!=j: rr[j]=field_norm_max
        hm.append(tuple(mr));hr.append(tuple(rr))
    r=diagonal_interval((mag_noise_std_min**2,)*3,(mag_noise_std_max**2,)*3)
    return IMat(tuple(hm),tuple(hr)),r


def shipping_acc_update_intervals(specific_force_norm_max: float,
                                  accel_noise_std_min: float,
                                  accel_noise_std_max: float) -> tuple[IMat,IMat]:
    """Conservative active-A21 H/R box for [J_att,0,0,0,J_aw,I_ba].

    J_att=-[f]x has off-diagonal magnitude <=|f|. J_aw=R_wb is enclosed
    entrywise in [-1,1]; J_ba=I. Lever-arm J_bg is absent in the commissioned
    theorem scope where lever arm is disabled.
    """
    vals=(specific_force_norm_max,accel_noise_std_min,accel_noise_std_max)
    if not all(math.isfinite(x) for x in vals) or specific_force_norm_max<=0 or not 0<accel_noise_std_min<=accel_noise_std_max:
        raise ValueError("valid accelerometer interval data required")
    hm=[];hr=[]
    for i in range(3):
        mr=[0.0]*N;rr=[0.0]*N
        for j in range(3):
            if i!=j: rr[j]=specific_force_norm_max
        for j in range(3): rr[15+j]=1.0
        mr[18+i]=1.0
        hm.append(tuple(mr));hr.append(tuple(rr))
    r=diagonal_interval((accel_noise_std_min**2,)*3,(accel_noise_std_max**2,)*3)
    return IMat(tuple(hm),tuple(hr)),r
