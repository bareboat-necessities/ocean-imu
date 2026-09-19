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


def shipping_prediction_intervals(*, dt_min: float, dt_max: float,
                                  tau_min: float, tau_max: float,
                                  omega_max: float, tau_bacc: float,
                                  gyro_white_density: float,
                                  gyro_bias_rw_density: float,
                                  aw_sigma_max: float,
                                  accel_bias_drive_density: float) -> tuple[IMat,IMat]:
    """Conservative literal A21 F/Q box for the shipping 21-state predictor.

    F follows the exact block structure used by time_update:
      [F_AA,0,0; 0,F_LL,0; 0,0,phi_b I],
    with all propagated cross-covariances handled by P-=FPF'+Q.
    Rotation entries are enclosed in [-1,1], Bstep by |B|<=dt_max.
    The integrated-OU influence coefficients are bounded by their causal
    integrals dt,dt^2/2,dt^3/6, uniformly in positive tau.

    Q uses rigorous PSD diagonal ceilings. For the OU chain,
    q_c=2 sigma_aw^2/tau <= q_c,max and the impulse-response bounds give
    q_c*(dt^3/3,dt^5/20,dt^7/252,dt) for (v,p,S,a). Off-diagonal entries use
    |Qij|<=sqrt(Qii Qjj). This is conservative but source-uniform.
    """
    vals=(dt_min,dt_max,tau_min,tau_max,omega_max,tau_bacc,
          gyro_white_density,gyro_bias_rw_density,aw_sigma_max,
          accel_bias_drive_density)
    if not all(math.isfinite(x) for x in vals) or min(vals) <= 0:
        raise ValueError("positive finite prediction bounds required")
    if dt_max<dt_min or tau_max<tau_min:
        raise ValueError("ordered timing/tau bounds required")

    fm=[[0.0]*N for _ in range(N)]; fr=[[0.0]*N for _ in range(N)]
    # Attitude rotation Rstep=Exp(-omega^ dt).  Do not throw away the
    # near-identity structure by enclosing every entry in [-1,1]: Rodrigues
    # gives |R_ii-1|<=1-cos(theta_max), |R_ij|<=sin(theta_max)+1-cos(theta_max).
    # This is the first mathematically justified dependency refinement because
    # theta_max=Omega_max*dt_max is only ~3.7e-3 rad.
    theta_max=omega_max*dt_max
    diag_rad=1.0-math.cos(theta_max)
    off_rad=math.sin(theta_max)+diag_rad
    for i in range(3):
        fm[i][i]=1.0;fr[i][i]=_out(diag_rad)
        for j in range(3):
            if i!=j: fr[i][j]=_out(off_rad)
    # Bstep=-int_0^dt Exp(-omega^ s) ds.  Its diagonal is near -dt and
    # off-diagonals are O(Omega dt^2); preserve that structure as an interval.
    dt_mid=.5*(dt_min+dt_max);dt_rad=.5*(dt_max-dt_min)
    b_off=.5*omega_max*dt_max*dt_max + (omega_max**2)*dt_max**3/6.0
    b_diag_extra=(omega_max**2)*dt_max**3/6.0
    for i in range(3):
        fm[i][3+i]=-dt_mid;fr[i][3+i]=_out(dt_rad+b_diag_extra)
        for j in range(3):
            if i!=j: fr[i][3+j]=_out(b_off)
    for i in range(3): fm[3+i][3+i]=1.0

    # Linear [v,p,S,a_w], group-first offsets.
    for a in range(3):
        v,p,s,aw=6+a,9+a,12+a,15+a
        fm[v][v]=fm[p][p]=fm[s][s]=1.0
        fm[p][v]=.5*(dt_min+dt_max);fr[p][v]=_out(.5*(dt_max-dt_min))
        lo2=.5*dt_min*dt_min;hi2=.5*dt_max*dt_max
        fm[s][v]=.5*(lo2+hi2);fr[s][v]=_out(.5*(hi2-lo2))
        fm[s][p]=.5*(dt_min+dt_max);fr[s][p]=_out(.5*(dt_max-dt_min))
        # OU homogeneous factor.
        plo=math.exp(-dt_max/tau_min);phi=math.exp(-dt_min/tau_max)
        fm[aw][aw]=.5*(plo+phi);fr[aw][aw]=_out(.5*(phi-plo))
        # Positive causal a_w influence coefficients.
        for row,hi in ((v,dt_max),(p,.5*dt_max**2),(s,dt_max**3/6.0)):
            fm[row][aw]=.5*hi;fr[row][aw]=_out(.5*hi)

    # Active accelerometer-bias OU factor.
    blo=math.exp(-dt_max/tau_bacc);bhi=math.exp(-dt_min/tau_bacc)
    for a in range(3):
        fm[18+a][18+a]=.5*(blo+bhi);fr[18+a][18+a]=_out(.5*(bhi-blo))

    # Q diagonal ceilings, then PSD Cauchy-Schwarz cross bounds inside blocks.
    qdiag=[0.0]*N
    qg=gyro_white_density**2; qbg=gyro_bias_rw_density**2
    for a in range(3):
        qdiag[a]=qg*dt_max+qbg*dt_max**3/3.0
        qdiag[3+a]=qbg*dt_max
    qc=2.0*aw_sigma_max**2/tau_min
    for a in range(3):
        qdiag[6+a]=qc*dt_max**3/3.0
        qdiag[9+a]=qc*dt_max**5/20.0
        qdiag[12+a]=qc*dt_max**7/252.0
        qdiag[15+a]=qc*dt_max
        qdiag[18+a]=accel_bias_drive_density**2*dt_max
    qm=[[0.0]*N for _ in range(N)];qr=[[0.0]*N for _ in range(N)]
    for i in range(N):
        qm[i][i]=.5*qdiag[i];qr[i][i]=_out(.5*qdiag[i])
    # Actual Q is block diagonal AA/LL/BA. Bound possible within-block cross terms.
    blocks=(range(0,6),range(6,18),range(18,21))
    for block in blocks:
        inds=list(block)
        for i in inds:
            for j in inds:
                if i!=j: qr[i][j]=_out(math.sqrt(qdiag[i]*qdiag[j]))
    return IMat(tuple(tuple(x) for x in fm),tuple(tuple(x) for x in fr)), IMat(
        tuple(tuple(x) for x in qm),tuple(tuple(x) for x in qr))


def _inverse3_exact(a: tuple[tuple[float,...],...]) -> tuple[tuple[float,...],...]:
    if len(a)!=3 or any(len(r)!=3 for r in a): raise ValueError("3x3 required")
    x=a
    det=(x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1])
         -x[0][1]*(x[1][0]*x[2][2]-x[1][2]*x[2][0])
         +x[0][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0]))
    if not math.isfinite(det) or det==0: raise ValueError("singular midpoint")
    adj=((x[1][1]*x[2][2]-x[1][2]*x[2][1],
          x[0][2]*x[2][1]-x[0][1]*x[2][2],
          x[0][1]*x[1][2]-x[0][2]*x[1][1]),
         (x[1][2]*x[2][0]-x[1][0]*x[2][2],
          x[0][0]*x[2][2]-x[0][2]*x[2][0],
          x[0][2]*x[1][0]-x[0][0]*x[1][2]),
         (x[1][0]*x[2][1]-x[1][1]*x[2][0],
          x[0][1]*x[2][0]-x[0][0]*x[2][1],
          x[0][0]*x[1][1]-x[0][1]*x[1][0]))
    return tuple(tuple(v/det for v in row) for row in adj)


def verified_inverse3_interval(s: IMat) -> tuple[IMat,dict]:
    """Entrywise inverse enclosure after the strict spectral residual check."""
    if s.shape!=(3,3): raise ValueError("3x3 innovation required")
    cert=innovation_inverse_spectral_certificate(s)
    if not cert["verified"]: raise ValueError("innovation inverse not verified")
    inv0=_inverse3_exact(s.mid)
    # ||S0^-1|| <= 1/a from the same certified midpoint eigen floor.
    a=1.0/cert["inverse_norm_bound"] + max(sum(row) for row in s.rad)
    er=max(sum(row) for row in s.rad)
    residual=er/a
    delta=(1.0/a)*residual/(1.0-residual)
    rad=tuple(tuple(_out(delta) for _ in range(3)) for _ in range(3))
    return IMat(inv0,rad),cert


def verified_gain_interval(p: IMat,h: IMat,r: IMat) -> tuple[IMat,dict]:
    """K=P H' S^-1 with a verified interval inverse of S."""
    s=innovation_covariance(p,h,r)
    sinv,cert=verified_inverse3_interval(s)
    pct=matmul(p,transpose(h))
    return matmul(pct,sinv),cert


def verified_joseph_update(p: IMat,h: IMat,r: IMat) -> tuple[IMat,dict]:
    """Complete verified 3-D shipping covariance correction."""
    k,cert=verified_gain_interval(p,h,r)
    return joseph_covariance(p,k,h,r),cert


def contains(outer: IMat, inner: IMat) -> bool:
    if outer.shape!=inner.shape: return False
    n,m=outer.shape
    for i in range(n):
        for j in range(m):
            olo=outer.mid[i][j]-outer.rad[i][j];ohi=outer.mid[i][j]+outer.rad[i][j]
            ilo=inner.mid[i][j]-inner.rad[i][j];ihi=inner.mid[i][j]+inner.rad[i][j]
            if ilo<olo or ihi>ohi: return False
    return True


def hull(a: IMat,b: IMat) -> IMat:
    if a.shape!=b.shape: raise ValueError("shape mismatch")
    n,m=a.shape;cm=[];cr=[]
    for i in range(n):
        mr=[];rr=[]
        for j in range(m):
            lo=min(a.mid[i][j]-a.rad[i][j],b.mid[i][j]-b.rad[i][j])
            hi=max(a.mid[i][j]+a.rad[i][j],b.mid[i][j]+b.rad[i][j])
            mr.append(.5*(lo+hi));rr.append(_out(.5*(hi-lo)))
        cm.append(tuple(mr));cr.append(tuple(rr))
    return IMat(tuple(cm),tuple(cr))


def inflate(a: IMat, relative: float, absolute: float=0.0) -> IMat:
    if not math.isfinite(relative) or relative<0 or not math.isfinite(absolute) or absolute<0:
        raise ValueError("nonnegative inflation required")
    r=[]
    for mr,rr in zip(a.mid,a.rad):
        r.append(tuple(_out(x+relative*abs(m)+absolute) for m,x in zip(mr,rr)))
    return IMat(a.mid,tuple(r))


def iterate_recurring_box(seed: IMat, word_map, *, max_iterations: int=32,
                          hull_inflation: float=.02) -> dict:
    """Iterate a verified interval word map to a self-containing covariance box.

    word_map must use interval kernels and verified innovation solves and return
    (image, innovation_certificates). No midpoint trajectory can promote this
    routine. Completion requires a final fresh image contained entrywise in the
    proposed box and a strictly positive interval spectral lower bound.
    """
    if seed.shape!=(N,N) or max_iterations<1:
        raise ValueError("valid 21-state seed/iteration count required")
    box=seed
    for it in range(max_iterations):
        image,certs=word_map(box)
        if any(not x.get("verified",False) for x in certs):
            return {"verified":False,"reason":"innovation inverse","iterations":it+1,
                    "box":box,"innovation_certificates":certs}
        if contains(box,image):
            lo,hi=spectral_box(box)
            return {"verified":lo>0.0,"reason":"contained" if lo>0 else "nonpositive spectral floor",
                    "iterations":it+1,"box":box,"image":image,
                    "spectral_lower":lo,"spectral_upper":hi,
                    "innovation_certificates":certs}
        box=inflate(hull(box,image),hull_inflation)
    image,certs=word_map(box)
    lo,hi=spectral_box(box)
    return {"verified":contains(box,image) and lo>0 and all(x.get("verified",False) for x in certs),
            "reason":"iteration limit","iterations":max_iterations,"box":box,"image":image,
            "spectral_lower":lo,"spectral_upper":hi,
            "innovation_certificates":certs}


def shipping_max_correction_step(p: IMat, *, f: IMat, q: IMat,
                                 h_s: IMat, r_s: IMat,
                                 h_acc: IMat, r_acc: IMat,
                                 h_mag: IMat, r_mag: IMat) -> tuple[IMat,list[dict]]:
    """Conservative maximal-correction A21 sample.

    Literal wrapper order is time_update (which may apply S), then accelerometer;
    updateMag is asynchronous afterward. For a covariance lower-bound search,
    applying every optional S and magnetic correction at every sample is a
    conservative maximal-information cadence. A final proof must justify this
    monotonic replacement in information form; this function does not itself
    promote that ordering argument.
    """
    out=predict_covariance(p,f,q);certs=[]
    for name,h,r in (("S",h_s,r_s),("acc",h_acc,r_acc),("mag",h_mag,r_mag)):
        out,cert=verified_joseph_update(out,h,r)
        certs.append({"name":name,**cert})
    return out,certs


def shipping_word_map(max_steps: int, **step_kwargs):
    """Return a recurring-box map built from the conservative sample step."""
    if max_steps<1: raise ValueError("positive step count required")
    def word(p: IMat):
        out=p;certs=[]
        for _ in range(max_steps):
            out,c=shipping_max_correction_step(out,**step_kwargs)
            certs.extend(c)
        return out,certs
    return word


def split_interval_matrix(a: IMat, entries: tuple[tuple[int,int],...]) -> tuple[IMat,...]:
    """Bisect selected interval entries, preserving all correlations not split.

    This is a dependency-control primitive for the verified Riccati enclosure.
    It never shrinks the union: the returned cells exactly cover the original
    entrywise box.  Use only a small set of dominant attitude/tuner entries per
    branch to avoid exponential explosion.
    """
    cells=[a]
    for i,j in entries:
        nxt=[]
        for x in cells:
            if not (0<=i<x.shape[0] and 0<=j<x.shape[1]):
                raise ValueError("split entry outside matrix")
            r=x.rad[i][j]
            if r==0.0:
                nxt.append(x);continue
            lo=x.mid[i][j]-r;hi=x.mid[i][j]+r;cut=.5*(lo+hi)
            for aa,bb in ((lo,cut),(cut,hi)):
                m=[list(row) for row in x.mid];q=[list(row) for row in x.rad]
                m[i][j]=.5*(aa+bb);q[i][j]=_out(.5*(bb-aa))
                nxt.append(IMat(tuple(tuple(row) for row in m),
                                tuple(tuple(row) for row in q)))
        cells=nxt
    return tuple(cells)


def branch_verified_joseph_update(p: IMat, h_cells: tuple[IMat,...],
                                  r: IMat) -> tuple[IMat,list[dict]]:
    """Union enclosure of Joseph updates over a finite H subdivision."""
    if not h_cells: raise ValueError("nonempty H subdivision required")
    outs=[];certs=[]
    for idx,h in enumerate(h_cells):
        out,cert=verified_joseph_update(p,h,r)
        outs.append(out);certs.append({"cell":idx,**cert})
    union=outs[0]
    for out in outs[1:]: union=hull(union,out)
    return union,certs


def adaptive_verified_update(p: IMat,h: IMat,r: IMat, *,
                             split_entries: tuple[tuple[int,int],...],
                             max_depth: int=3) -> tuple[IMat,list[dict]]:
    """Verify an update, subdividing H only when the innovation box is too wide.

    The full unsplit cell is attempted first. On failure, dominant uncertain H
    entries are bisected one level at a time. Every cell must verify; otherwise
    the update fails closed. This directly attacks interval dependency without
    trusting a midpoint branch.
    """
    try:
        out,cert=verified_joseph_update(p,h,r)
        return out,[{"depth":0,"cell":0,**cert}]
    except ValueError:
        pass
    cells=(h,)
    used=0
    for depth in range(1,max_depth+1):
        if used>=len(split_entries): break
        cells=tuple(y for x in cells for y in split_interval_matrix(x,(split_entries[used],)))
        used+=1
        try:
            out,certs=branch_verified_joseph_update(p,cells,r)
            return out,[{"depth":depth,**x} for x in certs]
        except ValueError:
            continue
    raise ValueError("innovation inverse not verified after H subdivision")


def recurring_box_over_cells(seed: IMat, word_maps: tuple, *,
                             max_iterations: int=24,
                             hull_inflation: float=.01) -> dict:
    """Verified recurring box for a finite union of schedule/attitude cells.

    Each word map represents one outward-covered branch cell. The invariant
    image is the hull of *all* cells, so no favorable midpoint schedule can
    certify the box.
    """
    if not word_maps: raise ValueError("nonempty branch-cell family required")
    box=seed
    for it in range(max_iterations):
        images=[];certs=[]
        try:
            for idx,w in enumerate(word_maps):
                image,c=w(box);images.append(image)
                certs.extend({"branch":idx,**x} for x in c)
        except ValueError:
            return {"verified":False,"reason":"branch innovation inverse",
                    "iterations":it+1,"box":box,"innovation_certificates":certs}
        image=images[0]
        for x in images[1:]: image=hull(image,x)
        if contains(box,image):
            lo,hi=spectral_box(box)
            return {"verified":lo>0.0,"reason":"contained" if lo>0 else "nonpositive spectral floor",
                    "iterations":it+1,"box":box,"image":image,
                    "spectral_lower":lo,"spectral_upper":hi,
                    "innovation_certificates":certs,"branch_count":len(word_maps)}
        box=inflate(hull(box,image),hull_inflation)
    lo,hi=spectral_box(box)
    return {"verified":False,"reason":"iteration limit","iterations":max_iterations,
            "box":box,"spectral_lower":lo,"spectral_upper":hi,
            "branch_count":len(word_maps)}


def exact_midpoint_seed(diagonal: tuple[float,...]) -> IMat:
    """Positive diagonal proposal seed; never a certificate by itself."""
    if len(diagonal)!=N or any((not math.isfinite(x) or x<=0) for x in diagonal):
        raise ValueError("21 positive diagonal seed entries required")
    return diagonal_interval(diagonal,diagonal)


def interval_failure_metrics(p: IMat,h: IMat,r: IMat) -> dict:
    """Non-promoting diagnostic required by the research protocol."""
    s=innovation_covariance(p,h,r)
    cert=innovation_inverse_spectral_certificate(s)
    lo,hi=spectral_box(p)
    return {"innovation_residual_ratio":cert["residual_norm_bound"],
            "innovation_verified":cert["verified"],
            "covariance_spectral_lower":lo,
            "covariance_spectral_upper":hi}


def _frob_mid(a: IMat) -> float:
    return math.sqrt(sum(x*x for row in a.mid for x in row))


def _frob_rad(a: IMat) -> float:
    return math.sqrt(sum(x*x for row in a.rad for x in row))


def midpoint_sigma_min_lower(a: IMat) -> float:
    """Gershgorin lower bound on sigma_min of the midpoint matrix."""
    at=transpose(IMat(a.mid,tuple(tuple(0.0 for _ in row) for row in a.mid)))
    am=IMat(a.mid,tuple(tuple(0.0 for _ in row) for row in a.mid))
    gram=matmul(at,am)
    lo,_=symmetric_interval_gershgorin(gram.mid,gram.rad)
    return math.sqrt(max(0.0,lo))


def psd_spectral_prediction_floor(p_lower: float, f: IMat,
                                  q_lower: float=0.0) -> dict:
    """PSD-structure-preserving prediction floor.

    For F=F0+dF, ||dF||2<=||rad(F)||F, Weyl gives
    sigma_min(F)>=sigma_min(F0)-||dF||. Therefore
    F P F'+Q >= (sigma_floor^2*p_lower+q_lower) I.
    Unlike entrywise covariance propagation this cannot manufacture a negative
    covariance eigenvalue merely from dependency.
    """
    if not all(math.isfinite(x) for x in (p_lower,q_lower)) or p_lower<=0 or q_lower<0:
        raise ValueError("positive covariance floor and nonnegative Q floor required")
    smid=midpoint_sigma_min_lower(f);dr=_frob_rad(f)
    sf=max(0.0,smid-dr)
    return {"midpoint_sigma_min_lower":smid,"F_radius_norm_upper":dr,
            "F_sigma_min_lower":sf,
            "predicted_covariance_floor":sf*sf*p_lower+q_lower}


def psd_spectral_correction_floor(p_lower: float,h: IMat,
                                  r_noise_lower: float) -> dict:
    """Information-form lower covariance bound for one 3-D correction.

    P+ = (P^-1+H'R^-1H)^-1 and ||H||2 is bounded by midpoint+radius
    Frobenius norms. This is conservative but PSD preserving and valid with all
    cross-covariances present.
    """
    if not all(math.isfinite(x) for x in (p_lower,r_noise_lower)) or p_lower<=0 or r_noise_lower<=0:
        raise ValueError("positive covariance/noise floors required")
    hn=_frob_mid(h)+_frob_rad(h)
    out=1.0/(1.0/p_lower+(hn*hn)/r_noise_lower)
    return {"H_norm_upper":hn,"posterior_covariance_floor":out}


def psd_spectral_word_floor(*, initial_floor: float, f: IMat,
                            q_floor: float,
                            corrections: tuple[tuple[IMat,float],...],
                            steps: int) -> dict:
    """Iterate a PSD-preserving scalar spectral lower recurrence.

    This is rigorous once q_floor is a certified lower eigenvalue of the
    literal full 21-state process covariance for every admitted prediction.
    It is intentionally not a replacement for multi-step controllability when
    q_floor is zero or numerically useless.
    """
    if not math.isfinite(initial_floor) or initial_floor<=0 or not math.isfinite(q_floor) or q_floor<0 or steps<1:
        raise ValueError("valid spectral word data required")
    p=initial_floor; prefix=[p]
    pred_meta=None
    for _ in range(steps):
        pred_meta=psd_spectral_prediction_floor(p,f,q_floor)
        p=pred_meta["predicted_covariance_floor"];prefix.append(p)
        for h,rlo in corrections:
            meta=psd_spectral_correction_floor(p,h,rlo)
            p=meta["posterior_covariance_floor"];prefix.append(p)
    return {"root_floor":initial_floor,"end_floor":p,
            "prefix_floor":min(prefix),"prefixes":tuple(prefix),
            "F_sigma_min_lower":pred_meta["F_sigma_min_lower"]}
