"""Source-uniform literal factor boxes for the OU-III rank-loss certificate.

The constructors here are deliberately Loewner-safe.  They use the shipping
parameter ranges and analytic integrated-OU transition, while process
covariances are dominated by trace(Q) I factors.  This is conservative but
preserves a valid common-source upper model without interval-Cholesky of an
entrywise covariance box.

No carried word is an input.
"""
from __future__ import annotations
import math
from .interval_riccati_21 import IMat
from .rank_loss_interval_factor import exact,schur_scalar_information


def _out(x): return math.nextafter(x,math.inf)


def scalar_interval(lo,hi):
    if not (math.isfinite(lo) and math.isfinite(hi) and lo<=hi):
        raise ValueError("ordered finite interval")
    return .5*(lo+hi),_out(.5*(hi-lo))


def _expm1_neg_bounds(xlo,xhi):
    # expm1(-x) is decreasing for x>=0.
    return math.expm1(-xhi),math.expm1(-xlo)


def psi_bounds(tlo,thi,taulo,tauhi):
    """Outward enclosure of tau^3*(.5 x^2-x-expm1(-x)), x=t/tau.

    psi is the integral of a positive impulse response, hence increasing in t.
    For tau dependence use a small interval subdivision in log(tau); callers
    may split further.  Endpoint/corner hull is enlarged by a derivative
    Lipschitz bound obtained from the causal integral psi=int_0^t
    (t-s)(1-exp(-s/tau)) ds: |d psi/d tau| <= t^3/(6 tau).
    """
    if not 0<=tlo<=thi or not 0<taulo<=tauhi: raise ValueError("positive box")
    def psi(t,tau):
        x=t/tau
        return tau**3*(.5*x*x-x-math.expm1(-x))
    vals=[psi(t,ta) for t in (tlo,thi) for ta in (taulo,tauhi)]
    lo=min(vals);hi=max(vals)
    # Rigorous dependency padding over tau interval around corner hull.
    pad=(thi**3/(6*taulo))*(tauhi-taulo)
    return math.nextafter(max(0.0,lo-pad),-math.inf),math.nextafter(hi+pad,math.inf)


def four_s_geometry_box(time_cells,tau_box):
    if len(time_cells)!=4: raise ValueError("four S cells required")
    tm=[];tr=[]; pm=[];pr=[]
    for lo,hi in time_cells:
        m,r=scalar_interval(lo,hi);tm.append(m);tr.append(r)
        plo,phi=psi_bounds(lo,hi,*tau_box)
        m,r=scalar_interval(plo,phi);pm.append(m);pr.append(r)
    rows=[];rads=[]
    for m,r,p,q in zip(tm,tr,pm,pr):
        # t^2/2 interval, monotone on t>=0.
        tlo=m-r;thi=m+r
        qlo=.5*tlo*tlo;qhi=.5*thi*thi
        qm,qr=scalar_interval(qlo,qhi)
        rows.append((1.0,m,qm,p));rads.append((0.0,r,qr,q))
    return IMat(tuple(rows),tuple(rads))


def four_s_residual_covariance_upper(*,time_horizon_s=1.05,
        sigma_aw_max=4.0,tau_min=.02,
        integral_noise_std_max=100.0,
        root_S_upper=1100.0,root_p_upper=8.1,root_v_upper=5.5,
        source_defect=1.001):
    """Loewner upper factor for the four S residual vector.

    This deliberately overbounds all common LIN/AW root/process uncertainty by
    one scalar lambda_max ceiling.  For any residual covariance R with
    trace(R)<=T, R<=T I.  Root amplitudes are theorem storage/physical ceilings,
    not carried covariances.  The OU process S impulse is <=t^3/6, so its
    variance is bounded by q_c int_0^T (s^3/6)^2 ds=q_c T^7/252.
    """
    T=time_horizon_s
    qc=2*sigma_aw_max**2/tau_min
    process=source_defect*qc*T**7/252.0
    # S residual at any selected epoch: S0+t p0+t^2/2 v0 + process + local noise.
    amp=root_S_upper+T*root_p_upper+.5*T*T*root_v_upper
    per=amp*amp+process+integral_noise_std_max**2
    # trace of 4x4 covariance <= four times max diagonal.
    lam=4.0*per
    return exact([[lam if i==j else 0.0 for j in range(4)] for i in range(4)])


def four_s_gamma_box(time_cells,tau_box,**kwargs):
    """Fail-closed source-uniform gamma_S box using Loewner residual ceiling."""
    v=four_s_geometry_box(time_cells,tau_box)
    v0=IMat(tuple(tuple(row[j] for j in range(3)) for row in v.mid),
            tuple(tuple(row[j] for j in range(3)) for row in v.rad))
    va=IMat(tuple((row[3],) for row in v.mid),
            tuple((row[3],) for row in v.rad))
    r=four_s_residual_covariance_upper(**kwargs)
    return schur_scalar_information(v0,va,r)


def default_four_s_cells():
    return ((0.0,.15),(.30,.45),(.60,.75),(.90,1.05))


def adaptive_four_s_gamma(*,tau_box=(.02,12.0),time_cells=None,max_depth=10,
                          target_width=.05):
    """Bisect tau/time cells until every leaf certifies or depth is exhausted."""
    if time_cells is None: time_cells=default_four_s_cells()
    leaves=[];best=math.inf
    def width(box): return box[1]-box[0]
    def rec(cells,tau,depth):
        nonlocal best
        z=four_s_gamma_box(cells,tau)
        if z["verified"]:
            best=min(best,z["lower"]);leaves.append({"verified":True,"lower":z["lower"],
                                                     "tau":tau,"cells":cells});return
        if depth>=max_depth:
            leaves.append({"verified":False,"lower":0.0,"tau":tau,"cells":cells,
                           "reason":z.get("reason","nonpositive interval Schur")});return
        # Split the widest normalized variable.  Tau is split geometrically.
        tw=math.log(tau[1]/tau[0])
        cwidth=[width(x)/.15 for x in cells]
        if tw>=max(cwidth):
            mid=math.sqrt(tau[0]*tau[1])
            rec(cells,(tau[0],mid),depth+1);rec(cells,(mid,tau[1]),depth+1)
        else:
            j=max(range(4),key=lambda i:cwidth[i]);lo,hi=cells[j];mid=.5*(lo+hi)
            a=list(cells);b=list(cells);a[j]=(lo,mid);b[j]=(mid,hi)
            rec(tuple(a),tau,depth+1);rec(tuple(b),tau,depth+1)
    rec(tuple(time_cells),tuple(tau_box),0)
    ok=all(x["verified"] for x in leaves)
    return {"verified":ok,"lower":best if ok else 0.0,"leaves":leaves,
            "leaf_count":len(leaves),"qualification":"OU3_FOUR_S_INTERVAL_V1"}


def magnetic_service_factor_box(events=None,ehb_indices=None,nuisance_indices=None):
    """Residualized magnetic certificate from literal interval event boxes.

    Each event is (H_m,S_m_actual,Phi_from_window_root), all as IMat.  If the
    boxes are absent, fail closed: the service floor alone is insufficient.
    """
    if events is None:
        return {"verified":False,"gamma_M_lower":0.0,
                "unshorted_service_mu":1.0,"window_s":1.0,
                "reason":"literal H_m/S_m_actual/Phi interval events not supplied"}
    if ehb_indices is None or nuisance_indices is None:
        raise ValueError("explicit E_hb and nuisance root columns required")
    from .magnetic_nuisance_interval import service_schur_information
    return service_schur_information(events,ehb_indices,nuisance_indices)


def rank_loss_factor_certificate(max_depth=10):
    s=adaptive_four_s_gamma(max_depth=max_depth)
    m=magnetic_service_factor_box()
    beta=0.0
    return {"qualification":"OU3_RANK_LOSS_LITERAL_BOX_V1",
            "four_S":{k:v for k,v in s.items() if k!="leaves"},
            "four_S_unresolved_boxes":sum(not x["verified"] for x in s["leaves"]),
            "magnetic":m,"gamma_S_lower":s["lower"],
            "gamma_M_lower":m["gamma_M_lower"],"beta_lower":beta,
            "source_uniform_verified":bool(s["verified"] and m["verified"] and beta>0),
            "theorem_closed":False}


def magnetic_certificate_from_operation_stratum(root_dim,ops,ehb_indices,nuisance_indices):
    """End-to-end literal operation boxes -> magnetic residualized certificate."""
    from .magnetic_literal_box_export import export_magnetic_event_boxes,event_tuples,validate_one_second_stratum
    events=export_magnetic_event_boxes(root_dim,ops)
    v=validate_one_second_stratum(events)
    if not v["verified"]:
        return {"verified":False,"gamma_M_lower":0.0,"reason":v["reason"]}
    z=magnetic_service_factor_box(event_tuples(events),ehb_indices,nuisance_indices)
    z["literal_event_count"]=len(events);z["event_stratum"]=v
    return z
