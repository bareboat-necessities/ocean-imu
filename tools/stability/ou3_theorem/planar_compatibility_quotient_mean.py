"""Covariance-metric audit of the planar homogeneous error-factor word.

FINITE DIAGNOSTIC ONLY. F, (I-KH), G are read from one shipping execution.
Their product is NOT the complete nonlinear mean/covariance derivative: in
particular dK*r, mean-dependent coefficients and exact injection differentials
are absent. Neither its spectra nor a removed physical tangent prove a tube.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from .planar_service_stream import records,expand
from .planar_service_audit import reset_matrix
from .compatibility_quotient import quotient_terminal_map

SCOPE="same-history F/(I-KH)/G homogeneous error factors; not complete nonlinear Jacobian"
G=9.80665

def word(path,start,end):
    if not (0<=start<end): raise ValueError("nonempty ordered word required")
    M=np.eye(21); P0=PN=None; pending=None; awaiting_reset=False; awaiting_post=False
    predictions=samples=corrections=resets=0; last_sample=start
    for kind,k,a in records(path):
        if k<=start: continue
        if k>end: break
        if kind==1:
            if k!=start+predictions+1 or pending is not None or awaiting_reset or awaiting_post:
                raise ValueError("incomplete prediction/correction chronology")
            if P0 is None: P0=expand(a[:225])
            M=expand(a[225:450])@M; predictions+=1
        elif kind==9:
            if pending is not None or awaiting_reset or awaiting_post:
                raise ValueError("unconsumed correction/reset")
            pending=(k,a)
        elif kind in (2,3,4):
            if pending is None or pending[0]!=k: raise ValueError("missing same-operation tangent")
            z=pending[1]
            for u,v in ((z[:288],a[:288]),(z[288:297],a[297:306]),(z[297:360],a[306:369]),(z[360:363],a[432:435])):
                if not np.array_equal(u,v): raise ValueError("tangent/correction operands differ")
            H=z[225:288].reshape(3,21); K=z[297:360].reshape(21,3)
            M=(np.eye(21)-K@H)@M
            pending=None; awaiting_reset=True; corrections+=1
        elif kind==5:
            if not awaiting_reset: raise ValueError("unpaired reset")
            M=reset_matrix(a[225:228])@M
            awaiting_reset=False; awaiting_post=True; resets+=1
        elif kind==8:
            if not awaiting_post: raise ValueError("unpaired post-reset")
            awaiting_post=False
        elif kind==7:
            if pending is not None or awaiting_reset or awaiting_post or k!=last_sample+1:
                raise ValueError("incomplete sample boundary")
            PN=expand(a[:225]); samples+=1; last_sample=k
    if (P0 is None or PN is None or pending is not None or awaiting_reset or awaiting_post
            or predictions!=end-start or samples!=predictions or corrections!=resets):
        raise ValueError("incomplete word")
    return M,P0,PN

def compatibility_line(sample):
    """First-order physical roll/BA tangent in the shipping left W->B chart.

    R_bw=Rx(alpha)Ry(psi); dR_wb R_wb'=-[Ry(-psi)e_x] d(alpha).
    db_a=g e_y d(alpha). Choosing the opposite orientation gives this line.
    At phase zero H_acc=[-[-g e_z]x, ..., I_ba] annihilates it.
    This is a tangent at alpha=0, not a nonlinear fibre/admission certificate.
    """
    psi=.02*math.sin(math.pi*(sample/200.)/10.)
    r=np.zeros(21); r[:3]=[math.cos(psi),0.,math.sin(psi)]; r[19]=-G
    return r

def legacy_line():
    r=np.zeros(21); r[0]=1.; r[19]=G
    return r

def high_precision_gains(M,P0,PN,r0,rN,dps=60):
    """High precision endpoint algebra on binary64 M, NOT an interval bound.

    The accumulated product and source binary32 arithmetic are not enclosed.
    Projectors avoid selecting a numerical complement or extra gauge vectors.
    """
    import mpmath as mp
    with mp.workdps(dps):
        A=mp.matrix(M.tolist()); L0=mp.cholesky(mp.matrix(P0.tolist())); LN=mp.cholesky(mp.matrix(PN.tolist()))
        T=LN**-1*A*L0
        u=L0**-1*mp.matrix(r0.tolist()); u=u/mp.norm(u)
        v=LN**-1*mp.matrix(rN.tolist()); v=v/mp.norm(v)
        N0=mp.eye(len(r0))-u*u.T; NN=mp.eye(len(rN))-v*v.T
        B=NN*T*N0
        return {"decimal_digits":dps,"full_metric_gain":mp.nstr(mp.sqrt(mp.eigsy(T.T*T,eigvals_only=True)[len(r0)-1]),40),
                "quotient_gain":mp.nstr(mp.sqrt(mp.eigsy(B.T*B,eigvals_only=True)[len(r0)-1]),40),
                "gauge_injection_norm":mp.nstr(mp.norm(NN*T*u),40),
                "scope":"endpoint algebra on rounded accumulated M; no directed rounding or cell bound"}

def diagnostic(path,start=40000,end=44000,high_precision=True):
    M,P0,PN=word(path,start,end); r0=compatibility_line(start); rN=compatibility_line(end)
    q=quotient_terminal_map(P0,r0,PN,rN,M,np.zeros(21))
    old=quotient_terminal_map(P0,legacy_line(),PN,legacy_line(),M,np.zeros(21))
    s=np.linalg.svd(q["M_Q"],compute_uv=False)
    L0=np.linalg.cholesky(P0); LN=np.linalg.cholesky(PN); T=np.linalg.solve(LN,M@L0)
    left,se,right=np.linalg.svd(M)
    groups=("attitude","BG","v","p","S","aw","BA")
    directions=[]
    for j in range(int(np.sum(se>1))):
        directions.append({"gain":float(se[j]),
            "right_group_norms":{g:float(np.linalg.norm(right[j,3*i:3*i+3])) for i,g in enumerate(groups)},
            "left_group_norms":{g:float(np.linalg.norm(left[3*i:3*i+3,j])) for i,g in enumerate(groups)},
            "right_coordinates":right[j].tolist(),"left_coordinates":left[:,j].tolist()})
    out={"qualification":"OU3_PLANAR_COMPATIBILITY_QUOTIENT_MEAN_V2",
        "result_type":"FINITE DIAGNOSTIC ONLY","word_samples":[start,end],"operator_scope":SCOPE,
        "stream_sha256":hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "compatibility_line_definition":"theta=Ry(-psi)e_x, b_a_y=-g in shipping left world-to-body chart",
        "compatibility_line_at_root":r0.tolist(),"compatibility_line_at_end":rN.tolist(),
        "quotient_dimension":q["quotient_dimension"],"quotient_gain_2norm":float(s[0]),
        "quotient_leading_singular_values":s[:8].tolist(),
        "gauge_to_transverse_injection_norm":q["gauge_injection_norm"],
        "legacy_positive_BAy_seed":{"quotient_gain":float(np.linalg.norm(old["M_Q"],2)),"gauge_injection_norm":old["gauge_injection_norm"],
            "physical_chart_conversion_verified":False},
        "full_euclidean_gain":float(se[0]),"full_covariance_metric_gain":float(np.linalg.norm(T,2)),
        "expanding_euclidean_directions":directions,
        "finite_expansion_classification":"metric mismatch for this homogeneous word; no instability or nonlinear contraction conclusion",
        "line_source":"physical MOVING family and shipping attitude chart; not SVD selection",
        "retained_transverse_equation":"xi_N=M_Q xi_0+b_Q+C_Q alpha_0",
        "physical_gauge_amplitude_bound":None,"uniform_b_Q":None,"uniform_c_Q":None,"uniform_q_P":None,"uniform_q_Q":None,
        "complete_nonlinear_mean_derivative_computed":False,
        "source_uniform_quotient_action_certified":False,"joint_cell_forward_invariant":False,
        "all_time_magnetic_service_verified":False,"theorem_closed":False}
    if high_precision: out["high_precision_endpoint_check"]=high_precision_gains(M,P0,PN,r0,rN)
    return out

def main():
    p=argparse.ArgumentParser(); p.add_argument("--stream",type=Path,required=True); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); o=diagnostic(a.stream); a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n"); print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
