#!/usr/bin/env python3
"""V11 perturbation closure for the OU-III P5 sample-1 accelerometer core.

V10 closes the canonical source-cell-0 sample-1 family at 7.017 rad by carrying
the exact combined perpendicular residual. Two covariance effects were still
kept outside that canonical algebra:

1. the cross-axis part of the small attitude/gyro-bias PSD remainder before the
   first accelerometer, plus the same source process remainder generated on the
   next 5 ms attitude prediction;
2. the sample-1 S=0 accepted covariance/reset branch (the solver-identity branch
   is the zero perturbation).

This producer keeps V10 as the nominal 1+2 block certificate and proves an
additive source-uniform perturbation budget. It does not rebuild the broad
sample-1 3x3 inverse.

V10 already includes the diagonal [t,t+eps] and [Y,Y+eps] widening. For the
remaining symmetric zero-diagonal attitude term, |E_ij|<=eps gives
||E_off||_2<=2 eps. For H0=[-[g e3]_x,I], the first innovation has the exact
floor lambda0>=p_aw+R_a, so resolvent identities bound the first attitude/a_w
gain perturbations. Those bounds explicitly charge the small off-axis first
attitude correction and the previously-zero first a_w,x correction into V10's
combined perpendicular residual.

The first covariance perturbation is transported with the exact derivative of
the Kalman covariance map,

    D Phi_P[E]=(I-KH) E (I-KH)^T,

using a source-cell bound on ||I-KH||, then through the shipping reset. A second
spectral eps is added for the next 5 ms attitude prediction.

For sample-1 S=0, a source-structured full covariance is used only to bound the
small P_[theta,aw],S cross block. Since S_S>=R_S,min I,

    ||Delta P_z|| <= ||P_zS||^2/R_S,min.

The same cross blocks and the propagated shipping S mean bound the S-induced
attitude and a_w mean corrections. The immediate S reset is added in operator
norm. The solver-identity S branch is the zero perturbation, hence is contained
by the same accepted-branch upper bound.

For the final sample-1 gain perturbation, S_acc>=R_a I is used directly. If
DeltaP is the reduced covariance perturbation and H is the nominal sample-1
Jacobian,

    ||K'-K|| <= (DeltaP||H|| + ||K|| ||H||^2 DeltaP)
                /(R_a-||H||^2 DeltaP),

whenever the denominator remains positive. V11 adds this gain perturbation and
the explicitly bounded residual perturbation to V10's directional correction.
It is fail-closed and does not yet perform the signed Cayley q<8 composition or
promote a complete source family/word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_identity, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_entry as ENTRY
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_structured_full_gain_v4 as V4
import ou3_p5_sample1_structured_full_gain_v5 as V5
import ou3_p5_sample1_structured_full_gain_v10 as V10
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN=FIRST.DEFAULT_DOMAIN
SCHEMA=11
RANGE=9.0


def _op(A)->float:
    return RG._op2_upper(A)


def _submatrix(A,rows,cols):
    return [[A[i][j] for j in cols] for i in rows]


def _first_reduced_H(g:float):
    z=FULL.I(0.0); one=FULL.I(1.0); gg=FULL.I(g)
    return [[z,gg,z,one,z,z],[-gg,z,z,z,one,z],[z,z,z,z,z,one]]


def _first_reduced_K(t:Interval,p:Interval,r:Interval,g:float):
    z=FULL.I(0.0); D=FULL.I(g*g)*t+p+r
    kt=FULL.I(g)*t/D; ka=p/D; kz=p/(p+r)
    return [[z,-kt,z],[kt,z,z],[z,z,z],[ka,z,z],[z,ka,z],[z,z,kz]]


def _first_A_norm(t:Interval,p:Interval,r:Interval,g:float)->float:
    H=_first_reduced_H(g); K=_first_reduced_K(t,p,r,g)
    A=matrix_identity(6); KH=matrix_mul(K,H)
    A=[[A[i][j]-KH[i][j] for j in range(6)] for i in range(6)]
    return _op(A)


def _first_psd_perturbation(*,t:Interval,Y:Interval,p:Interval,r:Interval,g:float,
                            eps:float,rho0:float,dhi:float,rt:Interval,rz:Interval,
                            alpha_hi:float,aw_pre:float)->dict:
    eoff=FULL.up(2.0*eps)
    lam=(p+r).lo
    if not lam>0.0: raise RuntimeError("first innovation floor lost")
    dS=FULL.up(g*g*eoff)
    lam2=FULL.down(lam*lam)
    dkth=FULL.up(FULL.up(g*eoff/lam)+FULL.up((g*t.hi)*dS/lam2))
    dkaw=FULL.up(p.hi*dS/lam2)
    dkfull=FULL.up(math.sqrt(FULL.up(dkth*dkth+FULL.up(dkaw*dkaw))))
    dd=FULL.up(dkth*rho0); daw=FULL.up(dkaw*rho0)

    D=FULL.I(g*g)*t+p+r; kawt=p/D; kz=p/(p+r)
    awt=kawt*rt; az=kz*rz
    dxaw=V5._norm2_upper(awt.abs_upper(),az.abs_upper())
    vec=FULL.up(g+FULL.up(alpha_hi*FULL.up(aw_pre+FULL.up(dxaw+daw))))
    drho=FULL.up(FULL.up(alpha_hi*daw)+FULL.up(dd*vec))

    h0=FULL.up(math.sqrt(FULL.up(g*g+1.0)))
    A0=_first_A_norm(t,p,r,g); Amax=FULL.up(A0+FULL.up(dkfull*h0))
    dPplus=FULL.up(FULL.up(Amax*Amax)*eoff)

    Tnom=FULL.up(math.sqrt(FULL.up(1.0+FULL.up(0.25*dhi*dhi))))
    Gactual=FULL.up(math.sqrt(FULL.up(1.0+FULL.up(0.25*FULL.up((dhi+dd)*(dhi+dd))))))
    dT=FULL.up(FULL.up(dd*Gactual)+FULL.up(0.5*dd))
    Pprior=max(Y.hi,t.hi,p.hi)
    dirterm=FULL.up(FULL.up(FULL.up(2.0*Tnom*dT)+FULL.up(dT*dT))*Pprior)
    after_reset=FULL.up(FULL.up(Tnom*Tnom*dPplus)+dirterm)
    dP1=FULL.up(after_reset+eps)
    return {
        "first_attitude_offdiagonal_operator_upper":eoff,
        "first_gain_theta_perturbation_upper":dkth,
        "first_gain_aw_perturbation_upper":dkaw,
        "first_offaxis_attitude_correction_upper_rad":dd,
        "first_aw_x_correction_upper_mps2":daw,
        "first_nominal_aw_correction_norm_upper_mps2":dxaw,
        "PSD_induced_sample1_residual_perturbation_upper_mps2":drho,
        "first_covariance_update_A_norm_upper":A0,
        "first_covariance_update_A_perturbed_norm_upper":Amax,
        "first_posterior_covariance_perturbation_upper":dPplus,
        "reset_gauge_transform_perturbation_upper":dirterm,
        "sample1_reduced_covariance_PSD_perturbation_upper":dP1,
    }


def _sample1_s_bounds(path:Path,src:dict,dom:dict,first:dict,vector:dict)->dict:
    h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"])
    rho0=float(first["source_cells"][0]["combined_useful_residual_norm_upper_mps2"])
    d0=float(first["source_cells"][0]["correction_norm_upper_rad"])
    H0=ENTRY._canonical_first_H(g)
    Racc=FULL._R_diag(float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"]))
    F,Q,_=FULL._transition_and_Q(src,dom)
    P0=FULL._initial_covariance(src,path)
    Pp=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,P0),matrix_transpose(F)),Q))
    Pp=ENTRY._canonicalize_first_attitude_covariance(Pp,path,h)
    Ppre,_=ENTRY._zero_residual_S_covariance(Pp,src)
    PHt,S=FULL._innovation(Ppre,H0,Racc); Sinv,backend=FULL._spd_inverse_enclosure(S,Racc); K=matrix_mul(PHt,Sinv)
    Pj=FULL._shipping_joseph(Ppre,K,S,PHt)
    db=[Interval(-FULL.up(d0),FULL.up(d0)) for _ in range(3)]
    Pr=FULL._reset_covariance(Pj,db)
    P1=FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F,Pr),matrix_transpose(F)),Q))

    Z=list(FULL.TH)+list(FULL.AW); SS=list(FULL.SS)
    PzS=_submatrix(P1,Z,SS); PtS=_submatrix(P1,list(FULL.TH),SS); PaS=_submatrix(P1,list(FULL.AW),SS); Pzz=_submatrix(P1,Z,Z)
    nzs=_op(PzS); nts=_op(PtS); nas=_op(PaS); nzz=_op(Pzz)
    rs2=src["R_S_filter_std"].square().lo
    if not rs2>0.0: raise RuntimeError("sample1 S covariance floor lost")

    FK=matrix_mul(F,K); MS=_submatrix(FK,SS,range(3)); xS=FULL.up(_op(MS)*rho0)
    ds=FULL.up(FULL.up(nts/rs2)*xS); daw=FULL.up(FULL.up(nas/rs2)*xS)
    decrement=FULL.up(FULL.up(nzs*nzs)/rs2)
    e=FULL.up(0.5*ds); reset=FULL.up(FULL.up(FULL.up(2.0*e)+FULL.up(e*e))*nzz)
    dP=FULL.up(decrement+reset)
    return {
        "first_accel_inverse_backend_for_S_bound":backend,
        "sample1_PzS_operator_upper":nzs,"sample1_PthetaS_operator_upper":nts,"sample1_PawS_operator_upper":nas,
        "sample1_reduced_P_operator_upper":nzz,"sample1_S_variance_floor":rs2,
        "sample1_estimator_S_mean_norm_upper":xS,"sample1_S_attitude_correction_upper_rad":ds,
        "sample1_S_aw_correction_upper_mps2":daw,"sample1_S_covariance_decrement_upper":decrement,
        "sample1_S_reset_covariance_perturbation_upper":reset,"sample1_S_total_reduced_covariance_perturbation_upper":dP,
    }


def _nominal_sample1_matrices(*,t:Interval,Y:Interval,p:Interval,r:Interval,g:float,
                              alpha:Interval,qaw:Interval,d:Interval,fy:Interval,fz:Interval):
    D=FULL.I(g*g)*t+p+r
    a=t*(p+r)/D; c0=-(FULL.I(g)*t*p/D); b=p*(FULL.I(g*g)*t+r)/D; bz=p*r/(p+r)
    Pth0=V4._diag3(a,a,Y); Paw0=V4._diag3(b,b,bz); C0=FULL._zero(3,3); C0[0][1]=-c0; C0[1][0]=c0
    L,Rx=V4._Ltheta(d); Pth=matrix_mul(matrix_mul(L,Pth0),matrix_transpose(L)); Paw_r=matrix_mul(matrix_mul(Rx,Paw0),matrix_transpose(Rx)); C_r=matrix_mul(matrix_mul(L,C0),matrix_transpose(Rx))
    Paw=[[alpha.square()*Paw_r[i][j]+(qaw if i==j else FULL.I(0.0)) for j in range(3)] for i in range(3)]
    C=[[alpha*C_r[i][j] for j in range(3)] for i in range(3)]
    P=FULL._zero(6,6)
    for i in range(3):
        for j in range(3):
            P[i][j]=Pth[i][j]; P[i][3+j]=C[i][j]; P[3+i][j]=C[j][i]; P[3+i][3+j]=Paw[i][j]
    Ht=V4._Htheta([FULL.I(0.0),fy,fz]); H=FULL._zero(3,6)
    for i in range(3):
        for j in range(3): H[i][j]=Ht[i][j]
        H[i][3+i]=FULL.I(1.0)
    R=FULL._zero(3,3)
    for i in range(3): R[i][i]=r
    S=FULL.matrix_symmetric_hull(matrix_add(matrix_mul(matrix_mul(H,P),matrix_transpose(H)),R))
    return P,H,S


def build(domain_path:Path=DEFAULT_DOMAIN,*,source_pieces:int=4,source_cell_index:int=0,
          p_pieces:int=24,tangent_pieces:int=24,axial_pieces:int=24)->dict:
    FULL3._install_backend(); path=Path(domain_path).resolve(); dom=json.loads(path.read_text(encoding="utf-8"))
    core=V10.build(path,source_pieces=source_pieces,source_cell_index=source_cell_index,p_pieces=p_pieces,tangent_pieces=tangent_pieces,axial_pieces=axial_pieces)
    first=FIRST.build(path,source_pieces=source_pieces); vec=VECTOR.build(); failures=[]
    failures += [f"V10: {x}" for x in V10.validate(core)] + [f"first: {x}" for x in FIRST.validate(first)] + [f"vector: {x}" for x in VECTOR.validate(vec)]
    if core.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10")!="PASS": failures.append("V10 prerequisite did not pass")
    src,phase=RG._source_phase_children(source_pieces)[source_cell_index]
    if phase!="due": failures.append("V11 focused perturbation requires first due source cell")
    if source_cell_index!=0: failures.append("V11 S perturbation helper currently certified for source cell 0")

    fr=first["source_cells"][source_cell_index]; p_all=Interval.outward_bounds(*map(float,fr["P_aw_variance_interval"]))
    rho0=float(fr["combined_useful_residual_norm_upper_mps2"]); aw_pre=float(fr["predicted_aw_error_norm_upper_mps2"])
    h=float(src["dt_s"]); g=float(dom["startup"]["gravity_mps2"])
    tilt,yaw,eps=RG._attitude_covariance_epsilon(path,h); t=Interval.outward_bounds(tilt,FULL.up(tilt+eps)); Y=Interval.outward_bounds(yaw,FULL.up(yaw+eps))
    r=FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F,Q,_=FULL._transition_and_Q(src,dom); alpha=F[15][15]; alpha_hi=alpha.hi; qaw=Q[15][15]
    pcells=SUB.parts(p_all.lo,p_all.hi,p_pieces); sb=_sample1_s_bounds(path,src,dom,first,vec)

    rows=[]; bad=None; worst=None; unclosed=0
    maxd=maxdp=maxdr=maxdk=maxlamloss=0.0; minlam=math.inf
    for base in core["rows"]:
        pi=int(base["p_cell"]); p=pcells[pi]
        rt=Interval.outward_bounds(*map(float,base["first_tangent_residual_magnitude_mps2"])); rz=Interval.outward_bounds(*map(float,base["first_axial_residual_mps2"])); d=Interval.outward_bounds(*map(float,base["first_attitude_correction_rad"]))
        D=FULL.I(g*g)*t+p+r; fy=-(alpha*(p/D)*rt); fz=FULL.I(g)+alpha*(p/(p+r))*rz
        psd=_first_psd_perturbation(t=t,Y=Y,p=p,r=r,g=g,eps=eps,rho0=rho0,dhi=max(0.0,d.hi),rt=rt,rz=rz,alpha_hi=alpha_hi,aw_pre=aw_pre)

        vecnorm=FULL.up(g+FULL.up(alpha_hi*(aw_pre+psd["first_nominal_aw_correction_norm_upper_mps2"]+psd["first_aw_x_correction_upper_mps2"])))
        drS=FULL.up(sb["sample1_S_aw_correction_upper_mps2"]+FULL.up(sb["sample1_S_attitude_correction_upper_rad"]*vecnorm))
        drho=FULL.up(psd["PSD_induced_sample1_residual_perturbation_upper_mps2"]+drS)
        dP=FULL.up(psd["sample1_reduced_covariance_PSD_perturbation_upper"]+sb["sample1_S_total_reduced_covariance_perturbation_upper"])

        Pn,Hn,Sn=_nominal_sample1_matrices(t=t,Y=Y,p=p,r=r,g=g,alpha=alpha,qaw=qaw,d=d,fy=fy,fz=fz)
        h1=_op(Hn); pnorm=_op(Pn)
        # Exact source fact: every real nominal innovation satisfies S>=R_a I.
        lam=r.lo
        dS=FULL.up(FULL.up(h1*h1)*dP); lamp=FULL.down(lam-dS)
        if not lamp>0.0:
            row={"reason":"perturbed innovation lower bound nonpositive","nominal_lambda_lower":lam,"innovation_perturbation_upper":dS,"p_cell":pi}
            rows.append(row); unclosed+=1
            if bad is None: bad=row
            continue
        k0=max(float(base["Ktheta_perpendicular_block_upper"]),float(base["Ktheta_parallel_block_upper"]))
        dC=FULL.up(dP*h1); dk=FULL.up(FULL.up(dC+FULL.up(k0*dS))/lamp)
        rho=float(base["sample1_full_residual_norm_upper_mps2"]); d10=float(base["combined_directional_correction_norm_upper_rad"])
        d11=FULL.up(d10+FULL.up(FULL.up(k0*drho)+FULL.up(dk*FULL.up(rho+drho))))
        closed=math.isfinite(d11) and d11<RANGE
        row={
            "p_cell":pi,"tangent_residual_cell":base["tangent_residual_cell"],"axial_residual_cell":base["axial_residual_cell"],
            "V10_directional_correction_upper_rad":d10,"PSD_residual_perturbation_upper_mps2":psd["PSD_induced_sample1_residual_perturbation_upper_mps2"],
            "S_residual_perturbation_upper_mps2":drS,"total_residual_perturbation_upper_mps2":drho,
            "PSD_reduced_covariance_perturbation_upper":psd["sample1_reduced_covariance_PSD_perturbation_upper"],
            "S_reduced_covariance_perturbation_upper":sb["sample1_S_total_reduced_covariance_perturbation_upper"],"total_reduced_covariance_perturbation_upper":dP,
            "nominal_reduced_covariance_operator_upper":pnorm,"sample1_H_operator_upper":h1,"nominal_innovation_lambda_lower":lam,
            "innovation_perturbation_upper":dS,"perturbed_innovation_lambda_lower":lamp,"sample1_gain_operator_perturbation_upper":dk,
            "V11_correction_norm_upper_rad":d11,"inside_9rad_range":closed,**psd,
        }
        rows.append(row); maxd=max(maxd,d11); maxdp=max(maxdp,dP); maxdr=max(maxdr,drho); maxdk=max(maxdk,dk); maxlamloss=max(maxlamloss,dS); minlam=min(minlam,lamp)
        if worst is None or d11>worst.get("V11_correction_norm_upper_rad",-math.inf): worst=row
        if not closed:
            unclosed+=1
            if bad is None: bad=row

    ok=bool(rows) and unclosed==0 and bad is None and not failures
    return {
        "schema":SCHEMA,"qualification":"OU3_P5_SAMPLE1_PSD_AND_S_PERTURBATION_CLOSURE_V11",
        "source_generated_not_trajectory_fit":True,"source_replay_used":False,"filter_changed":False,
        "V10_canonical_core_retained":True,"V10_exact_combined_perpendicular_residual_retained":True,
        "first_attitude_PSD_diagonal_remainder_already_in_V10":True,"first_attitude_PSD_cross_axis_remainder_included":True,
        "second_prediction_attitude_process_remainder_included":True,"sample1_S_covariance_update_included":True,
        "sample1_S_attitude_injection_included":True,"sample1_S_aw_mean_correction_included":True,
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation":True,"broad_sample1_3x3_interval_inverse_reintroduced":False,
        "temporal_force_slew_assumed":False,"complete_sample1_branch_closed_here":False,"signed_cayley_q8_composed_here":False,
        "q8_word_promoted_here":False,"whole_word_promoted_here":False,"N_H_words_set_here":False,
        "validated_deployed_quaternion_range_rad":RANGE,"attitude_covariance_remainder_spectral_upper":eps,
        "sample1_S_perturbation_bounds":sb,"evaluated_joint_cells":len(rows),"unclosed_joint_cells":unclosed,
        "max_total_residual_perturbation_upper_mps2":maxdr,"max_total_reduced_covariance_perturbation_upper":maxdp,
        "max_sample1_gain_operator_perturbation_upper":maxdk,"max_innovation_perturbation_upper":maxlamloss,
        "minimum_perturbed_innovation_lambda_lower":minlam,"max_V11_correction_norm_upper_rad":maxd,
        "first_unclosed_joint_cell":bad,"worst_joint_cell":worst,
        "P5_SAMPLE1_PSD_S_PERTURBATION_V11":"PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation":"SIGNED_CAYLEY_COMPOSE_SAMPLE1_S_AND_ACCELERATOR_INSIDE_Q8" if ok else "REFINE_DOMINANT_PSD_OR_S_PERTURBATION_TERM",
        "failures":failures,"rows":rows,
    }


def validate(d:dict)->list[str]:
    f=list(d.get("failures",[]))
    for k in ("source_generated_not_trajectory_fit","V10_canonical_core_retained","V10_exact_combined_perpendicular_residual_retained",
              "first_attitude_PSD_diagonal_remainder_already_in_V10","first_attitude_PSD_cross_axis_remainder_included","second_prediction_attitude_process_remainder_included",
              "sample1_S_covariance_update_included","sample1_S_attitude_injection_included","sample1_S_aw_mean_correction_included","sample1_S_solver_identity_branch_contained_as_zero_perturbation"):
        if d.get(k) is not True:f.append(f"{k} is not true")
    for k in ("source_replay_used","filter_changed","broad_sample1_3x3_interval_inverse_reintroduced","temporal_force_slew_assumed","complete_sample1_branch_closed_here",
              "signed_cayley_q8_composed_here","q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here"):
        if d.get(k) is not False:f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells",0))<=0:f.append("no V11 cells")
    for k in ("max_total_residual_perturbation_upper_mps2","max_total_reduced_covariance_perturbation_upper","max_sample1_gain_operator_perturbation_upper","max_V11_correction_norm_upper_rad"):
        if not math.isfinite(float(d.get(k,math.nan))):f.append(f"nonfinite {k}")
    st=d.get("P5_SAMPLE1_PSD_S_PERTURBATION_V11"); w=d.get("first_unclosed_joint_cell")
    if st=="PASS" and (w is not None or int(d.get("unclosed_joint_cells",-1))!=0):f.append("PASS retains V11 unclosed cell")
    if st=="NOT_ESTABLISHED" and w is None and not f:f.append("missing V11 witness")
    return f


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces",type=int,default=4); ap.add_argument("--source-cell-index",type=int,default=0)
    ap.add_argument("--p-pieces",type=int,default=24); ap.add_argument("--tangent-pieces",type=int,default=24); ap.add_argument("--axial-pieces",type=int,default=24)
    ap.add_argument("--output",type=Path,required=True); x=ap.parse_args()
    d=build(x.domain,source_pieces=x.source_pieces,source_cell_index=x.source_cell_index,p_pieces=x.p_pieces,tangent_pieces=x.tangent_pieces,axial_pieces=x.axial_pieces)
    vf=validate(d); d["validation_failures"]=vf; x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"status":d["P5_SAMPLE1_PSD_S_PERTURBATION_V11"],"cells":d["evaluated_joint_cells"],"unclosed":d["unclosed_joint_cells"],
        "eps":d["attitude_covariance_remainder_spectral_upper"],"S_bounds":d["sample1_S_perturbation_bounds"],"max_drho":d["max_total_residual_perturbation_upper_mps2"],
        "max_dP":d["max_total_reduced_covariance_perturbation_upper"],"max_dK":d["max_sample1_gain_operator_perturbation_upper"],
        "min_lambda":d["minimum_perturbed_innovation_lambda_lower"],"max_d":d["max_V11_correction_norm_upper_rad"],
        "first_unclosed":d["first_unclosed_joint_cell"],"worst":d["worst_joint_cell"],"next":d["next_obligation"],"validation_failures":vf},indent=2,sort_keys=True)); return 0 if not vf else 2

if __name__=="__main__": raise SystemExit(main())
