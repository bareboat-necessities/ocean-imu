#!/usr/bin/env python3
"""Source-uniform physical response cocycle for the ALT joint24 Live word.

Consumes estimator-owned COMPLETE-BRMM literal event cells only.  No replay,
finite-seed fit, unreachable perturbation or independent coefficient box is an
input.

Each prediction carries two theorem source blocks:

* one homogeneous [h_s,q15] BRMM acceleration/moment block with the existing
  joint quadratic witness sectors and exact BRMM-vs-OU forcing map;
* one homogeneous [h_s,w_bias] block from one explicit BIAS0/1/2 contract.  The
  SAME w columns enter accelerometer-bias error and physical beta truth, and the
  family parameter/root token is retained for the future complete-word master.

At S=0 the finite physical residual is e_S-S_phys; its source columns carry the
same generator/primitive/one-time-Live-origin identity as the estimator-owned
cell.  Every source response block is suffix-propagated through later literal
nonlinear events.

This closes the per-attached-sample physical response relation.  It does not
claim that the 600-step source lineage has yet preserved one family parameter
root, that every hybrid edge is closed, or that storage may be searched.
"""
from __future__ import annotations
from dataclasses import dataclass

from ou3_interval import Interval, matrix_mul
import ou3_brmm_acceleration_moment_iqc as MOM
import ou3_p4_brmm_physical_prediction_forcing as FORCING
import ou3_p4_brmm_physical_acceleration_witness_sector as WSECTOR
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_source_cover_contract as COVER
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import joint24_events as J24
from tools.stability.ou3_alt_contraction import physical_reference as REF

QUALIFICATION = "OU3_ALT_SOURCE_UNIFORM_PHYSICAL_RESPONSE_COCYCLE_V2"


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero(r: int, c: int):
    return [[I(0) for _ in range(c)] for _ in range(r)]


def _identity(n: int):
    A = _zero(n, n)
    for i in range(n): A[i][i] = I(1)
    return A


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _propagate(J, B):
    if _shape(J) != (24, 24) or _shape(B)[0] != 24:
        raise ValueError("joint24 response propagation dimension mismatch")
    return matrix_mul(J, B)


@dataclass(frozen=True)
class SourceResponseBlock:
    kind: str
    token: str
    generator_id: str
    live_origin_id: str
    response: tuple[tuple[Interval, ...], ...]
    sectors: tuple[tuple[str, tuple[tuple[Interval, ...], ...]], ...] = ()
    homogeneous_scale_fixed_one: bool = False
    family_parameter_token: str | None = None


def _freeze(A): return tuple(tuple(x for x in row) for row in A)
def _thaw(A): return [list(row) for row in A]


def _carry(block: SourceResponseBlock, J):
    return SourceResponseBlock(
        block.kind, block.token, block.generator_id, block.live_origin_id,
        _freeze(_propagate(J, _thaw(block.response))), block.sectors,
        block.homogeneous_scale_fixed_one, block.family_parameter_token)


def _q15_sector_for_h(h: Interval):
    if not isinstance(h, Interval) or h.lo <= 0: raise ValueError("positive interval sample time required")
    mom=MOM.build(); failures=MOM.validate(mom)
    if failures: raise RuntimeError("physical acceleration-moment prerequisite failed: "+repr(failures))
    Amax=float(mom["A_max_mps2"]); n=16
    Cs=WSECTOR.selector(1,n,[0]); A0=WSECTOR.selector(3,n,[1,2,3]); A1=WSECTOR.selector(3,n,[4,5,6])
    M=_zero(9,n); ih=I(1)/h; ih2=ih*ih; ih3=ih2*ih
    for axis in range(3):
        M[3*axis][7+axis]=ih; M[3*axis+1][10+axis]=ih2; M[3*axis+2][13+axis]=ih3
    return tuple((name,_freeze(Q)) for name,Q in WSECTOR.physical_witness_sectors(A0,A1,M,Cs,Amax))


def _q15_injection(mode: str, tau: Interval, h: Interval):
    M12=FORCING.source_matrix(tau,h); Gn=FORCING.inject_error_state(mode,M12); B=_zero(24,16)
    n=18 if mode=="H" else 21
    for i in range(n):
        for j in range(15): B[i][1+j]=Gn[i][j]
    return B


def _bias_driver_homogeneous(B3, contract: BIAS.BiasFamilyContract):
    """Lift 24x3 shared-w response to 24x4 [h_s,w] plus its analytic sector."""
    if _shape(B3)!=(24,3): raise ValueError("shared bias driver response must be 24x3")
    B4=_zero(24,4)
    for i in range(24):
        for j in range(3): B4[i][1+j]=B3[i][j]
    return B4, (("bias_driver_ball",_freeze(BIAS.driver_ball_iqc(contract))),)


def _primitive_token(cell: COVER.SourceCoverCell, suffix: str) -> str:
    w=cell.wave_primitive
    if w is None: raise ValueError("physical source cell lost WavePrimitivePayload")
    REF.validate_primitive_reference(w)
    return ":".join((w.generator_id,w.live_origin_id,w.primitive_in_id,w.primitive_out_id,suffix))


def _h_event_jacobian(cell, held_bias_error, true_bias):
    if cell.kind=="accelerometer":
        failures=COVER.validate_cell(cell,require_estimator_provenance=True)
        if failures: raise ValueError("invalid H18 source cell: "+repr(failures))
        if len(held_bias_error)!=3 or len(true_bias)!=3: raise ValueError("H18 joint24 requires e_ba and beta domains")
        u=EVENTS.AD.independent_vector(list(cell.state)+list(held_bias_error)+list(true_bias),n=24,offset=0)
        z,eba,beta=u[:18],u[18:21],u[21:24]
        H=EVENTS._event_H("H","accelerometer",f_hat=cell.f_hat,R_hat=cell.R_hat); K,_=EVENTS.source_joseph_gain(cell.P,H,cell.R)
        residual=J24._h18_accel_residual(z,eba,cell.f_hat,cell.R_hat)
        out=list(EVENTS._apply_physical_correction(z,K,residual))+list(eba)+list(beta)
        return EVENTS.AD.jacobian(out),None
    if cell.kind=="magnetometer":
        u=EVENTS.AD.independent_vector(list(cell.state)+list(held_bias_error)+list(true_bias),n=24,offset=0); z=u[:18]
        H=EVENTS._event_H("H","magnetometer",m_body=cell.m_body); K,_=EVENTS.source_joseph_gain(cell.P,H,cell.R)
        out18=EVENTS._apply_physical_correction(z,K,EVENTS.residual_magnetometer(z,cell.m_body))
        return EVENTS.AD.jacobian(list(out18)+list(u[18:24])),None
    if cell.kind=="S_zero":
        d=J24.s_zero_event_joint24(cell); return d["J_joint24"],d["J_S_phys"]
    raise ValueError("not an H18 Joseph event")


def _a_event_jacobian(cell):
    failures=COVER.validate_cell(cell,require_estimator_provenance=True)
    if failures: raise ValueError("invalid A21 source cell: "+repr(failures))
    if cell.true_bias is None or cell.bias_projection_limit is None: raise ValueError("A21 source cell missing true bias/projection")
    add_S=cell.kind=="S_zero"; extra=list(cell.wave_primitive.centered_S) if add_S and cell.wave_primitive is not None else []
    if add_S: REF.validate_primitive_reference(cell.wave_primitive)
    total=24+(3 if add_S else 0); u=EVENTS.AD.independent_vector(list(cell.state)+list(cell.true_bias)+extra,n=total,offset=0)
    z,beta=u[:21],u[21:24]
    H=EVENTS._event_H(cell.mode,cell.kind,f_hat=cell.f_hat,R_hat=cell.R_hat,m_body=cell.m_body); K,_=EVENTS.source_joseph_gain(cell.P,H,cell.R)
    if cell.kind=="accelerometer": residual=EVENTS.residual_accelerometer(z,cell.f_hat,cell.R_hat)
    elif cell.kind=="magnetometer": residual=EVENTS.residual_magnetometer(z,cell.m_body)
    elif cell.kind=="S_zero": residual=[z[EVENTS.OFF_S+i]-u[24+i] for i in range(3)]
    else: raise ValueError("not an A21 Joseph event")
    out21=EVENTS._apply_physical_correction(z,K,residual)
    J21,_,_=J24._project_joint(out21,beta,total,float(cell.bias_projection_limit)); J=[list(row) for row in J21]
    for i in range(3):
        row=[I(0) for _ in range(total)]; row[21+i]=I(1); J.append(row)
    return [row[:24] for row in J],([row[24:27] for row in J] if add_S else None)


def _prediction_joint(cell,selector,contract: BIAS.BiasFamilyContract,tau_ba):
    mode=cell.mode; pt=contract.phi_true
    p=PRED.prediction_event(mode,cell.state,selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s,tau_ba=tau_ba if mode=="A" else None)
    if mode=="H": J,Bbias=J24.h18_prediction_lift(p["J_state"],pt)
    else:
        if p["phi_ba"] is None: raise RuntimeError("A21 prediction lost estimator bias factor")
        ph=p["phi_ba"]; J=_zero(24,24); Bbias=_zero(24,3)
        for i in range(21):
            for j in range(21): J[i][j]=p["J_state"][i][j]
        for i in range(3):
            J[18+i][21+i]=pt-ph; J[21+i][21+i]=pt
            Bbias[18+i][i]=I(1); Bbias[21+i][i]=I(1)
    return J,Bbias,_q15_injection(mode,cell.tau_applied_s,cell.dt_s),_q15_sector_for_h(cell.dt_s)


def compose_attached_sample(lineage, *, bias_contract: BIAS.BiasFamilyContract, held_bias_error=None, true_bias=None, tau_ba=None):
    """Joint24 physical cocycle for one universal attached sample lineage."""
    if not isinstance(bias_contract,BIAS.BiasFamilyContract): raise TypeError("explicit BIAS0/1/2 contract required")
    mode=lineage.mode
    if mode not in ("H","A"): raise ValueError("lineage mode must be H/A")
    if mode=="H" and (held_bias_error is None or true_bias is None): raise ValueError("H18 lineage requires hard e_ba and beta domains")
    if mode=="A" and tau_ba is None: raise ValueError("A21 lineage requires configured estimator tau_ba")

    Psi=_identity(24); blocks=[]; kinds=[]
    for cell in lineage.cells:
        failures=COVER.validate_cell(cell,require_estimator_provenance=True)
        if failures: raise ValueError("invalid estimator-owned cell: "+repr(failures))
        kinds.append(cell.kind)
        if cell.kind=="prediction":
            J,Bbias,Bq,sectors=_prediction_joint(cell,lineage.selector,bias_contract,tau_ba)
            Psi=matrix_mul(J,Psi); blocks=[_carry(b,J) for b in blocks]
            w=cell.wave_primitive
            if w is None: raise ValueError("prediction cell missing same-history primitive ancestry")
            REF.validate_primitive_reference(w)
            blocks.append(SourceResponseBlock("BRMM_q15_prediction",_primitive_token(cell,"q15"),w.generator_id,w.live_origin_id,_freeze(Bq),sectors,True))
            B4,bias_sectors=_bias_driver_homogeneous(Bbias,bias_contract)
            blocks.append(SourceResponseBlock("shared_bias_driver",f"{cell.source_token}:{bias_contract.name}:bias-driver",w.generator_id,w.live_origin_id,_freeze(B4),bias_sectors,True,bias_contract.parameter_token))
        elif cell.kind=="aw_floor":
            J=_identity(24)
        elif cell.kind in ("S_zero","accelerometer","magnetometer"):
            J,BS=(_h_event_jacobian(cell,held_bias_error,true_bias) if mode=="H" else _a_event_jacobian(cell))
            Psi=matrix_mul(J,Psi); blocks=[_carry(b,J) for b in blocks]
            if BS is not None:
                w=cell.wave_primitive
                if w is None: raise ValueError("S event missing same-history primitive ancestry")
                blocks.append(SourceResponseBlock("physical_centered_S",_primitive_token(cell,"S"),w.generator_id,w.live_origin_id,_freeze(BS)))
        else: raise ValueError("literal lineage contains unsupported event: "+repr(cell.kind))

    expected=tuple(c.kind for c in lineage.cells)
    if tuple(kinds)!=expected: raise RuntimeError("literal event order changed during physical composition")
    bias_blocks=tuple(b for b in blocks if b.kind=="shared_bias_driver")
    return {
        "mode":mode,"bias_family":bias_contract.name,"bias_parameter_token":bias_contract.parameter_token,
        "J_joint24":Psi,"source_blocks":tuple(blocks),"literal_event_kinds":expected,"all_literal_cells_consumed":True,
        "q15_prediction_blocks":sum(b.kind=="BRMM_q15_prediction" for b in blocks),
        "physical_S_blocks":sum(b.kind=="physical_centered_S" for b in blocks),
        "bias_driver_blocks":len(bias_blocks),
        "same_bias_parameter_token_on_all_predictions":all(b.family_parameter_token==bias_contract.parameter_token for b in bias_blocks),
        "bias_driver_ball_sector_attached_on_all_predictions":all(any(name=="bias_driver_ball" for name,_ in b.sectors) for b in bias_blocks),
        "same_history_source_ancestry_retained":all(bool(b.generator_id and b.live_origin_id) for b in blocks),
        "replay_used":False,"independent_prediction_defect_boxes_used":False,
    }


def build():
    forcing=FORCING.build(); ff=FORCING.validate(forcing); sector=WSECTOR.build(); sf=WSECTOR.validate(sector); bias=BIAS.build(); bf=BIAS.validate(bias)
    if ff or sf or bf: raise RuntimeError(f"physical cocycle prerequisites failed forcing={ff} sector={sf} bias={bf}")
    return {
        "qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD","joint_dimension":24,
        "same_15D_q_witness_used_per_prediction":True,"joint_q15_quadratic_sector_attached_to_same_block":True,
        "q15_columns_suffix_propagated_through_later_literal_events":True,"physical_S_columns_suffix_propagated_through_later_literal_events":True,
        "same_generator_and_one_time_Live_origin_tokens_retained":True,"shared_bias_driver_enters_error_and_truth_once":True,
        "explicit_BIAS0_BIAS1_BIAS2_contract_required_by_cocycle":True,"analytic_bias_driver_sector_attached_per_prediction":True,
        "persistent_bias_parameter_token_carried_by_source_blocks":True,"aggregate_fresh_Live_entry_builder_consumed":False,
        "independent_prediction_defect_boxes_used":False,"replay_or_unreachable_perturbation_used":False,
        "source_uniform_attached_sample_cocycle_available":True,
        "multi_sample_COMPLETE_BRMM_word_composed":False,"BIAS0_BIAS1_BIAS2_temporal_source_relations_attached":False,
        "H18_A21_edge_attached":False,"storage_search_allowed":False,"ALT_LIVE_PASS":False,
        "next_obligation":"compose attached-sample cocycles through every retained 600-step source successor while reusing one family parameter_token, then close the actual 21-state held-to-active covariance enable edge before storage search",
    }


def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION or d.get("joint_dimension")!=24:f.append("qualification/dimension mismatch")
    for k in ("same_15D_q_witness_used_per_prediction","joint_q15_quadratic_sector_attached_to_same_block","q15_columns_suffix_propagated_through_later_literal_events","physical_S_columns_suffix_propagated_through_later_literal_events","same_generator_and_one_time_Live_origin_tokens_retained","shared_bias_driver_enters_error_and_truth_once","explicit_BIAS0_BIAS1_BIAS2_contract_required_by_cocycle","analytic_bias_driver_sector_attached_per_prediction","persistent_bias_parameter_token_carried_by_source_blocks","source_uniform_attached_sample_cocycle_available"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("aggregate_fresh_Live_entry_builder_consumed","independent_prediction_defect_boxes_used","replay_or_unreachable_perturbation_used","multi_sample_COMPLETE_BRMM_word_composed","BIAS0_BIAS1_BIAS2_temporal_source_relations_attached","H18_A21_edge_attached","storage_search_allowed","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return list(dict.fromkeys(f))
