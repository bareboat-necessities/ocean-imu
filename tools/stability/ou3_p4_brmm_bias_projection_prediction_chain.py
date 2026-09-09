#!/usr/bin/env python3
"""Exact A21 projection -> next BIAS1 prediction bias-subgraph for P4.

The local A21 projection prefix already exposes on one common coordinate z

    f_ba = projected physical accelerometer-bias error,
    beta = persistent physical true bias.

The qualified BIAS1 prediction lift evolves

    [e_b+; beta+] = A_b [f_ba; beta] + B_b [w; m_tau].

This module composes those maps exactly on one extended coordinate

    Z = [z ; w(3) ; m_tau(3)],

and embeds the BIAS1 joint-supply ball IQC on the same homogeneous h coordinate.
It therefore proves the projected error is the *actual next prediction input*
rather than a detached graph output.  No independent beta slot is introduced.

Only the six-dimensional bias/true-bias subgraph is closed here.  The remaining
18 shipping error coordinates still require the full finite-angle event
transport before production prefix LDLT can be assembled.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_add,matrix_mul
import ou3_p4_brmm_a21_projection_prefix as PROJP
import ou3_p4_brmm_bias1_prediction_lift as BIASP
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_bias1_joint_iss_supply as BIASISS
import ou3_p4_same_cell_correction_domain as CORR

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_A21_PROJECTION_TO_BIAS1_PREDICTION_CHAIN_V1'

def I(x):return Interval.point(float(x))
def shape(A):return len(A),len(A[0]) if A else 0
def zeros(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def vcat(A,B):
    if shape(A)[1]!=shape(B)[1]:raise ValueError('vcat width mismatch')
    return [list(r) for r in A]+[list(r) for r in B]
def extend(A,n):
    r,c=shape(A)
    if c>n:raise ValueError('cannot shrink map')
    return [list(row)+[I(0) for _ in range(n-c)] for row in A]
def selector(rows,n,offset):
    return [[I(1 if j==offset+i else 0) for j in range(n)] for i in range(rows)]

def compose_projection_to_prediction(prefix,bias_contract):
    n0=int(prefix['coordinate_dimension']);s0=n0;n=n0+6
    F=extend(prefix['projected_ba_map'],n);Beta=extend(prefix['beta_map'],n)
    X=vcat(F,Beta)  # [f_ba; beta]
    A6,B6=BIASP.build_matrices(bias_contract)
    S=selector(6,n,s0)
    Y=matrix_add(matrix_mul(A6,X),matrix_mul(B6,S))
    if shape(Y)!=(6,n):raise RuntimeError('bias prediction output map shape drifted')
    # Same homogeneous h from the projection prefix; supply occupies final 6 slots.
    supply_bound=float(bias_contract['joint_supply_norm_upper_per_prediction_mps2'])
    supply_iqc=HARD.ball_iqc(n,int(prefix['h_index']),tuple(range(s0,s0+6)),supply_bound)
    return {'coordinate_dimension':n,'source_offset':s0,'input_bias_true_map':X,
            'next_bias_true_map':Y,'supply_map':S,'supply_iqc':supply_iqc,
            'h_index':int(prefix['h_index']),'supply_bound':supply_bound}

def build():
    corr=CORR.build();cf=CORR.validate(corr);bc=BIASISS.build();bf=BIASISS.validate(bc)
    if cf or bf:raise RuntimeError(f'chain prerequisites failed correction={cf} bias={bf}')
    cell=PROJP._smoke_cell();p=PROJP.build_augmented_prefix(cell,corr);c=compose_projection_to_prediction(p,bc)
    n=c['coordinate_dimension'];Y=c['next_bias_true_map'];S=c['supply_map']
    # Structural checks: same w columns feed e_b and beta, m only e_b.
    w_shared=all(Y[i][c['source_offset']+i].contains(1) and Y[3+i][c['source_offset']+i].contains(1) for i in range(3))
    m_error_only=all(Y[i][c['source_offset']+3+i].contains(1) and Y[3+i][c['source_offset']+3+i].contains(0) for i in range(3))
    beta_input_present=any(not (Y[i][j].lo==0 and Y[i][j].hi==0) for i in range(6) for j in range(p['beta_offset'],p['beta_offset']+3))
    projected_input_present=any(not (Y[i][j].lo==0 and Y[i][j].hi==0) for i in range(3) for j in range(p['projected_ba_offset'],p['projected_ba_offset']+3))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'A21_projection_prefix_consumed':True,'projected_ba_is_next_prediction_error_input':projected_input_present,
      'same_persistent_beta_is_next_prediction_true_bias_input':beta_input_present,
      'BIAS1_prediction_A_B_matrices_consumed':True,'same_w_enters_next_error_and_true_bias':w_shared,
      'm_tau_enters_next_error_not_true_bias':m_error_only,'joint_supply_ball_IQC_embedded_on_same_h':shape(c['supply_iqc'])==(n,n),
      'independent_beta_slot_created_for_next_prediction':False,'projection_output_charged_as_exogenous_supply':False,
      'packet_count_multiplier_used':False,'trajectory_replay_used':False,
      'chain_coordinate_dimension':n,'joint_supply_bound_mps2':c['supply_bound'],
      'production_other_18_state_coordinates_composed_here':False,'production_full_24state_prefix_transport_closed_here':False,
      'production_prefix_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'combine this exact bias subgraph with the finite-angle transport of the other 18 error coordinates and binary32 state forcing, preserving the same event/root coordinate, then form the cumulative prefix master'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('A21_projection_prefix_consumed','projected_ba_is_next_prediction_error_input','same_persistent_beta_is_next_prediction_true_bias_input','BIAS1_prediction_A_B_matrices_consumed','same_w_enters_next_error_and_true_bias','m_tau_enters_next_error_not_true_bias','joint_supply_ball_IQC_embedded_on_same_h'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_beta_slot_created_for_next_prediction','projection_output_charged_as_exogenous_supply','packet_count_multiplier_used','trajectory_replay_used','production_other_18_state_coordinates_composed_here','production_full_24state_prefix_transport_closed_here','production_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not(math.isfinite(float(d.get('joint_supply_bound_mps2',math.nan))) and float(d['joint_supply_bound_mps2'])>=0):f.append('supply bound invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'projected_to_prediction':d['projected_ba_is_next_prediction_error_input'],'same_beta':d['same_persistent_beta_is_next_prediction_true_bias_input'],'same_w':d['same_w_enters_next_error_and_true_bias'],'production':d['production_full_24state_prefix_transport_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
