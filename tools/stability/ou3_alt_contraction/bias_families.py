#!/usr/bin/env python3
"""Theorem-facing BIAS0/1/2 contracts for the ALT joint24 word.

This module consumes the three existing conditional physical source families
separately.  It never calls the aggregate fresh-Live entry/supply builder.
For each family it retains one physical bias history over the word,

    beta+ = phi_true beta + w,

with the family-specific hard interval for phi_true and analytic hard bound on
the SAME three-component driver w.  The same w columns are used by beta and the
accelerometer-bias error recurrence in ``physical_lineage``.

The homogeneous driver sector is

    ||w||^2 <= W^2 h_s^2,   h_s=1,

and is an analytic theorem constraint, not a replay fit.  The family contract
also carries the hard true-bias norm bound needed by projection/basin closure.
A persistent ``parameter_token`` is deliberately part of the contract: the
multi-sample master must reuse it at every prediction instead of choosing a new
phi/root at each sample.  This file makes those source contracts available; it
does not by itself prove that the 600-step master preserved the token.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

from ou3_interval import Interval
import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2

QUALIFICATION = "OU3_ALT_ANALYTIC_BIAS_FAMILY_CONTRACTS_V1"


def I(x: float) -> Interval:
    return Interval.point(float(x))


@dataclass(frozen=True)
class BiasFamilyContract:
    name: str
    qualification: str
    phi_true: Interval
    driver_component_bound: float
    driver_norm_bound: float
    true_bias_component_bound: float
    true_bias_norm_bound: float
    parameter_token: str
    one_history_required: bool
    non_relaxing_limit_admitted: bool


def driver_ball_iqc(contract: BiasFamilyContract):
    """Q on [h_s,w_x,w_y,w_z] with z^T Q z >= 0."""
    W=float(contract.driver_norm_bound)
    if not (math.isfinite(W) and W>=0):
        raise ValueError("finite nonnegative driver norm bound required")
    Q=[[I(0) for _ in range(4)] for _ in range(4)]
    Q[0][0]=I(W*W)
    for i in range(3): Q[1+i][1+i]=I(-1)
    return Q


def true_bias_ball_iqc(contract: BiasFamilyContract):
    """Q on [h_s,beta_x,beta_y,beta_z] with ||beta|| <= B."""
    B=float(contract.true_bias_norm_bound)
    if not (math.isfinite(B) and B>=0):
        raise ValueError("finite nonnegative true-bias norm bound required")
    Q=[[I(0) for _ in range(4)] for _ in range(4)]
    Q[0][0]=I(B*B)
    for i in range(3): Q[1+i][1+i]=I(-1)
    return Q


def _contract(name, module):
    d=module.build(); failures=module.validate(d)
    if failures:
        raise RuntimeError(f"{name} source family invalid: {failures!r}")
    lo,hi=map(float,d['phi_true_interval'])
    one_history = bool(
        d.get('one_composite_history_required')
        or d.get('one_root_one_parameter_history_required')
        or d.get('one_drift_history_required')
    )
    if not one_history:
        raise RuntimeError(name+" lost one-history source requirement")
    c=BiasFamilyContract(
        name=name, qualification=str(d['qualification']),
        phi_true=Interval.outward_bounds(lo,hi),
        driver_component_bound=float(d['driver_increment_component_abs_upper_mps2']),
        driver_norm_bound=float(d['driver_increment_norm_upper_mps2']),
        true_bias_component_bound=float(d['true_bias_component_abs_upper_mps2']),
        true_bias_norm_bound=float(d['true_bias_norm_upper_mps2']),
        parameter_token=f"{name}:one-physical-history:phi-root-and-driver",
        one_history_required=True,
        non_relaxing_limit_admitted=bool(d.get('non_relaxing_limit_admitted',False)),
    )
    if not (0<c.phi_true.lo<=c.phi_true.hi<=1):
        raise RuntimeError(name+" invalid phi interval")
    return c


def contracts():
    return tuple(_contract(name,module) for name,module in (
        ('BIAS0',BIAS0),('BIAS1',BIAS1),('BIAS2',BIAS2)))


def build():
    cs=contracts()
    return {
        'qualification':QUALIFICATION,
        'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
        'families':tuple(c.name for c in cs),
        'three_families_invoked_separately':tuple(c.name for c in cs)==('BIAS0','BIAS1','BIAS2'),
        'one_persistent_parameter_token_per_family':all(c.one_history_required and bool(c.parameter_token) for c in cs),
        'analytic_driver_ball_IQC_available_for_each_family':all(driver_ball_iqc(c) for c in cs),
        'analytic_true_bias_ball_IQC_available_for_each_family':all(true_bias_ball_iqc(c) for c in cs),
        'BIAS2_non_relaxing_limit_retained':next(c for c in cs if c.name=='BIAS2').non_relaxing_limit_admitted,
        'aggregate_fresh_Live_entry_builder_consumed':False,
        'replay_or_finite_seed_bound_used':False,
        'independent_per_sample_beta_slots_used':False,
        'multi_sample_parameter_token_continuity_closed_here':False,
        'all_bias_families_attached_to_complete_word':False,
        'ALT_LIVE_PASS':False,
        'next_obligation':'use one selected family contract on every prediction of a complete source lineage, preserve its parameter_token and shared w response columns across all samples, and close all three family masters separately',
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if tuple(d.get('families',()))!=('BIAS0','BIAS1','BIAS2'):f.append('family set/order changed')
    for k in ('three_families_invoked_separately','one_persistent_parameter_token_per_family','analytic_driver_ball_IQC_available_for_each_family','analytic_true_bias_ball_IQC_available_for_each_family','BIAS2_non_relaxing_limit_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('aggregate_fresh_Live_entry_builder_consumed','replay_or_finite_seed_bound_used','independent_per_sample_beta_slots_used','multi_sample_parameter_token_continuity_closed_here','all_bias_families_attached_to_complete_word','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
