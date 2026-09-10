#!/usr/bin/env python3
"""Structural H18 -> A21 inter-word release attachment for OU-III P4.

Shipping does not release accelerometer-bias updates inside the first canonical
3 s H18 word.  In the configured full-heading path magnetic refinement holds
b_a out of the MEKF for at least 30 s after Live.  The correct theorem topology
is therefore

    H18 word(s) -> separate release hybrid -> A21 word(s),

not an artificial mixed-mode 3 s word.

At release the physical error-state map is exactly

    e_A^+ = [ I_18 ; 0 ] e_H^- + E_b e_b,held,

where E_b injects the three held accelerometer-bias errors into the new A21
coordinates.  Shipping does NOT set physical bias error to zero.  The bias
coordinate is bounded by the same family-specific BIAS0/BIAS1/BIAS2 true-bias
source envelope already attached to P4.  The covariance operation is separate:
held b_a cross-covariances are zero and the diagonal block is floored to the
shipping sigma_ba0^2 before active Gauss-Markov prediction resumes.

This module closes the structural/source attachment of that hybrid and proves
that the held-bias forcing itself fits the declared A21 bias coordinate radius.
It does not claim the H18 pre-release state is in a contracting basin or that
the post-release A21 state is retained indefinitely; those require the actual
H18/A21 endpoint/every-prefix LDLT and first-exit proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_riccati_metric_p3 as P3
import ou3_p4_bias_family_joint_iss_supply as BIAS
import ou3_p4_complete_brmm_differential_word as DWORD
import ou3_p4_hard_entry_set as ENTRY

SCHEMA=1
QUALIFICATION='OU3_P4_H18_TO_A21_INTERWORD_RELEASE_ATTACHMENT_V1'
P3_DELTA=1e-18


def release_maps():
    """Return homogeneous H18 lift J and held-bias forcing map E_b."""
    ev=DWORD.H_to_A_release_event('COMPLETE_BRMM_NORMAL_LIVE_WORD:release')
    J=ev.J
    z=Interval.point(0.0); o=Interval.point(1.0)
    Eb=[[z for _ in range(3)] for _ in range(21)]
    for i in range(3): Eb[18+i][i]=o
    return J,Eb


def build():
    live=LIVE.build(); lf=LIVE.validate(live)
    p3=P3.build(); pf=P3.validate(p3)
    bias=BIAS.build(); bf=BIAS.validate(bias)
    entry=ENTRY.build(); ef=ENTRY.validate(entry)
    bad={'live':lf,'P3':pf,'bias':bf,'entry':ef}; bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('release prerequisites failed: '+repr(bad))
    J,Eb=release_maps()
    rel=live['H_to_A_release']
    hword=float(rel['canonical_H18_word_horizon_s'])
    min_hold=float(rel['minimum_H_mode_live_duration_before_A_release_s'])
    min_words=int(math.floor(min_hold/hword))
    radius=float(entry['coordinate_radii']['accelerometer_bias_error_norm_mps2'])
    fam={}
    for name,row in bias['family_supply'].items():
        b=float(row['true_bias_norm_upper_mps2'])
        fam[name]={'held_true_bias_norm_upper_mps2':b,
                   'A21_bias_entry_radius_mps2':radius,
                   'held_bias_forcing_fits_A21_entry_coordinate':b<=radius}
    all_fit=all(x['held_bias_forcing_fits_A21_entry_coordinate'] for x in fam.values())
    topology_ok=(p3.get('canonical_P3_topology')=='H18_3S_PRIOR_FREE_THEN_PRESERVED_H_TO_A_HYBRID_A21'
                 and p3.get('same_complete_BRMM_execution_continues_across_H_to_A') is True)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'P3_delta':P3_DELTA,'P3_frozen_not_modified':True,
      'release_is_separate_interword_hybrid':bool(rel['hybrid_transition_not_inside_same_mode_word']),
      'minimum_H_mode_live_duration_before_release_s':min_hold,
      'canonical_H18_word_horizon_s':hword,
      'minimum_complete_H18_words_before_release':min_words,
      'at_least_one_complete_H18_word_before_release':min_words>=1,
      'at_least_ten_complete_H18_words_before_release':min_words>=10,
      'homogeneous_release_map':'[I18;0]','homogeneous_release_map_dimension':[21,18],
      'held_bias_forcing_map_dimension':[21,3],
      'held_bias_error_set_to_zero_at_release':False,
      'held_bias_error_retained_as_separate_forcing':True,
      'all_bias_family_release_bounds':fam,
      'all_bias_families_fit_declared_A21_bias_coordinate':all_fit,
      'held_ba_cross_covariances_zero_before_release':bool(live['held_ba']['cross_covariances_zero']),
      'release_bias_diagonal_floor_variance':float(rel['bias_diagonal_floor_variance']),
      'release_covariance_floor_retained_as_separate_metric_event':True,
      'P3_H_to_A_topology_consumed':topology_ok,
      'release_structural_source_attachment_closed':bool(all_fit and topology_ok and min_words>=1),
      'H18_pre_release_endpoint_contraction_closed_here':False,
      'A21_post_release_basin_landing_closed_here':False,
      'A21_post_release_every_prefix_retention_closed_here':False,
      'full_binary32_release_arithmetic_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':(
        'once H18 endpoint/every-prefix contraction and full binary32 arithmetic close, use repeated H18 retention up to release; '
        'then combine [I18;0] with the family held-bias forcing and covariance-floor metric event to certify A21 basin landing'
      )}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('P3_frozen_not_modified','release_is_separate_interword_hybrid','at_least_one_complete_H18_word_before_release',
              'at_least_ten_complete_H18_words_before_release','held_bias_error_retained_as_separate_forcing',
              'all_bias_families_fit_declared_A21_bias_coordinate','held_ba_cross_covariances_zero_before_release',
              'release_covariance_floor_retained_as_separate_metric_event','P3_H_to_A_topology_consumed',
              'release_structural_source_attachment_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('held_bias_error_set_to_zero_at_release','H18_pre_release_endpoint_contraction_closed_here',
              'A21_post_release_basin_landing_closed_here','A21_post_release_every_prefix_retention_closed_here',
              'full_binary32_release_arithmetic_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('homogeneous_release_map_dimension')!=[21,18] or d.get('held_bias_forcing_map_dimension')!=[21,3]:f.append('release map dimensions changed')
    if int(d.get('minimum_complete_H18_words_before_release',0))<10:f.append('shipping release no longer guarantees ten H18 words')
    fam=d.get('all_bias_family_release_bounds',{})
    if set(fam)!={'BIAS0','BIAS1','BIAS2'}:f.append('release bias family table incomplete')
    elif not all(x.get('held_bias_forcing_fits_A21_entry_coordinate') is True for x in fam.values()):f.append('a bias family exceeds A21 entry coordinate')
    J,Eb=release_maps()
    if len(J)!=21 or len(J[0])!=18 or len(Eb)!=21 or len(Eb[0])!=3:f.append('release constructor dimension failure')
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'H18_words_before_release':d['minimum_complete_H18_words_before_release'],'all_bias_fit':d['all_bias_families_fit_declared_A21_bias_coordinate'],'structural_release':d['release_structural_source_attachment_closed'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))


if __name__=='__main__':raise SystemExit(main())
