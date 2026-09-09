#!/usr/bin/env python3
"""Estimator-owned BRMM event-prefix -> augmented every-prefix LDLT.

Each literal prefix has two linked ancestries:

* ``image``: the same-signal estimator image that generated estimated
  (f,tau,sigma,T_S,R_S), candidate/staging and the current committed schedule;
* ``source_cell``: one literal Kalman/Joseph/reset event token owned by that
  image, with reachable state/P/geometry/R/bias/radial coordinates.

Several consecutive event prefixes may therefore share one estimator image, but
their event tokens must form a literal predecessor chain.  This is required for
endpoint AND every-prefix augmented LDLT.  The existing full interval ISS LDLT
backend is reused; no endpoint substitution or independent coefficient box is
accepted.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_joint_iss_augmented_master as ISS

SCHEMA=2
QUALIFICATION='OU3_P4_JOINT_BRMM_EVENT_PREFIX_AUGMENTED_LDLT_V2'

@dataclass(frozen=True)
class PrefixInput:
    image:JOINT.Image
    source_cell:COVER.SourceCoverCell
    master:Sequence[Sequence[Interval]]
    sectors:Sequence[Sequence[Sequence[Interval]]]
    multipliers:Sequence[float]
    source_map:Sequence[Sequence[Interval]]|None=None
    gamma_s:float=0.0
    fp_map:Sequence[Sequence[Interval]]|None=None
    gamma_n:float=0.0

def _shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):raise ValueError('ragged matrix')
    return r,c

def _same_interval(a,b):return isinstance(a,Interval) and isinstance(b,Interval) and a.lo==b.lo and a.hi==b.hi

def validate_sequence(cells):
    f=[]
    if not cells:return ['nonempty literal prefix sequence required']
    seen=set()
    for i,c in enumerate(cells):
        im=c.image;sc=c.source_cell
        vf=COVER.validate_cell(sc,require_estimator_provenance=True)
        if vf:f.extend(f'prefix {i}: '+x for x in vf)
        if not sc.source_token or sc.source_token in seen:f.append(f'prefix {i}: duplicate/missing event token')
        seen.add(sc.source_token)
        if i>0 and sc.predecessor_token!=cells[i-1].source_cell.source_token:f.append(f'prefix {i}: not literal previous-event ancestry')
        if sc.estimator_source_token!=im.source_token:f.append(f'prefix {i}: event coefficients detached from estimator image')
        if sc.estimator_predecessor_token!=im.predecessor_token:f.append(f'prefix {i}: estimator predecessor detached')
        if im.state.source_token!=im.source_token:f.append(f'prefix {i}: estimator state/image token detached')
        for name,x in (('f',im.frequency_hz),('tau',im.tau_target),('sigma',im.sigma_target_raw),('TS',im.pseudo_period_target),('RS',im.rs_target)):
            if not isinstance(x,Interval) or x.lo<=0:f.append(f'prefix {i}: invalid joint {name}')
        active=im.active_schedule_for_current_riccati
        if active.tau.lo<=0 or active.sigma.lo<=0 or active.rs_base.lo<=0:f.append(f'prefix {i}: invalid current active schedule')
        if not _same_interval(sc.tau_applied_s,active.tau):f.append(f'prefix {i}: cell tau detached from current active schedule')
        if not _same_interval(sc.sigma_aw_mps2,active.sigma):f.append(f'prefix {i}: cell sigma detached from current active schedule')
        if any(r.lo<=0 for r in im.actual_rs_std_xyz_for_current_riccati):f.append(f'prefix {i}: invalid actual per-axis RS')
        n,m=_shape(c.master)
        if n==0 or n!=m:f.append(f'prefix {i}: master must be nonempty square')
        if not c.sectors or len(c.sectors)!=len(c.multipliers):f.append(f'prefix {i}: sectors/multipliers missing')
        if c.source_map is None:f.append(f'prefix {i}: BIAS1/source map missing')
        if c.fp_map is None:f.append(f'prefix {i}: binary32 finite-precision map missing')
    return list(dict.fromkeys(f))

def certify(cells):
    f=validate_sequence(cells)
    if f:return {'closed':False,'validation_failures':f,'endpoint_closed':False,'every_prefix_closed':False,'prefixes':[]}
    rec=[]
    for i,c in enumerate(cells):
        ok,piv=ISS.certify_iss(c.master,c.sectors,c.multipliers,c.source_map,c.gamma_s,c.fp_map,c.gamma_n)
        rec.append({'prefix_length':i+1,'event_token':c.source_cell.source_token,'event_predecessor_token':c.source_cell.predecessor_token,'estimator_token':c.image.source_token,'ldlt_closed':bool(ok),'pivot_lowers':piv})
    every=all(x['ldlt_closed'] for x in rec);endpoint=bool(rec[-1]['ldlt_closed'])
    return {'closed':bool(every and endpoint),'validation_failures':[],'endpoint_closed':endpoint,'every_prefix_closed':every,'prefixes':rec}

def build():
    j=JOINT.build();jf=JOINT.validate(j);c=COVER.build();cf=COVER.validate(c)
    bad={k:v for k,v in (('joint',jf),('cover',cf)) if v}
    if bad:raise RuntimeError('BRMM event-prefix prerequisite failed: '+repr(bad))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'BRMM_attached_joint_image_required_at_every_prefix':True,
      'estimator_owned_source_cell_required_at_every_prefix':True,
      'literal_previous_prefix_ancestry_required':True,
      'multiple_events_may_share_one_estimator_image':True,
      'current_applied_schedule_retained_at_every_prefix':True,
      'raw_joint_f_tau_sigma_TS_RS_retained_at_every_prefix':True,
      'same_cell_augmented_master_sector_source_fp_maps_required':True,
      'full_interval_ISS_LDLT_backend_reused':True,'endpoint_cannot_substitute_for_every_prefix':True,
      'BIAS1_source_map_required_for_production':True,'binary32_fp_map_required_for_production':True,
      'production_complete_source_lineages_supplied_here':False,'production_endpoint_closed_here':False,'production_every_prefix_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'generate one PrefixInput per literal event for every retained BRMM/radial lineage; each source_cell must be built from its estimator image and carry reachable P/H/R/K/bias/projection plus nonempty BIAS1 and binary32 ISS maps, then require certify(...).closed for every lineage'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('BRMM_attached_joint_image_required_at_every_prefix','estimator_owned_source_cell_required_at_every_prefix','literal_previous_prefix_ancestry_required','multiple_events_may_share_one_estimator_image','current_applied_schedule_retained_at_every_prefix','raw_joint_f_tau_sigma_TS_RS_retained_at_every_prefix','same_cell_augmented_master_sector_source_fp_maps_required','full_interval_ISS_LDLT_backend_reused','endpoint_cannot_substitute_for_every_prefix','BIAS1_source_map_required_for_production','binary32_fp_map_required_for_production'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_complete_source_lineages_supplied_here','production_endpoint_closed_here','production_every_prefix_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'event_prefix_bridge':d['estimator_owned_source_cell_required_at_every_prefix'],'production_prefix':d['production_every_prefix_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())