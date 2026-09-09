#!/usr/bin/env python3
"""BRMM-attached same-history joint target -> every-prefix augmented LDLT.

Unlike the older structural estimator bridge, each prefix here carries an image
emitted by ``ou3_p4_joint_brmm_frontend_transition``: private Mahony, same-signal
stillness, tuner band/moments, previous-WPE f, raw target, current active
schedule, and WPE successor all share one ancestry token.  The existing full
interval ISS LDLT remains the backend.  Endpoint-only success is insufficient.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_joint_iss_augmented_master as ISS

SCHEMA=1
QUALIFICATION='OU3_P4_JOINT_BRMM_EVERY_PREFIX_AUGMENTED_LDLT_V1'

@dataclass(frozen=True)
class PrefixInput:
    image:JOINT.Image
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

def validate_sequence(cells):
    f=[]
    if not cells:return ['nonempty literal prefix sequence required']
    seen=set()
    for i,c in enumerate(cells):
        im=c.image
        if not im.source_token or im.source_token in seen:f.append(f'prefix {i}: duplicate/missing source token')
        seen.add(im.source_token)
        if not im.predecessor_token:f.append(f'prefix {i}: missing predecessor')
        if i>0 and im.predecessor_token!=cells[i-1].image.source_token:f.append(f'prefix {i}: not literal previous-prefix ancestry')
        if im.state.source_token!=im.source_token:f.append(f'prefix {i}: state/image token detached')
        for name,x in (('f',im.frequency_hz),('tau',im.tau_target),('sigma',im.sigma_target_raw),('TS',im.pseudo_period_target),('RS',im.rs_target)):
            if not isinstance(x,Interval) or x.lo<=0:f.append(f'prefix {i}: invalid joint {name}')
        if im.active_schedule_for_current_riccati.tau.lo<=0 or im.active_schedule_for_current_riccati.sigma.lo<=0 or im.active_schedule_for_current_riccati.rs_base.lo<=0:f.append(f'prefix {i}: invalid current active schedule')
        if any(r.lo<=0 for r in im.actual_rs_std_xyz_for_current_riccati):f.append(f'prefix {i}: invalid actual per-axis RS')
        n,m=_shape(c.master)
        if n==0 or n!=m:f.append(f'prefix {i}: master must be nonempty square')
        if not c.sectors or len(c.sectors)!=len(c.multipliers):f.append(f'prefix {i}: sectors/multipliers missing')
    return list(dict.fromkeys(f))

def certify(cells):
    f=validate_sequence(cells)
    if f:return {'closed':False,'validation_failures':f,'endpoint_closed':False,'every_prefix_closed':False,'prefixes':[]}
    rec=[]
    for i,c in enumerate(cells):
        ok,piv=ISS.certify_iss(c.master,c.sectors,c.multipliers,c.source_map,c.gamma_s,c.fp_map,c.gamma_n)
        rec.append({'prefix_length':i+1,'source_token':c.image.source_token,'predecessor_token':c.image.predecessor_token,'ldlt_closed':bool(ok),'pivot_lowers':piv})
    every=all(x['ldlt_closed'] for x in rec);endpoint=rec[-1]['ldlt_closed']
    return {'closed':bool(every and endpoint),'validation_failures':[],'endpoint_closed':endpoint,'every_prefix_closed':every,'prefixes':rec}

def build():
    j=JOINT.build();jf=JOINT.validate(j)
    if jf:raise RuntimeError('BRMM joint prerequisite failed: '+repr(jf))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','BRMM_attached_joint_image_required_at_every_prefix':True,'literal_previous_prefix_ancestry_required':True,'current_applied_schedule_retained_at_every_prefix':True,'raw_joint_f_tau_sigma_TS_RS_retained_at_every_prefix':True,'same_cell_augmented_master_sector_source_fp_maps_required':True,'full_interval_ISS_LDLT_backend_reused':True,'endpoint_cannot_substitute_for_every_prefix':True,'BIAS1_source_map_required_for_production':True,'binary32_fp_map_required_for_production':True,'production_complete_source_lineages_supplied_here':False,'production_endpoint_closed_here':False,'production_every_prefix_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'generate one PrefixInput per literal event prefix for every retained complete-BRMM/radial lineage, using the same source token for active schedule, P/H/R/K/bias/projection and both BIAS1/binary32 ISS maps, then require certify(...).closed for every lineage'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('BRMM_attached_joint_image_required_at_every_prefix','literal_previous_prefix_ancestry_required','current_applied_schedule_retained_at_every_prefix','raw_joint_f_tau_sigma_TS_RS_retained_at_every_prefix','same_cell_augmented_master_sector_source_fp_maps_required','full_interval_ISS_LDLT_backend_reused','endpoint_cannot_substitute_for_every_prefix','BIAS1_source_map_required_for_production','binary32_fp_map_required_for_production'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_complete_source_lineages_supplied_here','production_endpoint_closed_here','production_every_prefix_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'BRMM_prefix_bridge':d['BRMM_attached_joint_image_required_at_every_prefix'],'production_prefix':d['production_every_prefix_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
