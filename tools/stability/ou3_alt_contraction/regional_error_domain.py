#!/usr/bin/env python3
"""ALT binding of the declared full-scale regional Normal-Live error domain.

The P4 closure-domain JSON contains deterministic error-coordinate working scales
that are explicitly *not* covariance-confidence bounds and *not* a replay fit.
Those constants are useful to define the regional Live theorem, but ALT must not
call ``ou3_p4_hard_entry_set.build()`` because that builder also executes a
fresh-startup reachability proof.  Startup is a separate ALT phase.

This module reads only the declared working-domain constants, verifies that no
scale shrink is allowed, and exposes them as a conditional regional hypothesis.
It makes no claim that startup has entered the domain.
"""
from __future__ import annotations
import json,math
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'
QUALIFICATION='OU3_ALT_REGIONAL_NORMAL_LIVE_ERROR_DOMAIN_V1'
EXPECTED=(
 'attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m',
 'integral_displacement_norm_m_s','latent_acceleration_norm_mps2','accelerometer_bias_error_norm_mps2')

def build(path:Path=CLOSURE):
    root=json.loads(Path(path).read_text());h=root['hard_entry_search'];scale=float(h['minimum_certified_scale']);candidates=list(map(float,h['candidate_scale_factors']));r={k:math.nextafter(float(v),math.inf) for k,v in h['base_coordinate_radii'].items()}
    if scale!=1.0 or candidates!=[1.0] or h.get('requires_full_declared_scale') is not True:raise RuntimeError('declared regional error domain was shrunk')
    if tuple(r.keys())!=EXPECTED:raise RuntimeError('regional error-coordinate schema changed')
    if any(not math.isfinite(v) or v<=0 for v in r.values()):raise RuntimeError('regional error radius invalid')
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','scope':'conditional regional Normal-Live theorem domain','coordinate_radii':r,'full_declared_scale':scale,'proof_error_domain_shrunk':False,'trajectory_fit':False,'covariance_membership_used':False,'fresh_startup_reachability_consumed':False,'startup_membership_claimed_here':False,'regional_membership_is_theorem_hypothesis':True,'physical_S_bound_taken_from_this_error_radius':False,'legacy_S_error_radius_m_s':r['integral_displacement_norm_m_s'],'ALT_STARTUP_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('qualification/source mismatch')
    for k in ('regional_membership_is_theorem_hypothesis',):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('proof_error_domain_shrunk','trajectory_fit','covariance_membership_used','fresh_startup_reachability_consumed','startup_membership_claimed_here','physical_S_bound_taken_from_this_error_radius','ALT_STARTUP_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('full_declared_scale',0))!=1.0:f.append('regional scale changed')
    r=d.get('coordinate_radii',{})
    if tuple(r.keys())!=EXPECTED:f.append('regional coordinate schema mismatch')
    return f
