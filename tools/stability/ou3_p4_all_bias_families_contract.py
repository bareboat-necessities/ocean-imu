#!/usr/bin/env python3
"""Fail-closed bias-family quantifier for all future OU-III P4 proof attempts.

P4 is not allowed to prove only a convenient bias model.  Every future proof
attempt must explicitly quantify over all three admitted bias families:

    BIAS0, BIAS1, BIAS2.

BIAS1 already has a concrete repository semantic binding (the physical one-root
one-parameter driver family and its joint ISS lift).  BIAS0 hardware
qualification and BIAS2 remain distinct required families; until authoritative
repository definitions are attached to this contract they are intentionally
reported as semantically unbound.  They may NOT be silently replaced by BIAS1,
collapsed into a generic bias box, or omitted from a P4 attempt.

This module distinguishes two questions:

1. attempt coverage: did the proof attempt declare all BIAS0/1/2 families?
2. semantic closure: is each declared family bound to its authoritative proof
   source and propagated through source cover, projection, transport, ISS and
   every-prefix LDLT?

The first is mandatory immediately.  The second is a promotion prerequisite.
Therefore this contract itself validates while P4 remains fail-closed until all
three semantic bindings and downstream closures are true.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias1_joint_iss_supply as BIAS1_ISS

SCHEMA=1
QUALIFICATION='OU3_P4_ALL_BIAS0_BIAS1_BIAS2_MANDATORY_CONTRACT_V1'
REQUIRED_BIAS_FAMILIES=('BIAS0','BIAS1','BIAS2')


def require_all_bias_families(families):
    declared=tuple(dict.fromkeys(str(x).upper() for x in families))
    missing=[x for x in REQUIRED_BIAS_FAMILIES if x not in declared]
    extra=[x for x in declared if x not in REQUIRED_BIAS_FAMILIES]
    if missing:
        raise RuntimeError('P4 proof attempt omitted mandatory bias families: '+','.join(missing))
    if extra:
        raise RuntimeError('P4 proof attempt declared unknown bias families: '+','.join(extra))
    return list(REQUIRED_BIAS_FAMILIES)


def attempt_coverage(families):
    try:
        declared=require_all_bias_families(families)
        return {'closed':True,'declared':declared,'missing':[]}
    except RuntimeError:
        d=tuple(dict.fromkeys(str(x).upper() for x in families))
        return {'closed':False,'declared':list(d),'missing':[x for x in REQUIRED_BIAS_FAMILIES if x not in d]}


def promotion_bias_gate(*,families,semantic_bindings,source_cover,projection_transport,prefix_iss_ldlt):
    require_all_bias_families(families)
    tables=(semantic_bindings,source_cover,projection_transport,prefix_iss_ldlt)
    for table in tables:
        if set(table)!=set(REQUIRED_BIAS_FAMILIES):
            raise RuntimeError('all-bias promotion table must contain exactly BIAS0/BIAS1/BIAS2')
    return all(bool(table[x]) for table in tables for x in REQUIRED_BIAS_FAMILIES)


def build():
    b1=BIAS1.build();b1f=BIAS1.validate(b1)
    iss=BIAS1_ISS.build();issf=BIAS1_ISS.validate(iss)
    if b1f or issf:raise RuntimeError(f'BIAS1 binding prerequisite failed family={b1f} ISS={issf}')
    required=require_all_bias_families(REQUIRED_BIAS_FAMILIES)
    semantic={'BIAS0':False,'BIAS1':True,'BIAS2':False}
    source={'BIAS0':False,'BIAS1':False,'BIAS2':False}
    projection={'BIAS0':False,'BIAS1':False,'BIAS2':False}
    prefix={'BIAS0':False,'BIAS1':False,'BIAS2':False}
    promotion=promotion_bias_gate(families=required,semantic_bindings=semantic,
                                  source_cover=source,projection_transport=projection,
                                  prefix_iss_ldlt=prefix)
    omitted0=attempt_coverage(('BIAS1','BIAS2'))
    omitted1=attempt_coverage(('BIAS0','BIAS2'))
    omitted2=attempt_coverage(('BIAS0','BIAS1'))
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'required_bias_families':required,'required_bias_family_count':3,
      'all_future_P4_attempts_must_declare_BIAS0_BIAS1_BIAS2':True,
      'single_or_subset_bias_family_attempt_forbidden':True,
      'omitting_BIAS0_is_rejected':not omitted0['closed'],
      'omitting_BIAS1_is_rejected':not omitted1['closed'],
      'omitting_BIAS2_is_rejected':not omitted2['closed'],
      'BIAS1_authoritative_family_module_bound':True,
      'BIAS1_joint_ISS_supply_module_bound':True,
      'BIAS1_one_root_one_parameter_history_retained':bool(b1['one_root_one_parameter_history_required']),
      'BIAS0_authoritative_semantic_binding_closed_here':False,
      'BIAS2_authoritative_semantic_binding_closed_here':False,
      'BIAS0_may_be_substituted_by_BIAS1':False,'BIAS2_may_be_substituted_by_BIAS1':False,
      'generic_bias_box_may_replace_three_family_quantifier':False,
      'semantic_bindings':semantic,'source_cover_closure':source,
      'projection_transport_closure':projection,'every_prefix_ISS_LDLT_closure':prefix,
      'all_three_semantically_bound_here':all(semantic.values()),
      'all_three_source_cover_closed_here':all(source.values()),
      'all_three_projection_transport_closed_here':all(projection.values()),
      'all_three_every_prefix_ISS_LDLT_closed_here':all(prefix.values()),
      'all_bias_promotion_gate_closed':promotion,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'attach authoritative BIAS0 and BIAS2 definitions, then require each BIAS0/1/2 family independently through same-history source cover, projection/bias transport, finite-precision ISS maps, endpoint and every-prefix augmented LDLT before any P4 promotion'}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('required_bias_families')!=list(REQUIRED_BIAS_FAMILIES):f.append('required bias family set changed')
    for k in ('all_future_P4_attempts_must_declare_BIAS0_BIAS1_BIAS2','single_or_subset_bias_family_attempt_forbidden','omitting_BIAS0_is_rejected','omitting_BIAS1_is_rejected','omitting_BIAS2_is_rejected','BIAS1_authoritative_family_module_bound','BIAS1_joint_ISS_supply_module_bound','BIAS1_one_root_one_parameter_history_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('BIAS0_authoritative_semantic_binding_closed_here','BIAS2_authoritative_semantic_binding_closed_here','BIAS0_may_be_substituted_by_BIAS1','BIAS2_may_be_substituted_by_BIAS1','generic_bias_box_may_replace_three_family_quantifier','all_three_semantically_bound_here','all_three_source_cover_closed_here','all_three_projection_transport_closed_here','all_three_every_prefix_ISS_LDLT_closed_here','all_bias_promotion_gate_closed','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for name in REQUIRED_BIAS_FAMILIES:
        for table in ('semantic_bindings','source_cover_closure','projection_transport_closure','every_prefix_ISS_LDLT_closure'):
            if name not in d.get(table,{}):f.append(table+' missing '+name)
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'required':d['required_bias_families'],'BIAS0_bound':d['BIAS0_authoritative_semantic_binding_closed_here'],'BIAS1_bound':d['semantic_bindings']['BIAS1'],'BIAS2_bound':d['BIAS2_authoritative_semantic_binding_closed_here'],'promotion':d['all_bias_promotion_gate_closed'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
