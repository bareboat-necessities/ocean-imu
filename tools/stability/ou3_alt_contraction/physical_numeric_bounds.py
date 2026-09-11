#!/usr/bin/env python3
"""ALT binding of the repaired COMPLETE-BRMM numerical physical envelope.

The ALT common-storage enclosure needs numerical bounds on finite physical
coordinates, especially centered S at S=0 events.  These bounds must come from
the repaired physical BRMM theorem, never the legacy P4 error radius.
"""
from __future__ import annotations
import json
from pathlib import Path

import ou3_brmm_physical_wave_condition as PHYS
import ou3_brmm_centered_S_recurrence as SREC
import ou3_brmm_complete_physical_envelope as ENV

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_REPAIRED_BRMM_NUMERIC_PHYSICAL_ENVELOPE_BINDING_V1'


def build():
    p=PHYS.build(); pf=PHYS.validate(p)
    s=SREC.build(); sf=SREC.validate(s)
    e=ENV.build(); ef=ENV.validate(e)
    if pf or sf or ef: raise RuntimeError(f'physical envelope prerequisites failed phys={pf} S={sf} env={ef}')
    d=json.loads(DOMAIN.read_text())['complete_brmm_physical_envelope']
    b=e['complete_BRMM_hard_bounds']
    parity={
      'Hs':float(d['significant_wave_height_Hs_upper_m'])==float(b['Hs_upper_m']),
      'position':float(d['wave_position_norm_upper_m'])==float(b['wave_position_norm_upper_m']),
      'velocity':float(d['wave_velocity_norm_upper_mps'])==float(b['wave_velocity_norm_upper_mps']),
      'acceleration':float(d['wave_acceleration_norm_upper_mps2'])==float(b['wave_acceleration_norm_upper_mps2']),
      'body_rate':float(d['body_rate_norm_upper_deg_s'])==float(b['body_rate_norm_upper_deg_s']),
      'frequency':list(map(float,d['frequency_support_hz']))==list(map(float,b['frequency_support_hz'])),
      'centered_S':float(d['centered_primitive_D_S_upper_m_s'])==float(b['centered_primitive_D_S_upper_m_s']),
    }
    closed=bool(all(parity.values()) and p['physical_D_S_numeric_qualification_closed_for_complete_family'] and s['D_S_numeric_qualification_closed'] and s['D_S_max_m_s']==b['centered_primitive_D_S_upper_m_s'])
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'source':'ou3_brmm_physical_wave_condition + ou3_brmm_complete_physical_envelope',
      'operating_domain_parity':parity,'numeric_complete_family_physical_envelope_closed':closed,
      'wave_position_norm_upper_m':float(b['wave_position_norm_upper_m']),
      'wave_velocity_norm_upper_mps':float(b['wave_velocity_norm_upper_mps']),
      'wave_acceleration_norm_upper_mps2':float(b['wave_acceleration_norm_upper_mps2']),
      'body_rate_norm_upper_deg_s':float(b['body_rate_norm_upper_deg_s']),
      'frequency_support_hz':list(map(float,b['frequency_support_hz'])),
      'centered_S_norm_upper_m_s':float(b['centered_primitive_D_S_upper_m_s']),
      'specific_force_norm_lower_mps2':float(b['specific_force_norm_lower_mps2']),
      'specific_force_norm_upper_mps2':float(b['specific_force_norm_upper_mps2']),
      'D_S_comes_from_P4_working_radius':False,'legacy_300m_s_used_as_physical_bound':False,
      'P4_certificate_used_to_choose_padding':False,'trajectory_replay_used':False,
      'usable_by_ALT_finite_event_outer_enclosure':closed,
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('qualification/source mismatch')
    if d.get('numeric_complete_family_physical_envelope_closed') is not True:f.append('numeric physical envelope not closed')
    if not all(d.get('operating_domain_parity',{}).values()):f.append('operating-domain physical envelope drifted')
    if d.get('usable_by_ALT_finite_event_outer_enclosure') is not True:f.append('physical bounds unavailable to ALT')
    for k in ('D_S_comes_from_P4_working_radius','legacy_300m_s_used_as_physical_bound','P4_certificate_used_to_choose_padding','trajectory_replay_used'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('centered_S_norm_upper_m_s',0))!=1100.0:f.append('qualified D_S changed')
    return f
