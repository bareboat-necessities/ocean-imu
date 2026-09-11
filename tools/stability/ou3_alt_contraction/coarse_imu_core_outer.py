#!/usr/bin/env python3
"""Universal coarse outer for one Normal-Live IMU core, excluding async mag.

The shipping typed kernel has exactly one prediction and one accelerometer
update per valid Normal-Live IMU sample, and a boolean S_zero_due flag.  Thus
the homogeneous per-sample IMU core has exactly two literal event patterns:

    prediction -> accelerometer
    prediction -> S_zero -> accelerometer

Asynchronous magnetometer calls are *not* part of this finite family because
the theorem does not impose an event-count upper bound.  They are handled by
``asynchronous_event_star`` and must later be proved nonexpansive in the same
common storage.

This module reuses ``coarse_universal_local_outer`` only for the certified
coordinate-projection evaluation of the individual local shipping Jacobians.
It never consumes that module's old zero/one-magnetometer sample-family field.
"""
from __future__ import annotations

from ou3_interval import matrix_identity
from tools.stability.ou3_alt_contraction import coarse_universal_local_outer as LOCAL
from tools.stability.ou3_alt_contraction import asynchronous_event_star as MAGSTAR

QUALIFICATION='OU3_ALT_COARSE_UNIVERSAL_IMU_CORE_OUTER_V1'


def _compose(*maps):
    A=matrix_identity(24)
    from ou3_interval import matrix_mul
    for J in maps:A=matrix_mul(J,A)
    return A


def build():
    d=LOCAL.build();df=LOCAL.validate(d);star=MAGSTAR.build();sf=MAGSTAR.validate(star)
    if df or sf:raise RuntimeError(f'IMU-core prerequisites failed local={df} star={sf}')
    reports={}
    for mode,rows in d['reports'].items():
        reports[mode]={}
        for family,row in rows.items():
            finite=bool(row.get('prediction_finite') and row.get('S_zero_finite') and row.get('accelerometer_finite'))
            fam={'bias_family':family,'failures':[x for x in row['failures'] if not x.startswith('magnetometer:')],
                 'magnetometer_failures':[x for x in row['failures'] if x.startswith('magnetometer:')],
                 'imu_core_family_finite':finite,'imu_core_family':[]}
            if finite:
                P=row['maps']['prediction'];S=row['maps']['S_zero'];A=row['maps']['accelerometer'];Id=matrix_identity(24)
                fam['imu_core_family']=[_compose(P,s,A) for s in (Id,S)]
            reports[mode][family]=fam
    closed=all(r['imu_core_family_finite'] for rows in reports.values() for r in rows.values())
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'one_prediction_per_IMU':True,'one_required_accelerometer_per_IMU':True,
      'S_due_boolean_gives_zero_or_one_S_per_IMU':True,
      'async_magnetometer_excluded_from_finite_IMU_core_family':True,
      'async_magnetometer_count_upper_bound_assumed':False,
      'async_magnetometer_star_theorem_consumed':star['finite_star_composition_theorem_closed'],
      'reports':reports,'coarse_universal_IMU_core_family_materialized':closed,
      'async_magnetometer_nonexpansive_same_M_closed':False,
      'complete_async_sample_family_closed':False,'ALT_LIVE_PASS':False,
      'next_obligation':('induct the two-pattern IMU-core outer over 600 steps and search common M; with that same M separately certify every source-uniform magnetometer event nonexpansive' if closed else 'split whichever prediction/S/accelerometer certified projection failed; magnetic failures are separate star obligation')
    }


def summary(d):
    return {mode:{fam:{'failures':r['failures'],'magnetometer_failures':r['magnetometer_failures'],'imu_core_family_finite':r['imu_core_family_finite']} for fam,r in rows.items()} for mode,rows in d['reports'].items()}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('one_prediction_per_IMU','one_required_accelerometer_per_IMU','S_due_boolean_gives_zero_or_one_S_per_IMU','async_magnetometer_excluded_from_finite_IMU_core_family','async_magnetometer_star_theorem_consumed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('async_magnetometer_count_upper_bound_assumed','async_magnetometer_nonexpansive_same_M_closed','complete_async_sample_family_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
