#!/usr/bin/env python3
"""Count-free asynchronous-event closure for the ALT Live theorem.

Normal Live deliberately does not assume a magnetometer ODR as a theorem input.
The typed shipping kernel accepts a finite tuple of asynchronous magnetometer
events after an IMU sample and executes *all* of them.  Therefore an endpoint
cover that merely chooses zero/one magnetic update per IMU sample is not
universal.

The correct count-free composition is a Kleene-star argument.  Let V be the
same common storage used by the IMU-word certificate.  If every admitted
accepted magnetometer event satisfies

    V(M_m(x)) <= V(x)

on the required chart/domain, then by finite induction every finite sequence
M_{m_q} o ... o M_{m_1} is nonexpansive.  Rejected/not-due packets are identity
and satisfy the same inequality.  Thus arbitrary finite asynchronous insertions
do not require an event-count upper bound and cannot worsen the strict
contraction supplied by the IMU word.

This module closes only the logical composition rule.  The source-uniform
nonexpansive inequality for one physical magnetometer event remains a numerical
common-storage obligation and must be certified with the same M; no ODR or
finite packet-count assumption may replace it.
"""
from __future__ import annotations

import json
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
KERNEL=REPO/'tools/stability/ou3_brmm_complete_window_execution_kernel.py'
QUALIFICATION='OU3_ALT_ASYNC_MAGNETOMETER_KLEENE_STAR_V1'


def finite_star_rho(single_event_rho: float, count: int) -> float:
    """Composition bound for diagnostic/unit algebra; theorem uses rho<=1."""
    r=float(single_event_rho); q=int(count)
    if not (0.0 <= r <= 1.0) or q < 0: raise ValueError('require rho in [0,1], finite count >=0')
    return r**q


def build():
    live=json.loads(DOMAIN.read_text())['normal_live']; text=KERNEL.read_text()
    typed_tuple='magnetometer_events_after_imu: tuple[MagneticEvent, ...]' in text
    all_loop='for magnetic_event_index, event in enumerate(sample.magnetometer_events_after_imu):' in text
    no_odr=bool(live.get('magnetic_reference_handling'))
    # PE module/domain deliberately supplies recurrence, not a count cap.
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'typed_kernel_accepts_finite_async_event_tuple':typed_tuple,
      'typed_kernel_executes_every_async_event':all_loop,
      'magnetometer_event_count_upper_bound_assumed':False,
      'magnetometer_ODR_used_as_stability_assumption':False,
      'zero_or_one_mag_per_IMU_is_universal_cover':False,
      'rejected_or_not_due_event_map_is_identity':True,
      'finite_star_induction_rule':'V(M_m x)<=V(x) for every admitted m => V(M_q...M_1 x)<=V(x) for every finite q',
      'finite_star_composition_theorem_closed':bool(typed_tuple and all_loop and no_odr),
      'same_common_storage_required_for_mag_and_IMU_word':True,
      'single_source_uniform_magnetometer_nonexpansive_certificate_closed':False,
      'arbitrary_finite_async_magnetometer_sequence_closed':False,
      'ALT_LIVE_PASS':False,
      'next_obligation':'with the eventual common M, certify each admitted physical magnetometer event nonexpansive source-uniformly; then finite induction closes every asynchronous event tuple without any ODR/count assumption',
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('typed_kernel_accepts_finite_async_event_tuple','typed_kernel_executes_every_async_event','rejected_or_not_due_event_map_is_identity','finite_star_composition_theorem_closed','same_common_storage_required_for_mag_and_IMU_word'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('magnetometer_event_count_upper_bound_assumed','magnetometer_ODR_used_as_stability_assumption','zero_or_one_mag_per_IMU_is_universal_cover','single_source_uniform_magnetometer_nonexpansive_certificate_closed','arbitrary_finite_async_magnetometer_sequence_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
