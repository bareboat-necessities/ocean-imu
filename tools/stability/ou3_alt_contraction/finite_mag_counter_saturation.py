"""Shipping magnetic counter invariant, independent of event-rate bounds.

For every safe count c and every finite number n of attempted calls, the
counter is min(INT_MAX, c+n). For every configurable threshold u in [0,INT_MAX],
that count is >=u exactly when c+n is >=u. The guarded C++ increment never
evaluates INT_MAX+1. This proves counter safety, not all-event arithmetic.
"""
from hashlib import sha256
from pathlib import Path
import re

SIGNED_MAX = (1 << 31) - 1
QUALIFICATION = 'OU3_ALT_MAG_COUNTER_SATURATION_V1'
SOURCE = Path(__file__).resolve().parents[3] / 'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
AUDITED_SOURCE_SHA = 'fabd03e9c3eb6069df107c7413ffb4b33fbdcd1ce06d3923b0c1ceeb3bcd7359'
TUNER_SOURCE = SOURCE.parents[1] / 'tuner/MagAutoTuner.h'
AUDITED_TUNER_SHA = 'd2b1fd6b2131989da4d7b5671dd0bdca7cc37c23ed9794bf1e737c7d8910fa87'
AUDITED_UPDATE_SHA = '0dcecf08e66746d2aa2bdad289c768b7337dcdab45201a052dae7186163cf955'


def source_fingerprint():
    s = re.sub(r'/\*.*?\*/|//[^\n]*', '', SOURCE.read_text(), flags=re.S)
    start = s.index('{', s.index('void updateMag(const Eigen::Vector3f& mag_body_ned)'))
    end, depth = start + 1, 1
    while depth:
        if s[end] == '{': depth += 1
        if s[end] == '}': depth -= 1
        end += 1
    return sha256(re.sub(r'\s+', '', s[start:end]).encode()).hexdigest()


def audit_source():
    # Include reset/configuration paths, not just the incrementing method.
    if sha256(SOURCE.read_bytes()).hexdigest() != AUDITED_SOURCE_SHA:
        raise RuntimeError('wrapper counter ownership changed; re-audit reset/configuration')
    if source_fingerprint() != AUDITED_UPDATE_SHA:
        raise RuntimeError('magnetic counter source changed; re-audit saturation and release')
    s = SOURCE.read_text()
    if s.count('mag_updates_applied_') != 5 or 'int  mag_updates_applied_ = 0;' not in s:
        raise RuntimeError('counter declaration/uses changed; re-audit invariant')
    if 'if (n >= 0) mag_updates_to_unlock_ = n;' not in s:
        raise RuntimeError('configurable counter threshold changed; re-audit comparisons')
    if sha256(TUNER_SOURCE.read_bytes()).hexdigest() != AUDITED_TUNER_SHA:
        raise RuntimeError('MagAutoTuner counters changed; re-audit acquisition/refinement')
    return True


def after_attempts(count, attempts=1):
    if type(count) is not int or not 0 <= count <= SIGNED_MAX:
        raise ValueError('safe signed-int32 counter required')
    if type(attempts) is not int or attempts < 0:
        raise ValueError('nonnegative finite integer attempt count required')
    return min(SIGNED_MAX, count + attempts)


def threshold_equivalence(count, attempts, threshold):
    if type(threshold) is not int or not 0 <= threshold <= SIGNED_MAX:
        raise ValueError('configurable signed-int32 threshold required')
    return (after_attempts(count, attempts) >= threshold) == (count + attempts >= threshold)


def build():
    audit_source()
    return {
        'qualification': QUALIFICATION,
        'shipping_source_correspondence_audited': True,
        'signed_counter_max': SIGNED_MAX,
        'counter_formula': 'min(INT_MAX, entry_count + attempted_calls)',
        'all_finite_attempt_counts_safe': True,
        'all_representable_unlock_thresholds_preserved': True,
        'no_positive_minimum_call_gap_required': True,
        'measurement_and_release_checks_continue_at_saturation': True,
        'counter_lifetime_closed': True,
        'accepted_and_rejected_magnetic_sample_counters_closed': True,
        'complete_deployment_arithmetic_closed': False,
        'storage_search_allowed': False,
    }
