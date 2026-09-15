"""Exact finite-word totality obstruction and distinct real-caller count bound.

No billions-of-calls replay is needed. A compressed schedule and induction on
the literal counter projection prove the contradiction. Native tests check the
last two edges only; they do not pretend to reach INT_MAX from reset.
See docs/ou3-alt-deployment-prerequisite.md for the quantified argument.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import re

from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as CALLS
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

ROOT = Path(__file__).resolve().parents[3]
HEADER = ROOT / 'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SKETCH = ROOT / 'sensors/full_marine_ins/atomS3R_ins_kalman_ou3/atomS3R_ins_kalman_ou3.ino'
QUALIFICATION = 'OU3_ALT_FINITE_MAG_COUNTER_OBSTRUCTION_V1'


@dataclass(frozen=True)
class Burst:
    """All-time, locally finite schedule relative to an eligible gauged Live.

    Calls 1..size occur at first_time. Thereafter call size+j occurs at
    first_time+j*gap_max. Entry counter and horizon are theorem operands,
    not changes to the shipping state or physical source class.
    """
    entry_count: int = 0
    horizon: F = CALLS.CANONICAL_ALT_WINDOW
    schedule: CALLS.Schedule = CALLS.Schedule()

    def __post_init__(self):
        if type(self.entry_count) is not int or not 0 <= self.entry_count <= GATE.SIGNED_COUNTER_MAX:
            raise ValueError('safe signed int32 entry count required')
        if not isinstance(self.schedule, CALLS.Schedule):
            raise TypeError('current named call schedule required')
        object.__setattr__(self, 'horizon', F(self.horizon))
        if self.horizon <= 0:
            raise ValueError('positive finite horizon required')

    @property
    def size(self):
        return GATE.SIGNED_COUNTER_MAX + 1 - self.entry_count

    @property
    def first_time(self):
        # The master explicitly admits MAG at the fresh physical endpoint,
        # before transition 1. Do not invent a between-grid source endpoint.
        return F(0)

    def time(self, ordinal):
        if type(ordinal) is not int or ordinal < 1:
            raise ValueError('positive call ordinal required')
        return self.first_time + max(0, ordinal - self.size) * self.schedule.gap_max

    def projected_count_after_defined_prefix(self, calls):
        """Induction consequence if those calls return without earlier UB.

        At each eligible return the source increments exactly once. This
        expression is valid only through INT_MAX, not a wraparound model.
        """
        if type(calls) is not int or not 0 <= calls < self.size:
            raise ValueError('prefix must stop before the undefined increment')
        return self.entry_count + calls

    def calls_through(self, time):
        t = F(time)
        if t < self.first_time:
            return 0
        return self.size + int((t - self.first_time) // self.schedule.gap_max)


def _body(source, signature, occurrence=0):
    """Extract a named brace body after removing comments (audited sources)."""
    source = re.sub(r'/\*.*?\*/|//[^\n]*', '', source, flags=re.S)
    matches = list(re.finditer(re.escape(signature), source))
    if occurrence >= len(matches):
        raise RuntimeError('audited source signature missing: ' + signature)
    start = source.index('{', matches[occurrence].end())
    depth = 1
    end = start + 1
    while depth:
        if source[end] == '{': depth += 1
        if source[end] == '}': depth -= 1
        end += 1
    return re.sub(r'\s+', '', source[start:end])


# Bind the manual source-order argument to exactly the inspected functions.
# These hashes ignore comments/spacing; a semantic edit requires a fresh audit.
AUDITED = {
    'inner_updateMag': 'b9fcaab91834f76bb826c406aaebb33eb0a67406b17801e1331ca67ccf4e6078',
    'outer_updateMag': 'e19383628b13245cb0134f6d953ae7bcaee0517352079b3593fe16fcf4350ca4',
    'caller_updateFilter': 'f8e240d00ee7296b0f53f9707740c36a63030463b631ea9c46637a56656132d4',
    'caller_fresh_gate': '612ef91785a565a2b89b6cac2a77f5d93e60819ec487091c732ea05c48a26db9',
    'caller_dt': '7d387fd41e167c4b92ff79527aa2ad0282fb3afd03cb5113f83fbf01066d3e7b',
    'caller_reset': '04b78bdef917741c8def2ad9017f53f6d040b89637e3caa0f3881b973b1cd9f9',
}


def source_fingerprints():
    h, s = HEADER.read_text(), SKETCH.read_text()
    bodies = {
        'inner_updateMag': _body(h, 'void updateMag(const Eigen::Vector3f& mag_body_ned)', 0),
        'outer_updateMag': _body(h, 'void updateMag(const Eigen::Vector3f& mag_body_ned)', 1),
        'caller_updateFilter': _body(s, 'void updateFilter_(const ImuSample& s)'),
        'caller_fresh_gate': _body(s, 'bool updateMagFreshGate_(bool mag_ok, uint32_t now_ms)'),
        'caller_dt': _body(s, 'float computeFusionDtFromSampleTimestamp_(const ImuSample& s)'),
        'caller_reset': _body(s, 'void resetFusion_()'),
    }
    return {k: sha256(v.encode()).hexdigest() for k, v in bodies.items()}


def audit_source():
    if source_fingerprints() != AUDITED:
        raise RuntimeError('magnetic counter/caller source changed; re-audit the totality argument')
    h = HEADER.read_text()
    if h.count('mag_updates_applied_') != 4 or 'int  mag_updates_applied_ = 0;' not in h:
        raise RuntimeError('magnetic counter declaration/uses changed; re-audit induction')
    # Every write to this private counter is the declaration, reset to zero,
    # or the audited unchecked increment. The reset is not a per-call action.
    writes = re.findall(r'mag_updates_applied_\s*(?:\+\+|=[^;]*)\s*;', h)
    if [re.sub(r'\s+', '', w) for w in writes] != [
            'mag_updates_applied_++;', 'mag_updates_applied_=0;', 'mag_updates_applied_=0;']:
        raise RuntimeError('magnetic counter writes changed; re-audit induction')
    s = SKETCH.read_text()
    if s.count('fusion_.updateMag(') != 1 or s.count('fusion_.update(') != 1:
        raise RuntimeError('caller filter call sites changed; re-audit invocation bound')
    return True


def caller_count_bound(imu_invocations):
    """Reset-rooted count <= invocations for this sketch, if execution is defined.

    One sequential updateFilter_ invocation calls update once and updateMag at
    most once. Gates can remove increments. This says nothing about elapsed
    time, universal startup termination, other API callers, or target libm.
    """
    audit_source()
    if type(imu_invocations) is not int or imu_invocations < 0:
        raise ValueError('nonnegative invocation count required')
    return imu_invocations


def build():
    audit_source()
    b = Burst()
    CALLS.check_prefix(CALLS.Clock(0), b.schedule, time=b.first_time)
    if b.projected_count_after_defined_prefix(b.size - 1) != GATE.SIGNED_COUNTER_MAX:
        raise AssertionError('counter induction endpoint mismatch')
    return {
        'qualification': QUALIFICATION,
        'source_correspondence_audited': True,
        'schedule': b.schedule.assumption_id,
        'entry_count': b.entry_count,
        'burst_calls': b.size,
        'burst_physical_offset_s': b.first_time,
        'finite_horizon_s': b.horizon,
        'count_before_last_call': GATE.SIGNED_COUNTER_MAX,
        'required_last_integer_successor': GATE.SIGNED_COUNTER_MAX + 1,
        'locally_finite_time_unbounded_schedule_constructed': True,
        'same_timestamp_calls_use_same_source_sample': True,
        'earlier_undefined_execution_also_defeats_totality': True,
        'total_defined_execution_for_current_schedule_possible': False,
        'conditional_on_eligible_gauged_Live_entry': True,
        'actual_sketch_mag_calls_per_imu_invocation_upper': 1,
        'actual_sketch_conditional_30602_invocation_count_upper': caller_count_bound(CLOCK.MAX_STEPS),
        'actual_sketch_35ms_gate_is_unconditional_min_gap': False,
        'actual_sketch_measured_dt_is_canonical_exact_5ms': False,
        'actual_caller_substituted_for_generic_schedule': False,
        'target_compiler_Eigen_libm_correspondence_closed': False,
        'shipping_filter_instability_claimed': False,
        'storage_search_allowed': False,
    }
