"""Explicit conditional ALT service class; membership never proves stability.

Parameters are theorem hypotheses, not inferred from callback cadence or fitted
trace extrema. Unconditional COMPLETE-BRMM and the ungauged theorem survive.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class InformativeServiceAssumption:
    max_gap_s: float
    superword_s: float
    heading_response_floor: float
    pair_information_floor: float

    def __post_init__(self):
        values = (self.max_gap_s, self.superword_s, self.heading_response_floor,
                  self.pair_information_floor)
        if any(not math.isfinite(x) or x <= 0 for x in values):
            raise ValueError('strictly positive finite service parameters required')
        if self.max_gap_s > self.superword_s / 2:
            raise ValueError('superword must span at least two maximum service gaps')

    def declaration(self):
        return {
            'qualification': 'OU3_ALT_CONDITIONAL_INFORMATIVE_SERVICE_V1',
            'parameters': dict(vars(self)),
            'assumption_explicitly_selected': True,
            'starts_after_actual_live_and_heading_acquisition': True,
            'quantifier': 'every overlapping superword within the declared service regime',
            'accepted_events_required': True,
            'transport': 'actual whitened H times preceding full event product, in one root heading/bias frame',
            'frame': 'orthonormal heading and axial gyro-bias coordinates in declared rad and rad/s units',
            'gramian_requirement': 'sum B_i.T B_i >= pair_information_floor * I2',
            'boundary_gaps_included': True,
            'isolated_finite_windows_imply_recurrence': False,
            'call_cadence_implies_membership': False,
            'unconditional_source_class_restricted': False,
            'uniform_metric_coercivity_certified': False,
            'source_uniform_rho_certified': False,
            'storage_search_allowed': False,
            'ALT_STARTUP_PASS': False, 'ALT_LIVE_PASS': False,
            'ALT_END_TO_END_PASS': False,
        }

    def audit_native_window(self, row, dt_s):
        """Point audit only; accepts parser-produced native rows, not certificates."""
        duration = (row['end_sample']-row['start_sample'])*dt_s
        vals = (duration, row['max_observed_magnetic_gap_s'],
                row['minimum_transported_heading_response'],
                row['transported_heading_bias_min_eigenvalue'])
        if any(not math.isfinite(x) or x < 0 for x in vals):
            raise ValueError('finite nonnegative native measurements required')
        # A shorter interval may satisfy a stronger instance, but is not an
        # audit of the requested superword. Binary32 sample time is explicit.
        duration_matches = math.isclose(duration, self.superword_s, rel_tol=0, abs_tol=1e-6)
        passed = (duration_matches and row['accepted_magnetic_events'] >= 2
                  and vals[1] <= self.max_gap_s
                  and vals[2] >= self.heading_response_floor
                  and vals[3] >= self.pair_information_floor)
        return {'finite_window_membership_point_pass': bool(passed),
                'duration_matches': duration_matches,
                'infinite_recurrence_certified': False,
                'source_uniform_service_certified': False,
                'stability_implied_by_membership': False}
