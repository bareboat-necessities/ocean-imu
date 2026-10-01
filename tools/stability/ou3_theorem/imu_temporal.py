"""Conditional, deterministic SLOW + FAST error calculus (proof only).

No horizon/cap is selected here. Numerical device qualification is OPEN in
constants.json. The fast signal is the hold of CALIBRATED DELIVERED samples;
analog cancellation before filtering/decimation does not qualify this signal.
See docs/ou3-imu-two-timescale.md, SF1--SF10. Floating evaluation of these
identities is diagnostic, not an interval or all-time hardware certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


class OpenTemporalQualification(ValueError):
    """A required fast temporal envelope has not been supplied."""


def _nonnegative(*values: float) -> None:
    if not all(math.isfinite(x) and x >= 0 for x in values):
        raise ValueError('finite nonnegative bounds required')


@dataclass(frozen=True)
class FastWindow:
    """Every window 0<T<=horizon_s obeys |integral b_f|<=min(B_f T, cap).

    `cap` has units error*seconds (m/s for accel, radians for gyro).
    This is an externally qualified parameter, NOT a covariance/noise sigma.
    """
    horizon_s: float
    cap: float

    def __post_init__(self) -> None:
        _nonnegative(self.horizon_s, self.cap)
        if self.horizon_s == 0:
            raise ValueError('positive finite fast horizon required')

    def validate_amplitude(self, amplitude: float) -> None:
        _nonnegative(amplitude)
        product = amplitude * self.horizon_s
        if not math.isfinite(product):
            raise ValueError('fast horizon/amplitude product overflow')
        if (amplitude == 0 and self.cap != 0) or (amplitude > 0 and self.cap >= product):
            raise ValueError('fast qualification must improve on the amplitude-only horizon box')

    def accumulation(self, duration_s: float, amplitude: float) -> float:
        """Valid tiled outer envelope for ANY duration; never resets the signal."""
        self.validate_amplitude(amplitude)
        _nonnegative(duration_s)
        if duration_s == 0 or amplitude == 0:
            return 0.0
        n, rem = divmod(duration_s, self.horizon_s)
        return min(amplitude * duration_s, n * self.cap + min(amplitude * rem, self.cap))

    def cell_amplitude(self, cell_s: float, amplitude: float) -> float:
        """Exact maximal value on one complete, constant delivered-sample cell."""
        self.validate_amplitude(amplitude)
        _nonnegative(cell_s)
        if cell_s == 0:
            raise ValueError('positive complete sample-cell duration required')
        return min(amplitude, self.cap / min(cell_s, self.horizon_s))


def slow_change(amplitude: float, rate: float, h_s: float) -> float:
    _nonnegative(amplitude, rate, h_s)
    return min(2 * amplitude, rate * h_s)


def two_epoch_envelope(slow_amplitude: float, slow_rate: float,
                       fast_amplitude: float, h_s: float, *,
                       fast: FastWindow | None,
                       first_cell_s: float, second_cell_s: float) -> dict:
    """SF3: two epochs on DISTINCT complete cells, or the identical epoch h=0.

    The exact fast-class supremum is attained by two oppositely signed cells
    with zero fast error elsewhere. Unknown temporal parameters do not certify
    that the amplitude-only outer relaxation is reachable. A fixed predecessor
    can further shrink this envelope; preserve its full reachable set.
    """
    _nonnegative(slow_amplitude, slow_rate, fast_amplitude, h_s,
                 first_cell_s, second_cell_s)
    if min(first_cell_s, second_cell_s) <= 0:
        raise ValueError('positive complete cells required')
    if 0 < h_s < first_cell_s:
        raise ValueError('distinct sample epochs cannot lie inside the first hold cell')
    ds = slow_change(slow_amplitude, slow_rate, h_s)
    if h_s == 0:
        df = 0.0
    elif fast is None:
        df = None
    else:
        df = (fast.cell_amplitude(first_cell_s, fast_amplitude)
              + fast.cell_amplitude(second_cell_s, fast_amplitude))
    return {'slow_change': ds, 'fast_change': df,
            'total_change': None if df is None else ds + df,
            'amplitude_only_outer': 0.0 if h_s == 0 else ds + 2 * fast_amplitude,
            'temporal_qualification': 'OPEN' if fast is None else 'CONDITIONAL',
            'certifies_physical_history': False}


def averaged_fast_difference(fast: FastWindow | None, amplitude: float,
                             averaging_s: float, separation_s: float) -> float:
    """SF4: overlapping boxcar identity, including both boundary intervals."""
    _nonnegative(amplitude, averaging_s, separation_s)
    if averaging_s == 0:
        raise ValueError('positive averaging interval required')
    if fast is None:
        raise OpenTemporalQualification('both the horizon and accumulation cap are required')
    return 2 * min(fast.accumulation(averaging_s, amplitude),
                   fast.accumulation(separation_s, amplitude)) / averaging_s


def fast_weighted_outer(coefficients, steps_s: Sequence[float], amplitude: float,
                        fast: FastWindow | None) -> float:
    """SF5 signed Abel bound on sum A_k b_f,k, with matrix/rotating weights.

    A_k already includes the actual sample weight, frame and transported gain.
    Dividing by dt BEFORE differencing is essential. The endpoint term cannot
    be dropped. A finite supplied profile makes this a conditional algebraic
    bound, not a declaration that a particular device is qualified.
    """
    import numpy as np
    if fast is None:
        raise OpenTemporalQualification('cannot replace missing temporal evidence by independent boxes')
    a = np.asarray(coefficients, dtype=float)
    dt = np.asarray(steps_s, dtype=float)
    if (a.ndim != 3 or a.shape[0] != len(dt) or a.shape[2] != 3 or len(dt) == 0
            or a.shape[1] == 0 or not np.isfinite(a).all()
            or not np.isfinite(dt).all() or np.any(dt <= 0)):
        raise ValueError('one finite output-by-3 matrix per positive complete cell required')
    q = a / dt[:, None, None]
    elapsed = np.cumsum(dt)
    opnorm = lambda x: float(np.linalg.norm(x, 2))
    bound = opnorm(q[-1]) * fast.accumulation(float(elapsed[-1]), amplitude)
    for k in range(1, len(dt)):
        bound += opnorm(q[k-1] - q[k]) * fast.accumulation(float(elapsed[k-1]), amplitude)
    pointwise = sum(opnorm(ak) * fast.cell_amplitude(float(h), amplitude) for ak,h in zip(a,dt))
    return min(bound, pointwise)


def slow_weighted_outer(coefficients, steps_between_s: Sequence[float],
                        amplitude: float, rate: float) -> float:
    """SF5 linked predecessor form; a sound outer bound, not independent truth."""
    import numpy as np
    _nonnegative(amplitude, rate)
    a = np.asarray(coefficients, dtype=float)
    dt = np.asarray(steps_between_s, dtype=float)
    if (a.ndim != 3 or a.shape[0] != len(dt)+1 or a.shape[2] != 3 or a.shape[1] == 0
            or not np.isfinite(a).all() or not np.isfinite(dt).all() or np.any(dt <= 0)):
        raise ValueError('finite matrix weights and positive inter-epoch steps required')
    opnorm = lambda x: float(np.linalg.norm(x, 2))
    bound = amplitude * opnorm(a.sum(axis=0))
    for k,h in enumerate(dt, 1):
        bound += slow_change(amplitude, rate, float(h)) * opnorm(a[k:].sum(axis=0))
    return min(bound, amplitude * sum(opnorm(ak) for ak in a))


def audit_fast_prefix(times_s, values, amplitude: float, fast: FastWindow | None,
                      *, tolerance: float = 0.0) -> dict:
    """Audit ALL placed windows on a finite ZOH prefix, never an all-time cert.

    N values require N+1 timestamps; the final timestamp ends the last cell.
    The norm is convex on each (start,end) polygon. Its maximum occurs at
    knot pairs or at a knot intersecting end-start=H. Testing the enlarged
    knot +/- H candidate grid therefore covers all windows, not just one mean
    or windows aligned to proof-word boundaries. O(N^2), for small proof traces.
    """
    import numpy as np
    _nonnegative(amplitude, tolerance)
    t, v = np.asarray(times_s, dtype=float), np.asarray(values, dtype=float)
    if (t.ndim != 1 or len(t) < 2 or v.shape != (len(t)-1, 3)
            or not np.isfinite(t).all() or not np.isfinite(v).all() or np.any(np.diff(t) <= 0)):
        raise ValueError('finite N-by-3 values on N positive complete hold cells required')
    if len(t) > 4096:
        raise ValueError('proof-prefix audit limited to 4095 cells; no uncharged decimation is allowed')
    amp_pass = bool(np.all(np.linalg.norm(v, axis=1) <= amplitude + tolerance))
    out = {'amplitude_pass':amp_pass, 'finite_prefix_pass':None,
           'temporal_qualification':'OPEN', 'all_time_certified':False,
           'max_short_window_integral':None, 'interval_s':[float(t[0]),float(t[-1])]}
    if fast is None:
        return out
    fast.validate_amplitude(amplitude)
    grid = np.unique(np.concatenate((t,t+fast.horizon_s,t-fast.horizon_s)))
    grid = grid[(grid >= t[0]) & (grid <= t[-1])]
    u = np.vstack((np.zeros(3), np.cumsum(np.diff(t)[:,None]*v, axis=0)))
    vals = np.column_stack([np.interp(grid,t,u[:,i]) for i in range(3)])
    largest = 0.0
    for i,s in enumerate(grid):
        # Include the horizon boundary despite a one-ulp subtraction error.
        stop = np.searchsorted(grid, np.nextafter(s+fast.horizon_s, np.inf), side='right')
        if stop > i+1:
            largest = max(largest, float(np.linalg.norm(vals[i+1:stop]-vals[i],axis=1).max()))
    out.update(finite_prefix_pass=amp_pass and largest <= fast.cap+tolerance,
               temporal_qualification='FINITE_PREFIX_ONLY', max_short_window_integral=largest)
    return out
