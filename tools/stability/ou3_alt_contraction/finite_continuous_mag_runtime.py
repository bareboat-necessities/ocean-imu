"""Conditional finite-real continuous magnetic calibration and coupled apply.

Supplying source: ContinuousMagHardIronEstimator.h and the outer OU-III wrapper.
The exact 3-by-3 relations retain one set of sufficient statistics, including
solve failures and the wrapper mutations preceding an unsuccessful application.
No numerical estimate, field reference or gain may be supplied independently.

This evaluator is rational. Exponential values and eigensolver outcomes are
explicit, same-operand witnesses, NOT deployment/libm certificates. Squared
norm predicates are equivalent to the shipping predicates over finite reals;
float/double casts and nonfinite branches remain open. The source-uniform bound
in the supplying note uses convexity, not the regression sample values.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction as F

ZERO = (F(0),) * 3
IDENTITY = tuple(tuple(F(i == j) for j in range(3)) for i in range(3))
ZERO_MATRIX = (ZERO,) * 3


def rational(x):
    if isinstance(x, (float, bool)):
        raise TypeError('use an explicit rational, not an implicit floating value')
    return F(x)


def vector(v):
    v = tuple(rational(x) for x in v)
    if len(v) != 3:
        raise ValueError('three coordinates required')
    return v


def matrix(a):
    a = tuple(vector(row) for row in a)
    if len(a) != 3:
        raise ValueError('three rows required')
    return a


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def scale(a, k):
    return tuple(k*x for x in a)


def transpose(a):
    return tuple(zip(*a))


def mv(a, v):
    return tuple(dot(row, v) for row in a)


def mm(a, b):
    return tuple(tuple(dot(row, col) for col in transpose(b)) for row in a)


def solve(a, rhs):
    """Exact equation solve, not an emulation of Eigen's rounded LDLT."""
    rows = [list(row)+[v] for row, v in zip(matrix(a), vector(rhs))]
    for j in range(3):
        pivot = next((i for i in range(j, 3) if rows[i][j]), None)
        if pivot is None:
            raise NotImplementedError('singular/nonfinite deployed LDLT branch remains open')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        p = rows[j][j]
        rows[j] = [x/p for x in rows[j]]
        for i in range(3):
            if i != j:
                k = rows[i][j]
                rows[i] = [x-k*y for x, y in zip(rows[i], rows[j])]
    return tuple(row[3] for row in rows)


def require_rotation(a):
    a = matrix(a)
    if mm(transpose(a), a) != IDENTITY:
        raise ValueError('same normalized proxy rotation must be orthogonal')
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
           - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    if det != 1:
        raise ValueError('proper proxy rotation required')
    return a


@dataclass(frozen=True)
class Decay:
    dt: F
    tau: F
    value: F

    def __post_init__(self):
        for name in ('dt', 'tau', 'value'):
            object.__setattr__(self, name, rational(getattr(self, name)))
        if self.dt <= 0 or self.tau <= 0 or not 0 <= self.value <= 1:
            raise ValueError('positive exp arguments and a nonexpansive decay required')

    def at(self, dt, tau):
        if (self.dt, self.tau) != (dt, tau):
            raise ValueError('exponential witness detached from the same clock/configuration')
        return self.value


@dataclass(frozen=True)
class Spectrum:
    vectors: tuple
    values: tuple

    def __post_init__(self):
        u, d = matrix(self.vectors), vector(self.values)
        if mm(transpose(u), u) != IDENTITY or tuple(sorted(d)) != d:
            raise ValueError('orthonormal ordered eigenpair witness required')
        object.__setattr__(self, 'vectors', u)
        object.__setattr__(self, 'values', d)

    def for_matrix(self, a):
        ud = tuple(tuple(self.vectors[i][j]*self.values[j] for j in range(3))
                   for i in range(3))
        if mm(ud, transpose(self.vectors)) != a:
            raise ValueError('eigenpair witness detached from the same sufficient statistics')
        return self.values[0]


@dataclass(frozen=True)
class Config:
    # Wrapper defaults, not the estimator class's different standalone defaults.
    memory: F = F(600)
    solve_period: F = F(1)
    min_weight: F = F(500)
    min_information: F = F(1, 10)
    ridge: F = F(1, 2000)
    relative_ridge: F = F(1, 4)
    max_bias_fraction: F = F(7, 20)
    max_residual: F = F(3)
    min_norm: F = F(1, 1000)

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, rational(getattr(self, name)))
        if self.min_weight < 0 or self.min_norm < 0 or self.max_bias_fraction < 0:
            raise ValueError('nonnegative physical calibration gates required')


@dataclass(frozen=True)
class Estimate:
    bias: tuple = ZERO
    field: tuple = ZERO
    information: F = F(0)
    effective_weight: F = F(0)
    residual_squared: F | None = None
    valid: bool = False

    def __post_init__(self):
        object.__setattr__(self, 'bias', vector(self.bias))
        object.__setattr__(self, 'field', vector(self.field))
        for name in ('information', 'effective_weight'):
            object.__setattr__(self, name, rational(getattr(self, name)))
        if self.residual_squared is not None:
            x = rational(self.residual_squared)
            if x < 0:
                raise ValueError('nonnegative squared residual required')
            object.__setattr__(self, 'residual_squared', x)
        if not isinstance(self.valid, bool):
            raise TypeError('literal estimate-valid branch required')


@dataclass(frozen=True)
class State:
    weight: F = F(0)
    rotation_sum: tuple = ZERO_MATRIX
    level_sum: tuple = ZERO
    body_sum: tuple = ZERO
    square_sum: F = F(0)
    elapsed: F = F(0)
    estimate: Estimate = Estimate()

    def __post_init__(self):
        for name in ('weight', 'square_sum', 'elapsed'):
            x = rational(getattr(self, name))
            if x < 0:
                raise ValueError('nonnegative sufficient statistic required')
            object.__setattr__(self, name, x)
        object.__setattr__(self, 'rotation_sum', matrix(self.rotation_sum))
        object.__setattr__(self, 'level_sum', vector(self.level_sum))
        object.__setattr__(self, 'body_sum', vector(self.body_sum))
        if not isinstance(self.estimate, Estimate):
            raise TypeError('persistent calibration estimate required')


@dataclass(frozen=True)
class Result:
    state: State
    sample_accepted: bool
    solve_due: bool
    branch: str


def _fit(s, cfg, *, eigen_success, spectrum):
    out = Estimate(effective_weight=s.weight)
    if s.weight < cfg.min_weight:
        if eigen_success is not None or spectrum is not None:
            raise ValueError('underweight solve consumes no eigensolver witness')
        return out, 'weight'
    if not isinstance(eigen_success, bool):
        raise TypeError('literal same-matrix eigensolver outcome required')
    if not eigen_success:
        if spectrum is not None:
            raise ValueError('failed eigensolver consumes no eigenpairs')
        return out, 'eigensolver_rejected'
    if not isinstance(spectrum, Spectrum):
        raise TypeError('successful eigensolve requires same-matrix eigenpairs')
    a = tuple(scale(row, 1/s.weight) for row in s.rotation_sum)
    w, m = scale(s.level_sum, 1/s.weight), scale(s.body_sum, 1/s.weight)
    ata = mm(transpose(a), a)
    normal = tuple(tuple(IDENTITY[i][j]-ata[i][j] for j in range(3)) for i in range(3))
    info = s.weight*spectrum.for_matrix(normal)
    out = replace(out, information=info)
    if info < cfg.min_information:
        return out, 'information'
    ridge = max(F(0), cfg.ridge) + max(F(0), cfg.relative_ridge)*sum(normal[i][i] for i in range(3))/3
    regularized = tuple(tuple(normal[i][j]+(ridge if i == j else F(0))
                              for j in range(3)) for i in range(3))
    rhs = sub(m, mv(transpose(a), w))
    b = solve(regularized, rhs)
    field = sub(w, mv(a, b))
    if s.square_sum/s.weight <= cfg.min_norm**2:
        return out, 'field_norm'
    if dot(b, b) > cfg.max_bias_fraction**2*s.square_sum/s.weight:
        return out, 'bias_norm'
    sse = (s.square_sum - 2*dot(field, s.level_sum) - 2*dot(b, s.body_sum)
           + 2*dot(field, mv(s.rotation_sum, b))
           + s.weight*(dot(field, field)+dot(b, b)))
    out = replace(out, residual_squared=max(F(0), sse/s.weight))
    if cfg.max_residual > 0 and out.residual_squared > cfg.max_residual**2:
        return out, 'residual'
    return replace(out, bias=b, field=field, valid=True), 'accepted'


def update(s: State, cfg: Config, *, dt, rotation, raw_body,
           decay: Decay | None = None, eigen_success=None, spectrum: Spectrum | None = None):
    """Literal finite-real accumulator/solve relation after wrapper clock service."""
    if not isinstance(s, State) or not isinstance(cfg, Config):
        raise TypeError('continuous calibration state/config required')
    dt, raw = rational(dt), vector(raw_body)
    if dt <= 0 or dot(raw, raw) <= cfg.min_norm**2:
        if any(x is not None for x in (decay, eigen_success, spectrum)):
            raise ValueError('rejected sample consumes no decay or solve witnesses')
        return Result(s, False, False, 'sample_rejected')
    rot = require_rotation(rotation)
    if cfg.memory > 0:
        if not isinstance(decay, Decay):
            raise TypeError('same-dt/memory exponential witness required')
        lam = decay.at(dt, cfg.memory)
    else:
        if decay is not None:
            raise ValueError('cumulative-memory branch consumes no exponential witness')
        lam = F(1)
    nxt = State(lam*s.weight+1,
                tuple(add(scale(old, lam), new) for old, new in zip(s.rotation_sum, rot)),
                add(scale(s.level_sum, lam), mv(rot, raw)),
                add(scale(s.body_sum, lam), raw),
                lam*s.square_sum+dot(raw, raw), s.elapsed+dt, s.estimate)
    due = nxt.elapsed >= max(F(1, 50), cfg.solve_period)
    if not due:
        if eigen_success is not None or spectrum is not None:
            raise ValueError('not-due solve consumes no eigensolver operands')
        return Result(nxt, True, False, 'not_due')
    estimate, branch = _fit(nxt, cfg, eigen_success=eigen_success, spectrum=spectrum)
    return Result(replace(nxt, elapsed=F(0), estimate=estimate), True, True, branch)


def level_reference(s: State, bias, min_norm):
    """Same-statistics reference, including the shipping nonzero-norm guard."""
    b = vector(bias)
    if s.weight <= F(1, 10**6):
        return None
    ref = scale(sub(s.level_sum, mv(s.rotation_sum, b)), 1/s.weight)
    return ref if dot(ref, ref) > rational(min_norm)**2 else None


@dataclass(frozen=True)
class HorizontalNorm:
    xy: tuple
    value: F

    def __post_init__(self):
        xy = tuple(rational(x) for x in self.xy)
        v = rational(self.value)
        if len(xy) != 2 or v < 0 or v*v != dot(xy, xy):
            raise ValueError('same-input exact horizontal norm required')
        object.__setattr__(self, 'xy', xy)
        object.__setattr__(self, 'value', v)


@dataclass(frozen=True)
class Applied:
    applied: tuple = ZERO
    startup_bias: tuple = ZERO
    anchor_bias: tuple | None = None
    anchor_reference: tuple | None = None
    last_time: F | None = None

    def __post_init__(self):
        for name in ('applied', 'startup_bias'):
            object.__setattr__(self, name, vector(getattr(self, name)))
        if (self.anchor_bias is None) != (self.anchor_reference is None):
            raise ValueError('calibration anchor fields latch together')
        for name in ('anchor_bias', 'anchor_reference'):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, vector(getattr(self, name)))
        if self.last_time is not None:
            object.__setattr__(self, 'last_time', rational(self.last_time))

    @property
    def total_bias(self):
        return add(self.startup_bias, self.applied)


@dataclass(frozen=True)
class ApplyResult:
    state: Applied
    reference: tuple
    wrote_reference: bool
    branch: str


def apply(s: Applied, stats: State, reference, *, time, sample_dt=F(1, 200),
          fraction=F(1), slew_tau=F(45), enabled=True, live=True,
          reference_set=True, reference_valid=True, refinement_enabled=True,
          refinement_done=False, proxy_min_norm=F(1, 1000),
          decay: Decay | None = None, new_norm: HorizontalNorm | None = None,
          anchor_norm: HorizontalNorm | None = None):
    """Coupled offset/reference write AFTER possible refinement, BEFORE measurement."""
    if not isinstance(s, Applied) or not isinstance(stats, State):
        raise TypeError('persistent applied calibration and sufficient statistics required')
    ref = vector(reference)
    flags = (enabled, live, reference_set, reference_valid, refinement_enabled, refinement_done)
    if any(not isinstance(x, bool) for x in flags):
        raise TypeError('literal wrapper branches required')
    blocked = (not enabled or not live or not reference_set or not reference_valid
               or (refinement_enabled and not refinement_done) or not stats.estimate.valid)
    if blocked:
        if any(x is not None for x in (decay, new_norm, anchor_norm)):
            raise ValueError('blocked apply consumes no arithmetic witnesses')
        return ApplyResult(s, ref, False, 'blocked')
    time, sample_dt, fraction, tau = map(rational, (time, sample_dt, fraction, slew_tau))
    if time < 0 or sample_dt <= 0:
        raise ValueError('nonnegative clock and positive fallback dt required')
    if s.anchor_bias is None:
        s = replace(s, anchor_bias=s.applied, anchor_reference=ref)
    dt = time-s.last_time if s.last_time is not None and time > s.last_time else sample_dt
    s = replace(s, last_time=time)  # persists even if either subsequent reference gate fails
    if tau > F(1, 1000):
        if not isinstance(decay, Decay):
            raise TypeError('same-apply-dt/slew exponential witness required')
        alpha = 1-decay.at(dt, tau)
    else:
        if decay is not None:
            raise ValueError('instantaneous slew consumes no exponential witness')
        alpha = F(1)
    target = scale(stats.estimate.bias, fraction)
    candidate = add(s.applied, scale(sub(target, s.applied), alpha))
    ln = level_reference(stats, candidate, proxy_min_norm)
    la = level_reference(stats, s.anchor_bias, proxy_min_norm)
    if ln is None or la is None:
        if new_norm is not None or anchor_norm is not None:
            raise ValueError('invalid level reference consumes no horizontal norms')
        return ApplyResult(s, ref, False, 'level_reference')
    if not isinstance(new_norm, HorizontalNorm) or new_norm.xy != ln[:2]:
        raise ValueError('new horizontal norm detached from current calibration statistics')
    if not isinstance(anchor_norm, HorizontalNorm) or anchor_norm.xy != la[:2]:
        raise ValueError('anchor horizontal norm detached from current calibration statistics')
    h = s.anchor_reference[0]+new_norm.value-anchor_norm.value
    z = s.anchor_reference[2]+ln[2]-la[2]
    if h <= rational(proxy_min_norm):
        return ApplyResult(s, ref, False, 'horizontal_reference')
    return ApplyResult(replace(s, applied=candidate), (h, F(0), z), True, 'accepted')


def deterministic_bounds(*, world=F(75), hard_iron=F(5), residual=F(2), fraction=F(1)):
    """Consequences of the documented source-bound induction, not source admission.

    Default startup offset is zero and the continuous anchor is taken before
    any continuous application. A different initializer is not covered by these
    constants; arbitrary source/moment snapshots do not satisfy the induction.
    """
    world, hard_iron, residual, fraction = map(rational, (world, hard_iron, residual, fraction))
    if min(world, hard_iron, residual, fraction) < 0:
        raise ValueError('nonnegative deterministic envelopes required')
    raw = world+hard_iron+residual
    fitted = F(7, 20)*raw
    applied = fraction*fitted
    active_reference = raw+applied
    return {
        'raw_norm_uT': raw,
        'accepted_fitted_bias_norm_uT': fitted,
        'applied_bias_norm_uT': applied,
        'active_reference_norm_uT': active_reference,
        'corrected_measurement_norm_uT': raw+applied,
        'physical_to_active_residual_norm_uT': raw+active_reference+applied,
    }


def readiness():
    return {
        'continuous_raw_proxy_statistics_and_solve_order_materialized': True,
        'fit_not_free_but_solved_from_same_statistics': True,
        'rejected_solve_replaces_previously_valid_estimate': True,
        'applied_offset_and_reference_use_same_statistics': True,
        'anchor_and_apply_clock_survive_later_rejection': True,
        'default_real_calibration_bound_induction_available': True,
        'same_history_startup_statistics_qualified': False,
        'deployment_exp_cast_solver_nonfinite_branches_closed': False,
        'complete_word_finite_identity': False,
        'ALT_LIVE_PASS': False,
        'ALT_STARTUP_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
