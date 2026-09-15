"""Literal finite startup decision; its input predicates still need source proofs.

In particular the timeout requires initialized proxy + aligned gravity, but
neither TunerReady nor magnetic north. The deadline arithmetic below retains
the binary32 evaluation order of maybeHandOffToMekf_. A decision is a conditional
program relation, not evidence that a physical history reaches its inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK


@dataclass(frozen=True)
class Config:
    with_mag: bool = True
    minimum: F = F(8)
    timeout: F = F(150)
    settle: F = F(0)
    magnetic_window: F = F(15)
    fallback: F = F(30)
    gravity_hold: F = F(2)

    def __post_init__(self):
        if type(self.with_mag) is not bool:
            raise TypeError('literal with_mag flag required')
        for name in ('minimum', 'timeout', 'settle', 'magnetic_window', 'fallback', 'gravity_hold'):
            value = F(getattr(self, name))
            if value < 0 or not B.is_binary32(value):
                raise ValueError('finite nonnegative binary32 startup configuration required')
            object.__setattr__(self, name, value)

    @property
    def acquisition_deadline(self):
        return (B.rn32(B.rn32(self.settle + B.rn32(2 * max(self.magnetic_window, F(1))))
                        + self.fallback) if self.with_mag else F(0))

    @property
    def deadline(self):
        return max(self.timeout, self.acquisition_deadline)


@dataclass(frozen=True)
class Decision:
    config: Config
    wrapper_time: F
    begun: bool
    proxy_initialized: bool
    tuner_stage: str
    aligned_branch: bool
    gravity_good: F
    north_reference: bool

    def __post_init__(self):
        if not isinstance(self.config, Config):
            raise TypeError('carried startup control configuration required')
        for name in ('begun', 'proxy_initialized', 'aligned_branch', 'north_reference'):
            if type(getattr(self, name)) is not bool:
                raise TypeError('literal startup branch predicates required')
        for name in ('wrapper_time', 'gravity_good'):
            value = F(getattr(self, name))
            if value < 0 or not B.is_binary32(value):
                raise ValueError('finite nonnegative binary32 startup clock required')
            object.__setattr__(self, name, value)
        if self.tuner_stage not in ('Cold', 'TunerWarm', 'TunerReady'):
            raise ValueError('handoff starts from a pre-Live stage')

    @property
    def by_quality(self):
        return bool(self.begun and self.proxy_initialized
                    and self.wrapper_time >= self.config.minimum
                    and self.aligned_branch and self.gravity_good >= self.config.gravity_hold
                    and (not self.config.with_mag or self.north_reference)
                    and self.tuner_stage == 'TunerReady')

    @property
    def by_timeout(self):
        return bool(self.begun and self.proxy_initialized
                    and self.wrapper_time >= self.config.deadline and self.aligned_branch)

    @property
    def handoff_due(self):
        return self.by_quality or self.by_timeout


def require_bridge(decision, frontend):
    """Attach the conditional decision to the same stage/proxy/clock prefix.

    Gravity/north predicate production remains an explicit master obligation.
    Legacy TunerReady callers retain their conditional boundary implication.
    Non-TunerReady callers must provide the actual timeout branch relation.
    """
    stage = frontend.tuner.stage
    if decision is None:
        if stage != 'TunerReady':
            raise ValueError('non-TunerReady goLive requires an attached timeout decision')
        return
    if not isinstance(decision, Decision) or not decision.handoff_due:
        raise ValueError('startup handoff decision is not due')
    if decision.tuner_stage != stage or decision.proxy_initialized != frontend.tuner.vertical.initialized:
        raise ValueError('handoff decision detached from carried frontend stage/proxy')
    if decision.wrapper_time != CLOCK.clock_at_real_grid_time(frontend.tuner.time):
        raise ValueError('handoff decision detached from carried physical clock')


def readiness():
    return {
        'quality_and_timeout_predicates_materialized': True,
        'timeout_does_not_require_TunerReady_or_north': True,
        'timeout_retains_initialized_proxy_and_aligned_branch': True,
        'default_magnetic_acquisition_deadline_is_60_s': Config().acquisition_deadline == 60,
        'source_uniform_control_predicate_production_closed': False,
        'nonfinite_control_branches_qualified': False,
        'universal_startup_deadline_closed': False,
        'ALT_STARTUP_PASS': False,
    }
