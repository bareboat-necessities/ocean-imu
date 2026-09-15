"""Literal post-log-update WPE usability latch on the finite binary32 graph.

The comparison reads the updated log-period, elapsed time and the same lambda.
The period exp remains a target-libm witness, including a nonfinite outcome;
this module establishes the control relation, not numerical libm accuracy.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
AUDITED_SOURCE_SHA='91c15942f24e08f1c39fa131f632d0974e6ab536a93a29a5af7650024ee41b7e'


@dataclass(frozen=True)
class PeriodWitness:
    log_period:F|None
    period_exp:F|None  # None denotes the nonfinite return class.
    def __post_init__(self):
        for name in ('log_period','period_exp'):
            x=getattr(self,name)
            if x is not None:
                x=F(x)
                if not B.is_binary32(x): raise ValueError('binary32 WPE usability operand required')
                object.__setattr__(self,name,x)
        if self.log_period is None and self.period_exp is not None:
            raise ValueError('nonfinite log must take the no-exp NAN getter branch')


@dataclass(frozen=True)
class Decision:
    before:bool
    after:bool
    branch:str
    moment_start:F|None=None
    usable_floor:F|None=None
    moment_history:F|None=None


def update(before,*,produced_period,log_period,elapsed,lambda_,witness=None):
    if type(before) is not bool or type(produced_period) is not bool:
        raise TypeError('literal WPE usability/production flags required')
    if not produced_period or before:
        if witness is not None: raise ValueError('unexecuted usability getter consumes no witness')
        return Decision(before,before,'latched' if before else 'no-valid-period')
    if not isinstance(witness,PeriodWitness):
        raise TypeError('post-update WPE usability getter witness required')
    if witness.log_period!=log_period:
        raise ValueError('usability getter detached from post-update log state')
    if witness.period_exp is None or witness.period_exp<=0:
        return Decision(False,False,'invalid-period')
    t,lam=F(elapsed),F(lambda_)
    if not B.is_binary32(t) or not B.is_binary32(lam) or lam<=0:
        raise ValueError('same finite binary32 elapsed/lambda operands required')
    start=B.div(3,lam); floor=B.div(4,lam); history=B.sub(t,start)
    ready=t>=floor and history>=witness.period_exp
    return Decision(False,ready,'takeover' if ready else 'waiting',start,floor,history)


def readiness():
    if sha256(SOURCE.read_bytes()).hexdigest()!=AUDITED_SOURCE_SHA:
        raise RuntimeError('WPE source changed; re-audit usability and arithmetic order')
    return {'post_update_binary32_usability_predicates_materialized':True,
            'usable_latch_preserved_without_reexecuting_getter':True,
            'target_period_exp_correspondence_closed':False,
            'source_uniform_takeover_deadline_closed':False}
