"""Exact-shadow versus machine-rooted applied ActiveParameters displacement.

The tuner candidate state is not the state used by prediction immediately.
Shipping installs tau, stationary Sigma_aw, pseudo-S cadence and Live R_S only
at goLive or a later pending boundary.  The exact finite Live word currently
uses the exact-real commit shadow, while deployment owns a global compiler-mode
machine TuneState.  This module keeps their applied parameter discrepancy as an
explicit same-boundary supply instead of identifying them.

No smallness bound is asserted here.  In particular, scheduler timing effects of
``delta pseudo_period`` and coefficient effects of ``delta tau/Sigma/R_S`` must
be propagated through the finite Live event relation before storage search.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_runtime_parameters as A

QUALIFICATION='OU3_ALT_ACTIVE_PARAMETER_MACHINE_REAL_JOIN_V1'


def _submat(a,b):
    return tuple(tuple(F(a[i][j])-F(b[i][j]) for j in range(3)) for i in range(3))


@dataclass(frozen=True)
class Supply:
    tau:F
    Sigma_aw:tuple
    pseudo_period:F
    R_S:tuple|None
    def __post_init__(self):
        object.__setattr__(self,'tau',F(self.tau)); object.__setattr__(self,'pseudo_period',F(self.pseudo_period))
        s=tuple(tuple(F(x) for x in row) for row in self.Sigma_aw)
        if len(s)!=3 or any(len(row)!=3 for row in s): raise ValueError('3x3 Sigma_aw supply required')
        object.__setattr__(self,'Sigma_aw',s)
        if self.R_S is not None:
            r=tuple(tuple(F(x) for x in row) for row in self.R_S)
            if len(r)!=3 or any(len(row)!=3 for row in r): raise ValueError('3x3 R_S supply required')
            object.__setattr__(self,'R_S',r)


@dataclass(frozen=True)
class Join:
    exact:A.ActiveParameters
    machine:A.ActiveParameters
    supply:Supply
    mode:str
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.exact,A.ActiveParameters) or not isinstance(self.machine,A.ActiveParameters):
            raise TypeError('exact and machine ActiveParameters required')
        if self.mode not in ('separate','fma'): raise ValueError('invalid global compiler mode')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong active-parameter join qualification')
        if self.supply.tau!=self.machine.tau-self.exact.tau:
            raise ValueError('tau active-parameter supply detached')
        if self.supply.Sigma_aw!=_submat(self.machine.Sigma_aw,self.exact.Sigma_aw):
            raise ValueError('Sigma_aw active-parameter supply detached')
        if self.supply.pseudo_period!=self.machine.pseudo_period-self.exact.pseudo_period:
            raise ValueError('pseudo-period active-parameter supply detached')
        if (self.machine.R_S is None)!=(self.exact.R_S is None):
            raise ValueError('machine/exact R_S activation branches differ')
        expected=None if self.machine.R_S is None else _submat(self.machine.R_S,self.exact.R_S)
        if self.supply.R_S!=expected: raise ValueError('R_S active-parameter supply detached')


def join(exact:A.ActiveParameters,machine:A.ActiveParameters,mode:str):
    if not isinstance(exact,A.ActiveParameters) or not isinstance(machine,A.ActiveParameters):
        raise TypeError('exact and machine ActiveParameters required')
    rs=None if machine.R_S is None or exact.R_S is None else _submat(machine.R_S,exact.R_S)
    return Join(exact,machine,Supply(machine.tau-exact.tau,
        _submat(machine.Sigma_aw,exact.Sigma_aw),machine.pseudo_period-exact.pseudo_period,rs),mode)


def readiness():
    return {
      'same_boundary_exact_and_machine_applied_parameters_joined':True,
      'tau_stationary_Sigma_pseudo_period_and_RS_displacements_exposed':True,
      'candidate_state_not_confused_with_applied_active_state':True,
      'active_parameter_displacement_bound_closed':False,
      'pseudo_period_displacement_scheduler_effect_closed':False,
      'tau_Sigma_displacement_prediction_effect_closed':False,
      'RS_displacement_measurement_effect_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
