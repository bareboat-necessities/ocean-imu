"""Carry exact-candidate -> machine-input displacement into SpectralMSE R_S.

The local SpectralMSE deployment theorem is intentionally rooted at binary32
``tau`` and ``sigma`` operands.  The physical/frontend candidate is exact-real.
Those are not the same object.  This layer composes the already source-owned tau
step and sigma machine/real join into one binary32 ``TargetState`` for the R_S
machine graph while retaining every input displacement from the exact candidate.

The resulting local spectral join then measures machine arithmetic/libm error
relative to the exact mathematical SpectralMSE relation at those SAME machine
inputs.  Therefore two effects remain separate:

  (1) exact physical candidate -> binary32 tau/sigma/input displacement; and
  (2) binary32-input exact spectral relation -> deployed pow/sqrt/cache graph.

No Lipschitz/smallness bound between the two spectral targets is asserted here.
That source-uniform propagation remains a later finite-master obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as SPEC

QUALIFICATION='OU3_ALT_SPECTRAL_INPUT_JOIN_V1'


@dataclass(frozen=True)
class InputSupply:
    frequency:F
    variance_wave:F
    tau_target:F
    sigma_target:F
    def __post_init__(self):
        for n in ('frequency','variance_wave','tau_target','sigma_target'):
            object.__setattr__(self,n,F(getattr(self,n)))


@dataclass(frozen=True)
class Join:
    exact_target:C.TargetState
    machine_target:C.TargetState
    tau_step:TAU.TauStep
    sigma_join:SIG.Join
    input_supply:InputSupply
    spectral:SPEC.Join
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.exact_target,C.TargetState) or not isinstance(self.machine_target,C.TargetState):
            raise TypeError('exact and machine-input TargetState required')
        if not isinstance(self.tau_step,TAU.TauStep) or not isinstance(self.sigma_join,SIG.Join):
            raise TypeError('source-owned tau step and sigma join required')
        if not isinstance(self.input_supply,InputSupply) or not isinstance(self.spectral,SPEC.Join):
            raise TypeError('input supply and local spectral join required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong spectral input join qualification')
        if self.exact_target!=self.sigma_join.exact_target:
            raise ValueError('spectral exact target detached from sigma exact target')
        expected=C.TargetState(self.tau_step.frequency,self.sigma_join.machine.var_wave,
                               self.tau_step.tau_target,self.sigma_join.machine.sigma_target)
        if self.machine_target!=expected:
            raise ValueError('spectral machine TargetState detached from tau/sigma machine operands')
        if self.spectral.target!=self.machine_target:
            raise ValueError('local spectral graph detached from machine input target')
        s=self.input_supply
        if s.frequency!=self.machine_target.frequency-self.exact_target.frequency:
            raise ValueError('spectral frequency input supply detached')
        if s.variance_wave!=self.machine_target.variance_wave-self.exact_target.variance_wave:
            raise ValueError('spectral variance input supply detached')
        if s.tau_target!=self.machine_target.tau_target-self.exact_target.tau_target:
            raise ValueError('spectral tau input supply detached')
        if s.sigma_target!=self.machine_target.sigma_target-self.exact_target.sigma_target:
            raise ValueError('spectral sigma input supply detached')
        if s.variance_wave!=self.sigma_join.supply.variance_wave or s.sigma_target!=self.sigma_join.supply.sigma_target:
            raise ValueError('spectral sigma-side inputs do not reuse same sigma join supplies')


def join(tau_step:TAU.TauStep,sigma_join:SIG.Join,*,pow_result,sqrt_result,bits=96):
    if not isinstance(tau_step,TAU.TauStep) or not isinstance(sigma_join,SIG.Join):
        raise TypeError('TauStep and same-source sigma join required')
    exact=sigma_join.exact_target
    machine=C.TargetState(tau_step.frequency,sigma_join.machine.var_wave,
                          tau_step.tau_target,sigma_join.machine.sigma_target)
    supply=InputSupply(machine.frequency-exact.frequency,
                       machine.variance_wave-exact.variance_wave,
                       machine.tau_target-exact.tau_target,
                       machine.sigma_target-exact.sigma_target)
    spectral=SPEC.join(sigma_join.deployment_cfg,machine,
                       pow_result=pow_result,sqrt_result=sqrt_result,bits=bits)
    return Join(exact,machine,tau_step,sigma_join,supply,spectral)


def readiness():
    return {
      'source_owned_tau_step_and_sigma_machine_real_join_composed':True,
      'machine_SpectralMSE_target_built_only_from_same_tau_sigma_machine_operands':True,
      'exact_candidate_to_machine_frequency_variance_tau_sigma_supplies_exposed':True,
      'sigma_variance_and_target_supplies_reused_without_duplication':True,
      'machine_input_displacement_kept_separate_from_local_spectral_libm_roundoff':True,
      'spectral_input_displacement_to_exact_candidate_RS_bound_closed':False,
      'target_libm_correspondence_closed':False,
      'source_uniform_spectral_input_supply_bounds_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
