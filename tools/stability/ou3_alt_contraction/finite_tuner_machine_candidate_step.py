"""One same-sample coherent machine TuneState candidate successor for ALT.

This is the product transition missing between the exact frontend candidate and
the common next-IMU pending boundary.  It does not independently accept tau,
sigma and R_S successors.  Instead, for each global compiler history it composes

  source-owned TauStep
    -> qualified common tau/sigma alpha
    -> same-source binary32 sigma target
    -> sigma_applied successor
    -> SpectralMSE machine-input join
    -> same-tau R_S alpha
    -> RS_applied successor

and only then forms one ``finite_tuner_machine_tunestate_product.State`` with the
exact shipping candidate's pending bit.

The exact candidate remains the common physical/frontend shadow.  All machine
input displacements and local libm supplies stay visible in the supplied joins;
this module does not bound or discard them.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as COMMON
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SIGJOIN
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_input_join as SPECIN
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as RSALPHA
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_candidate_uniform_bounds as SUPPLY

QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_CANDIDATE_STEP_V1'


@dataclass(frozen=True)
class ModeResult:
    sigma_join:SIGJOIN.Join
    common_alpha:COMMON.Qualified
    spectral_input:SPECIN.Join
    rs_alpha:RSALPHA.Join
    mode:str
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('invalid compiler mode')
        if not isinstance(self.sigma_join,SIGJOIN.Join) or not isinstance(self.common_alpha,COMMON.Qualified):
            raise TypeError('sigma join and common alpha required')
        if not isinstance(self.spectral_input,SPECIN.Join) or not isinstance(self.rs_alpha,RSALPHA.Join):
            raise TypeError('spectral-input and R_S-alpha joins required')
        tau=self.common_alpha.step
        if self.spectral_input.tau_step!=tau:
            raise ValueError('mode spectral input detached from same common-alpha tau step')
        if self.spectral_input.sigma_join!=self.sigma_join:
            raise ValueError('mode spectral input detached from same sigma join')
        if self.rs_alpha.cfg!=self.sigma_join.deployment_cfg:
            raise ValueError('mode R_S alpha detached from deployment config')
        if self.rs_alpha.tau_target!=tau.tau_target:
            raise ValueError('mode R_S alpha detached from same tau target')


@dataclass(frozen=True)
class Result:
    product:PRODUCT.Result
    sigma:SIG.Result
    rs:RS.Result
    separate:ModeResult
    fma:ModeResult
    exact_candidate:C.CandidateResult
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.product,PRODUCT.Result) or not isinstance(self.sigma,SIG.Result) or not isinstance(self.rs,RS.Result):
            raise TypeError('whole product, sigma and R_S results required')
        if not isinstance(self.separate,ModeResult) or not isinstance(self.fma,ModeResult):
            raise TypeError('both compiler-mode result ledgers required')
        if not isinstance(self.exact_candidate,C.CandidateResult): raise TypeError('exact candidate required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine candidate qualification')
        if self.separate.mode!='separate' or self.fma.mode!='fma': raise ValueError('compiler histories crossed')
        if self.product.sigma!=self.sigma or self.product.rs!=self.rs:
            raise ValueError('whole TuneState product detached from sigma/R_S results')
        if self.product.state.pending!=self.exact_candidate.pending_after:
            raise ValueError('machine pending bit detached from exact shipping candidate control branch')


def _exact_target(candidate:C.CandidateResult):
    return C.TargetState(candidate.frequency,candidate.variance_wave,
                         candidate.tau_target,candidate.sigma_target)


def step(previous:PRODUCT.State,exact_candidate:C.CandidateResult,
         candidate_cfg:C.CandidateConfig,deployment_cfg:D.DeploymentConfig,
         tau:TAU.StepResult,*,dt,
         separate_sigma_join:SIGJOIN.Join,fma_sigma_join:SIGJOIN.Join,
         separate_spectral_pow,separate_spectral_sqrt,
         fma_spectral_pow,fma_spectral_sqrt,
         separate_rs_exp_decay,fma_rs_exp_decay,uniform_supplies=False):
    if not isinstance(previous,PRODUCT.State) or not isinstance(exact_candidate,C.CandidateResult):
        raise TypeError('machine TuneState predecessor and exact candidate required')
    if not isinstance(tau,TAU.StepResult): raise TypeError('same-sample dual tau StepResult required')
    if tau.separate_step.previous!=previous.tau.separate or tau.fma_step.previous!=previous.tau.fma:
        raise ValueError('dual tau result detached from whole machine TuneState predecessor')
    if type(uniform_supplies) is not bool: raise TypeError('literal candidate supply qualification required')
    if uniform_supplies:
        SUPPLY.require_state(previous); SUPPLY.require_config(deployment_cfg)
        if F(dt)!=SUPPLY.W.DT: raise ValueError('candidate finite supplies require canonical machine dt')
        for t,sj,pw,sw in ((tau.separate_step,separate_sigma_join,separate_spectral_pow,separate_spectral_sqrt),
                          (tau.fma_step,fma_sigma_join,fma_spectral_pow,fma_spectral_sqrt)):
            if t.exp_profile!=SUPPLY.LIBM.EXP_PROFILE:
                raise ValueError('bounded candidate requires target common-alpha profile')
            SUPPLY.require_inputs(deployment_cfg,tau=t.tau_target,sigma=sj.machine.sigma_target,
                                  pow_result=pw,sqrt_result=sw)
    expected=_exact_target(exact_candidate)
    for name,j in (('separate',separate_sigma_join),('fma',fma_sigma_join)):
        if not isinstance(j,SIGJOIN.Join): raise TypeError(name+' same-source sigma join required')
        if j.exact_target!=expected:
            raise ValueError(name+' sigma join detached from common exact candidate')
        if j.candidate_cfg!=candidate_cfg or j.deployment_cfg!=deployment_cfg:
            raise ValueError(name+' sigma join detached from declared configs')

    h=F(dt)
    sa=COMMON.qualify(tau.separate_step,candidate_cfg,deployment_cfg,dt=h)
    fa=COMMON.qualify(tau.fma_step,candidate_cfg,deployment_cfg,dt=h)
    sigma=SIG.step(previous.sigma,
        separate_target=separate_sigma_join.machine,separate_alpha=sa,
        fma_target=fma_sigma_join.machine,fma_alpha=fa)

    ss=SPECIN.join(tau.separate_step,separate_sigma_join,
                   pow_result=separate_spectral_pow,sqrt_result=separate_spectral_sqrt)
    fs=SPECIN.join(tau.fma_step,fma_sigma_join,
                   pow_result=fma_spectral_pow,sqrt_result=fma_spectral_sqrt)
    sra=RSALPHA.join(deployment_cfg,tau_target=tau.separate_step.tau_target,dt=h,
                     exp_decay=separate_rs_exp_decay,
                     exp_profile=SUPPLY.LIBM.EXP_PROFILE if uniform_supplies else 'legacy-enclosure')
    fra=RSALPHA.join(deployment_cfg,tau_target=tau.fma_step.tau_target,dt=h,
                     exp_decay=fma_rs_exp_decay,
                     exp_profile=SUPPLY.LIBM.EXP_PROFILE if uniform_supplies else 'legacy-enclosure')
    rs=RS.step(previous.rs,separate_target=ss.spectral,separate_alpha=sra,
               fma_target=fs.spectral,fma_alpha=fra)

    whole=PRODUCT.compose_after_sample(previous,tau=tau,sigma=sigma,rs=rs,
                                       pending_after=exact_candidate.pending_after)
    sm=ModeResult(separate_sigma_join,sa,ss,sra,'separate')
    fm=ModeResult(fma_sigma_join,fa,fs,fra,'fma')
    result=Result(whole,sigma,rs,sm,fm,exact_candidate)
    if uniform_supplies: SUPPLY.require_result(previous,result)
    return result


def readiness():
    return {
      'one_exact_frontend_candidate_anchors_both_global_compiler_histories':True,
      'each_mode_tau_step_feeds_its_same_common_tau_sigma_alpha':True,
      'each_mode_sigma_machine_join_feeds_its_sigma_applied_successor':True,
      'each_mode_same_tau_sigma_machine_inputs_feed_SpectralMSE_target':True,
      'each_mode_RS_alpha_is_rooted_at_its_same_tau_target_and_config':True,
      'tau_sigma_RS_scalar_ledgers_advance_once_before_pending_bit_is_attached':True,
      'one_pending_bit_equals_literal_exact_candidate_control_branch':True,
      'machine_input_and_libm_supplies_retained_not_collapsed':True,
      'upstream_frontend_binary32_correspondence_closed':False,
      'all_target_libm_correspondence_closed':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
