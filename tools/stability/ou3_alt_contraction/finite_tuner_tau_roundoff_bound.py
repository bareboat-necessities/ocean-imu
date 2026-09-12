"""Analytic binary32 roundoff bound for the shipping tau EMA edge.

On the source-locked scalar domain |p_f|<=13, |t|<=12 and the canonical
5 ms update, shipping's adaptive horizon is >=0.05 s.  Hence x=dt/adapt<=0.1,
the carried exp witness has e>=1-x, and both exact and rounded EMA alphas are
<1/8.  This makes the product small enough to use much tighter operation-wise
RNE cells than a generic <32 magnitude bound.

With source target error <=2^-21:
  * target-minus-previous subtraction roundoff <=2^-21;
  * alpha subtraction roundoff <=2^-28 because alpha<1/8;
  * alpha*delta multiplication roundoff <=2^-24 because |product|<2;
  * final add/FMA roundoff <=2^-21 because the result stays below 16.

The propagated discrepancy is <2^-20 s for BOTH the FMA and separate mul/add
shapes.  We therefore publish 2^-20 s as the uniform local tau supply bound.
This bound uses no target-libm correctness: e_f remains the same explicit
binary32 exp witness already carried by the tau edge.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET

PREVIOUS_ABS_MAX=F(13); TARGET_ABS_MAX=F(12); DELTA_ABS_MAX=F(13)
TARGET_ERROR_MAX=F(1,1<<21); SUB_ROUND=F(1,1<<21)
ALPHA_ROUND=F(1,1<<28); MUL_ROUND=F(1,1<<24); FINAL_ROUND=F(1,1<<21)
ALPHA_MAX=F(1,8)
FMA_DERIVED_MAX=ALPHA_MAX*(TARGET_ERROR_MAX+SUB_ROUND)+DELTA_ABS_MAX*ALPHA_ROUND+FINAL_ROUND
SEPARATE_DERIVED_MAX=FMA_DERIVED_MAX+MUL_ROUND
UNIFORM_RESIDUAL_MAX=F(1,1<<20)


@dataclass(frozen=True)
class BoundCertificate:
    supply: TAU.TauRoundoffSupply
    target_error:F
    bound:F=UNIFORM_RESIDUAL_MAX
    def __post_init__(self):
        te=abs(F(self.target_error)); b=F(self.bound)
        object.__setattr__(self,'target_error',te); object.__setattr__(self,'bound',b)
        if b!=UNIFORM_RESIDUAL_MAX: raise ValueError('tau roundoff theorem bound is source-fixed')
        if te>TARGET_ERROR_MAX: raise ValueError('tau target binary32/exact discrepancy exceeds certified cell')
        if abs(self.supply.residual_fma)>FMA_DERIVED_MAX: raise ValueError('FMA tau residual exceeds analytic derived bound')
        if abs(self.supply.residual_separate)>SEPARATE_DERIVED_MAX: raise ValueError('separate tau residual exceeds analytic derived bound')


def _domain(binary:TAU.TauStep, exact_target):
    if not isinstance(binary,TAU.TauStep): raise TypeError('TauStep required')
    t=F(exact_target)
    if abs(binary.previous)>PREVIOUS_ABS_MAX: raise ValueError('previous tau outside certified scalar domain')
    if abs(t)>TARGET_ABS_MAX: raise ValueError('exact tau target outside certified scalar domain')
    a=F(1)-binary.exp_decay
    if not 0<=a<ALPHA_MAX or not 0<=binary.alpha<ALPHA_MAX:
        raise ValueError('tau EMA alpha outside certified canonical-5ms domain')
    return t


def certify(binary:TAU.TauStep, exact, *, exact_target):
    t=_domain(binary,exact_target); supply=TAU.roundoff_supply(binary,exact)
    if F(exact.tau_target)!=t: raise ValueError('exact target argument detached from exact candidate')
    return BoundCertificate(supply,binary.tau_target-t)


def certify_source_target(binary:TAU.TauStep,target:TARGET.TargetPair):
    if not isinstance(target,TARGET.TargetPair): raise TypeError('source-locked TargetPair required')
    t=_domain(binary,target.exact_target)
    if binary.frequency!=target.clamped_frequency: raise ValueError('binary tau frequency detached from source-locked target pair')
    if binary.tau_target!=target.binary32_target: raise ValueError('binary tau target detached from source-locked target pair')
    a=F(1)-binary.exp_decay
    real_next=binary.previous+a*(t-binary.previous)
    supply=TAU.TauRoundoffSupply(real_next,binary.next_separate,binary.next_fma,
                                  binary.next_separate-real_next,binary.next_fma-real_next)
    return BoundCertificate(supply,target.error)


def readiness():
    target=TARGET.readiness()
    return {
      'canonical_5ms_alpha_strictly_below_one_eighth':True,
      'operationwise_normal_binary32_rounding_cells_closed':True,
      'tau_FMA_roundoff_supply_uniform_bound_conditional_on_scalar_domain':True,
      'tau_separate_mul_add_roundoff_supply_uniform_bound_conditional_on_scalar_domain':True,
      'uniform_tau_roundoff_bound_seconds':UNIFORM_RESIDUAL_MAX,
      'compiler_contraction_selection_needed_for_bound':False,
      'target_libm_correctness_needed_for_roundoff_bound':False,
      'source_locked_target_pair_can_define_roundoff_supply_without_free_exact_candidate':True,
      'shipping_tau_target_exact_vs_binary32_cell_bound_source_locked':target['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'],
      'shipping_tau_predecessor_domain_inductively_closed':False,
      'source_uniform_tau_roundoff_supply_bound_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
