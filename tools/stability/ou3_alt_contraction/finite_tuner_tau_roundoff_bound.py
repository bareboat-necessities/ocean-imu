"""Analytic binary32 roundoff bound for the shipping tau EMA edge.

This is deliberately a *conditional* deployment theorem. It proves a single
history-independent numerical bound once the tau edge is known to be in the
shipping scalar domain; it does not by itself prove that every startup/Live
prefix reaches and stays in that domain.

For the deployed shape

    d_f = RN32(t_f - p_f)
    a_f = RN32(1 - e_f)
    y_sep = RN32(p_f + RN32(a_f*d_f))
    y_fma = RN32(p_f + a_f*d_f)

compare against the canonical exact shadow

    y = p_f + (1-e_f) (t - p_f).

On the admitted domain |p_f| <= 13, |t| <= 12, 0 <= a_f <= 1, all ordinary
intermediates have magnitude < 32. A deliberately coarse absolute allowance
R=2^-19 therefore dominates one RNE binary32 rounding cell throughout this
normal range. If |t_f-t| <= R, then

  |d_f-(t-p_f)| <= 2R,
  |a_f-(1-e_f)| <= R,
  |a_f*d_f-(1-e_f)(t-p_f)| <= 15R,

and the final FMA is <=16R while separate mul/add is <=17R. We publish the
round number 32R = 2^-14 s, leaving almost 2x algebraic margin.

``certify_source_target`` is the theorem-facing path: its exact target comes
from ``finite_shipping_tau_target_binary32.TargetPair``, so no arbitrary exact
CandidateResult is needed to define the deployment residual.

No target-libm correctness is used: e_f is the same explicit exp witness already
carried by the tau edge. Compiler contraction need not be selected because the
bound covers both represented evaluation shapes.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET

ROUND_CELL = F(1, 1 << 19)
PREVIOUS_ABS_MAX = F(13)
TARGET_ABS_MAX = F(12)
TARGET_ERROR_MAX = ROUND_CELL
FMA_DERIVED_MAX = 16 * ROUND_CELL
SEPARATE_DERIVED_MAX = 17 * ROUND_CELL
UNIFORM_RESIDUAL_MAX = 32 * ROUND_CELL  # 2^-14 s


@dataclass(frozen=True)
class BoundCertificate:
    supply: TAU.TauRoundoffSupply
    target_error: F
    bound: F = UNIFORM_RESIDUAL_MAX
    def __post_init__(self):
        te=abs(F(self.target_error)); b=F(self.bound)
        object.__setattr__(self,'target_error',te); object.__setattr__(self,'bound',b)
        if b != UNIFORM_RESIDUAL_MAX:
            raise ValueError('tau roundoff theorem bound is source-fixed')
        if te > TARGET_ERROR_MAX:
            raise ValueError('tau target binary32/exact discrepancy exceeds certified cell')
        if abs(self.supply.residual_fma) > FMA_DERIVED_MAX:
            raise ValueError('FMA tau residual exceeds analytic derived bound')
        if abs(self.supply.residual_separate) > SEPARATE_DERIVED_MAX:
            raise ValueError('separate tau residual exceeds analytic derived bound')
        if max(abs(self.supply.residual_fma),abs(self.supply.residual_separate)) > b:
            raise ValueError('tau roundoff supply exceeds uniform theorem bound')


def _domain(binary:TAU.TauStep, exact_target):
    if not isinstance(binary,TAU.TauStep): raise TypeError('TauStep required')
    t=F(exact_target)
    if abs(binary.previous) > PREVIOUS_ABS_MAX:
        raise ValueError('previous tau outside certified scalar domain')
    if abs(t) > TARGET_ABS_MAX:
        raise ValueError('exact tau target outside certified scalar domain')
    if not 0 <= F(1)-binary.exp_decay <= 1 or not 0 <= binary.alpha <= 1:
        raise ValueError('tau EMA alpha outside convex shipping domain')
    return t


def certify(binary:TAU.TauStep, exact, *, exact_target):
    """Legacy bridge against the exact tuner CandidateResult."""
    t=_domain(binary,exact_target)
    supply=TAU.roundoff_supply(binary,exact)
    if F(exact.tau_target) != t:
        raise ValueError('exact target argument detached from exact candidate')
    return BoundCertificate(supply,binary.tau_target-t)


def certify_source_target(binary:TAU.TauStep, target:TARGET.TargetPair):
    """Certify deployment supply from the same source-locked target expression."""
    if not isinstance(target,TARGET.TargetPair): raise TypeError('source-locked TargetPair required')
    t=_domain(binary,target.exact_target)
    if binary.frequency != target.clamped_frequency:
        raise ValueError('binary tau frequency detached from source-locked target pair')
    if binary.tau_target != target.binary32_target:
        raise ValueError('binary tau target detached from source-locked target pair')
    a=F(1)-binary.exp_decay
    real_next=binary.previous+a*(t-binary.previous)
    supply=TAU.TauRoundoffSupply(real_next,binary.next_separate,binary.next_fma,
                                  binary.next_separate-real_next,
                                  binary.next_fma-real_next)
    return BoundCertificate(supply,target.error)


def readiness():
    target=TARGET.readiness()
    return {
      'normal_binary32_rounding_cell_bound_for_abs_intermediate_lt_32':True,
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
