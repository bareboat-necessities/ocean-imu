"""Two shipping Qaxis covariance ``std::exp(-x)`` calls, source-bound.

In the IntegratedOUChain<3> general branch shipping first calls the nested
IntegratedOUChain<2>::process_covariance, which evaluates ``std::exp(-x)``, and
then evaluates ``std::exp(-x)`` again for the additional S row/column.  These
are not the OU transition's make_prims() call and must not be replaced by its
``alpha`` merely because the mathematical argument is the same.

The two results are therefore retained independently.  Each must be an actual
binary32 value and lie in a rigorous real enclosure for the SAME binary32 x
already produced by the Qaxis branch graph.  No bit identity between the two
calls, or with the OU transition exp, is assumed.  The pinned newlib target graph now has an all-input error bound;
its instruction semantics and final-link attachment remain separate deployment
obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as QB
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_common/KalmanOUCoreMath.h'
QUALIFICATION='OU3_ALT_QAXIS_TWO_EXP_BINARY32_V1'


@dataclass(frozen=True)
class ExpPair:
    branch: QB.Branch
    marginal:F
    final:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.branch,QB.Branch): raise TypeError('binary32 Qaxis branch required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong Qaxis exp qualification')
        m=F(self.marginal); f=F(self.final)
        if not B.is_binary32(m) or not B.is_binary32(f):
            raise ValueError('Qaxis exp witnesses must be binary32 values')
        if self.branch.small:
            raise ValueError('small-series Qaxis branch executes no covariance std::exp calls')
        lo,hi,_,_=EXP.enclosure(self.branch.x)
        if not lo<=m<=hi: raise ValueError('nested Qaxis exp detached from SAME binary32 x')
        if not lo<=f<=hi: raise ValueError('final Qaxis exp detached from SAME binary32 x')
        object.__setattr__(self,'marginal',m); object.__setattr__(self,'final',f)


def pair(branch:QB.Branch,*,marginal,final): return ExpPair(branch,F(marginal),F(final))


def general_branch_exp_error_budget():
    """Uniform sufficient libm budget, NOT evidence that target libm meets it.

    On t=RN32(.01)<=x<=1, alternating Taylor sums give
    exp(-x)-(1-x) >= x²/2-x³/6 >= t²/3,
    (1-x+x²/2)-exp(-x) >= x³/6-x⁴/24 >= t³/8.
    Thus a separately proved absolute target exp error <=2^-24 leaves BOTH
    original real-envelope inequalities valid, including the threshold word.
    This also covers canonical Qaxis x<=.25. No equality of repeated calls
    or correctly-rounded target exp assumption is introduced.
    """
    t=QB.SMALL_THRESHOLD
    budget=F(1,2**24)
    lower=t*t/3
    upper=t*t*t/8
    if min(lower,upper)<=budget:
        raise AssertionError('Qaxis source envelope has insufficient libm slack')
    return {'argument_interval':(t,F(1)),
            'sufficient_absolute_exp_error':budget,
            'lower_boundary_slack_after_error':lower-budget,
            'upper_boundary_slack_after_error':upper-budget,
            'source_uniform_sufficient_error_budget_proved':True,
            'target_libm_satisfies_budget_proved':False}


def pinned_target_graph_envelope():
    """Compose the actual target exp graph bound with Qaxis real slack.

    This is an all-input library lemma with shared scalar/link premises,
    not an assertion that supplied witnesses are the final firmware output.
    """
    from tools.stability.ou3_alt_contraction import target_qaxis_exp as TARGET
    proof=TARGET.error_certificate()
    from tools.stability.ou3_alt_contraction import target_wpe_libm as LIBM
    profile=LIBM.profile_correspondence()
    budget=general_branch_exp_error_budget()
    domain=proof['argument_x_interval']
    if domain[0]!=QB.SMALL_THRESHOLD or domain[1]!=F(1,4):
        raise AssertionError('target exp certificate detached from canonical Qaxis domain')
    if proof['total_absolute_exp_error']>=budget['sufficient_absolute_exp_error']:
        raise AssertionError('target exp graph exceeds Qaxis source envelope slack')
    return {
        'all_input_pinned_target_graph_exp_error_proved':True,
        'Qaxis_pinned_libm_approximation_correspondence_closed':(
            profile['pinned_scalar_profile_attached'] and profile['pinned_math_namespace_attached']),
        'target_graph_results_satisfy_original_Qaxis_envelope_under_scalar_link_premises':True,
        'argument_interval':domain,
        'absolute_exp_error':proof['total_absolute_exp_error'],
        'shared_scalar_and_final_link_premises_discharged_here':False,
        'named_RNE_scalar_profile_attached':profile['pinned_scalar_profile_attached'],
        'whole_firmware_compiler_and_link_correspondence_closed':False,
    }


def _source_shape_matches():
    s=SOURCE.read_text()
    # One call is in IntegratedOUChain<T,2>; the other is in T<3>'s general
    # branch. make_prims contains another exp call, so require at least three
    # total rather than pretending the covariance pair is the only use.
    return s.count('std::exp(-x)')>=3 and 'IntegratedOUChain<T,2>::process_covariance(tau_eff, h, sigma2, marginal);' in s


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_nested_and_final_Qaxis_covariance_exp_calls_present':_source_shape_matches(),
      'Qaxis_nested_exp_retained_separately_from_OU_transition_exp':True,
      'Qaxis_final_exp_retained_separately_from_OU_transition_exp':True,
      'Qaxis_nested_and_final_exp_bit_identity_assumed':False,
      'both_Qaxis_exp_results_bound_to_same_binary32_x_real_enclosure':True,
      'Qaxis_general_branch_sufficient_libm_error_budget':general_branch_exp_error_budget(),
      'Qaxis_pinned_target_graph_real_envelope':pinned_target_graph_envelope(),
      'Qaxis_pinned_libm_approximation_correspondence_closed':
          pinned_target_graph_envelope()['Qaxis_pinned_libm_approximation_correspondence_closed'],
      'Qaxis_exp_libm_binary32_correspondence_closed':False,
      'ALT_LIVE_PASS':False,
    }
