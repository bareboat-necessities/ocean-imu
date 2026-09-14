"""Two shipping Qaxis covariance ``std::exp(-x)`` calls, source-bound.

In the IntegratedOUChain<3> general branch shipping first calls the nested
IntegratedOUChain<2>::process_covariance, which evaluates ``std::exp(-x)``, and
then evaluates ``std::exp(-x)`` again for the additional S row/column.  These
are not the OU transition's make_prims() call and must not be replaced by its
``alpha`` merely because the mathematical argument is the same.

The two results are therefore retained independently.  Each must be an actual
binary32 value and lie in a rigorous real enclosure for the SAME binary32 x
already produced by the Qaxis branch graph.  No bit identity between the two
calls, or with the OU transition exp, is assumed.  Target-libm error remains a
separate deployment obligation.
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
      'Qaxis_exp_libm_binary32_correspondence_closed':False,
      'ALT_LIVE_PASS':False,
    }
