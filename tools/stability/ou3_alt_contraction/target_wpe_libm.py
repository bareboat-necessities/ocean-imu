"""All-input newlib expf/logf approximation bounds for the WPE profile.

Analyzes the pinned ef_exp.o/ef_log.o algorithm, not the host libm. The full
binary32 operation error is included in exp relative 2^-20 and log absolute
2^-14; correct rounding of transcendentals is unnecessary. Fixed-coefficient
polynomials are compared algebraically with convergent series. All domain
covering is analytic, with no sampled input histories or subdivision.

Common deployment premises: RNE scalar ADD/SUB/MUL; one-ulp division;
MADD/MSUB have either fused or separate RNE semantics; normal finite inputs
and outputs of these operations; the final firmware selects these objects.
Target ISA and link qualification discharge those premises elsewhere.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path
import struct
import subprocess

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import target_qaxis_exp as Q

QUALIFICATION = 'OU3_ALT_PINNED_NEWLIB_WPE_EXP_LOG_ALL_INPUT_V1'
EXP_DOMAIN = (F(-60), F(60))
LOG_DOMAIN = (F(1, 2**62), F(2**16))
LOG_MEMBER_SHA256 = 'd2cda42cc88f4125cce59dd1b45adc0cb26a5e6c86bc90933009afef4baf83be'
LOG_SOURCE_SHA256 = '5857e5cc44d7e8808a4378540418f3e4420f290582dc0ffa6a6eb4d4489df340'


def _word(w):
    return F(struct.unpack('>f', struct.pack('>I', w))[0])


LN2_HI, LN2_LO, INV_LN2 = map(_word, (0x3f317180, 0x3717f7d1, 0x3fb8aa3b))
LG_WORDS = (0x3f2aaaab, 0x3ecccccd, 0x3e924925, 0x3e638e29,
            0x3e3a3325, 0x3e1cd04f, 0x3e178897)
LG = tuple(map(_word, LG_WORDS))
U = Q.U


class _Error:
    """Magnitude of ideal expression plus forward absolute rounding error.

    Subtraction uses the addition bound. This deliberately sacrifices
    cancellation only where the final proven budget has sufficient margin.
    Two rounded multiply/add operations dominate either contracted choice.
    """
    def __init__(self, magnitude, error=0):
        self.m, self.e = F(magnitude), F(error)
        if min(self.m, self.e) < 0:
            raise ValueError('nonnegative error bounds required')

    def __add__(self, other):
        if not isinstance(other, _Error):
            other = _Error(abs(other))
        m, e = self.m+other.m, self.e+other.e
        return _Error(m, e+U*(m+e))

    __sub__ = __add__

    def __mul__(self, other):
        if not isinstance(other, _Error):
            other = _Error(abs(other))
        m = self.m*other.m
        e = self.m*other.e+other.m*self.e+self.e*other.e
        return _Error(m, e+U*(m+e))


def ln2_enclosure():
    # log(2)=2*atanh(1/3). Positive geometric tail after j=29.
    lower = 2*sum(F(1, (2*j+1)*3**(2*j+1)) for j in range(30))
    upper = lower+F(2, 61*3**61)*F(9, 8)
    return lower, upper


def exp_certificate():
    ln_lo, ln_hi = ln2_enclosure()
    split_error = max(abs(LN2_HI+LN2_LO-ln_lo), abs(LN2_HI+LN2_LO-ln_hi))
    inverse_error = max(abs(INV_LN2-1/ln_lo), abs(INV_LN2-1/ln_hi))
    argument_max, k_max = F(60), 88
    # k=trunc(RN(a*invln2+sign(a)/2)); both compiler contraction choices.
    # Account for the constant error and both rounding sites before trunc.
    k_round_error = argument_max*inverse_error+2*U*(argument_max*INV_LN2+1)
    assert argument_max/ln_lo+F(1, 2)+k_round_error < k_max
    # The separate +/-one reduction branch uses binary32 half/1.5ln2
    # thresholds. Their deviations also fit the same k selection budget.
    half_threshold = _word(0x3eb17218)
    one_half_threshold = _word(0x3f851592)
    assert half_threshold/ln_hi > F(1, 2)-k_round_error
    assert one_half_threshold/ln_lo < F(3, 2)+k_round_error
    # Unlike a generic product at magnitude60, k*ln2HI is exact: every
    # possible integer k in this finite branch has a binary32 product.
    assert all(B.is_binary32(k*LN2_HI) for k in range(-k_max, k_max+1))
    reduced_real_max = (F(1, 2)+k_round_error)*ln_hi
    hi_max = reduced_real_max+k_max*(ln_hi-LN2_HI)
    lo_max = k_max*LN2_LO
    residual_round_error = U*(hi_max*(1+U)+lo_max*(1+U))
    argument_error = (k_max*split_error+U*hi_max+U*lo_max
                      + residual_round_error)
    a = F(7, 20)
    assert reduced_real_max+argument_error < a
    assert half_threshold < a

    # The reduced real exp approximant is 1+2z/(2-c). Denominator>=1.6.
    denominator = [F(2), F(-1)]
    for p in Q.P:
        denominator += [p, F(0)]
    denominator.pop()
    numerator = denominator.copy()
    numerator[1] += 2
    residual = Q._polymul(denominator, [F(1, factorial(i)) for i in range(15)])
    for i, q in enumerate(numerator):
        residual[i] -= q
    approximation_error = (sum(abs(v)*a**i for i, v in enumerate(residual))/F(8, 5)
                           + 2*a**15/factorial(15))
    x = _Error(a)
    t = x*x
    p = _Error(abs(Q.P[-1]))
    for coefficient in reversed(Q.P[:-1]):
        p = _Error(abs(coefficient))+t*p
    c = x+t*p
    assert c.m+c.e < F(2, 5)
    n, d = x*c, _Error(2)+c
    d_min = F(8, 5)
    # Quotient magnitude<1/8: one division ulp <=2^-27.
    assert (n.m+n.e)/(d_min-d.e) < F(1, 8)
    quotient_error = (n.e/(d_min-d.e)
                      + n.m*d.e/(d_min*(d_min-d.e))+F(1, 2**27))
    q = _Error(n.m/d_min, quotient_error)
    k_zero = _Error(1)+(q+x)
    k_nonzero = _Error(1)+((_Error(lo_max)+q)+_Error(hi_max))
    local_error = max(k_zero.e, k_nonzero.e)
    # The nonzero-k graph reconstructs hi-lo instead of z=RN(hi-lo),
    # so retain that compensation error explicitly. exp(z)>=1-|z|>=.65.
    # For |delta|<1/2, |exp(delta)-1|<=2|delta|. Exponent-bit scaling
    # is exact here: |k|<=88 and .65<y<1.5 imply normal finite output.
    assert argument_error < F(1, 2)
    relative_error = ((local_error+residual_round_error+approximation_error)
                      /F(13, 20)+2*argument_error)
    # |a|<2^-23 selects return 1+a. This covers zero and subnormal a,
    # including input flush to zero: both possible returns differ from
    # exp(a) by at most 2^-22 relative (well within the declared budget).
    tiny_branch_relative_error = F(1, 2**22)
    relative_error = max(relative_error, tiny_branch_relative_error)
    budget = F(1, 2**20)
    assert relative_error < budget
    return {
        'argument_interval': EXP_DOMAIN,
        'integer_range_reduction_bound': k_max,
        'reduced_argument_magnitude_bound': a,
        'range_reduction_argument_error': argument_error,
        'rational_approximation_error': approximation_error,
        'operation_absolute_error': local_error,
        'total_relative_exp_error': relative_error,
        'required_relative_exp_error': budget,
        'positive_relative_error_margin': budget-relative_error,
        'all_input_exp_error_bound_proved_under_scalar_premises': True,
        'transcendental_correct_rounding_assumed': False,
    }


def log_certificate():
    ln_lo, ln_hi = ln2_enclosure()
    split_error = max(abs(LN2_HI+LN2_LO-ln_lo), abs(LN2_HI+LN2_LO-ln_hi))
    k_max = 62
    # Exact bit normalization changes exponent and possibly halves the
    # mantissa; all admitted log inputs are normal. No transcendental
    # range-reduction approximation is involved.
    threshold = 0x800000-(0x95f64 << 3)
    normalized_min = (1+F(threshold, 2**23))/2
    normalized_max = 1+F(threshold-1, 2**23)
    assert F(7, 10) < normalized_min < normalized_max < F(3, 2)
    f_max, s_max = F(1, 2), F(1, 5)
    # f=x-1 is exact by Sterbenz. The denominator 2+f is in[1.7,2.5];
    # one division ulp below |s|=.2 is <=2^-25.
    denominator_error = U*F(5, 2)
    s_error = (f_max*denominator_error
               /(F(17, 10)*(F(17, 10)-denominator_error))+F(1, 2**25))
    # Both algebraic log branches equal 2s+sR, R=sum Lg_i*s^(2i).
    # Compare with log(x)=2*atanh(s); no table/sample approximation.
    approximation_error = (
        sum(abs(p-F(2, 2*i+3))*s_max**(2*i+3) for i, p in enumerate(LG))
        + 2*s_max**17/(17*(1-s_max*s_max)))
    f, s = _Error(f_max), _Error(s_max, s_error)
    z, w = s*s, (s*s)*(s*s)
    t1 = w*(_Error(LG[1])+w*(_Error(LG[3])+w*_Error(LG[5])))
    t2 = z*(_Error(LG[0])+w*(_Error(LG[2])+w*(_Error(LG[4])+w*_Error(LG[6]))))
    r = t1+t2
    hfsq = _Error(F(1, 2))*f*f
    kh, kl = _Error(k_max)*_Error(LN2_HI), _Error(k_max)*_Error(LN2_LO)
    branch_errors = (
        (f+(hfsq+s*(hfsq+r))).e,
        (kh+((hfsq+(s*(hfsq+r)+kl))+f)).e,
        (f+s*(f+r)).e,
        (kh+((s*(f+r)+kl)+f)).e,
    )
    # The small-f branch tests only mantissas ix=0 and the last15 words.
    # Enumerating these integer branch selectors proves |f|<=2^-20;
    # this is not sampling the continuum or motion history.
    small_mantissas = [0]+list(range(2**23-15, 2**23))
    for ix in small_mantissas:
        assert ((15+ix)&0x7fffff) < 16
        normalized = (1+F(ix, 2**23))/(2 if ix >= threshold else 1)
        assert abs(normalized-1) <= F(1, 2**20)
    tiny = F(1, 2**20)
    tf = _Error(tiny)
    third = B.rn32(F(1, 3))
    tr = tf*tf*(_Error(F(1, 2))+_Error(third)*tf)
    tiny_approximation_error = (tiny**4/(4*(1-tiny))
                                + abs(third-F(1, 3))*tiny**3)
    tiny_round_error = max((tf+tr).e, (kh+((tr+kl)+tf)).e, (kh+kl).e)
    total = max(max(branch_errors)+approximation_error,
                tiny_round_error+tiny_approximation_error)+k_max*split_error
    budget = F(1, 2**14)
    assert total < budget
    return {
        'argument_interval': LOG_DOMAIN,
        'normalized_argument_interval': (normalized_min, normalized_max),
        'normalization_exponent_absolute_bound': k_max,
        'atanh_polynomial_approximation_error': approximation_error,
        'large_branch_operation_absolute_errors': branch_errors,
        'tiny_branch_operation_absolute_error': tiny_round_error,
        'total_absolute_log_error': total,
        'required_absolute_log_error': budget,
        'positive_absolute_error_margin': budget-total,
        'all_input_log_error_bound_proved_under_scalar_premises': True,
        'transcendental_correct_rounding_assumed': False,
    }


def certificate():
    return {
        'qualification': QUALIFICATION,
        'newlib_source_commit': Q.NEWLIB_COMMIT,
        'exp': exp_certificate(),
        'log': log_certificate(),
        'scalar_instruction_and_division_premises_discharged_here': False,
        'whole_firmware_link_resolution_qualified_here': False,
        'WPE_sqrt_error_qualification_supplied_here': False,
    }


def audit(libm: Path, ar: Path, objdump: Path, source: Path):
    if hashlib.sha256(source.read_bytes()).hexdigest() != LOG_SOURCE_SHA256:
        raise ValueError('unqualified pinned logf source bytes')
    exp_report = Q.audit(libm, ar, objdump)
    member = subprocess.check_output([str(ar), 'p', str(libm), 'libm_a-ef_log.o'])
    if hashlib.sha256(member).hexdigest() != LOG_MEMBER_SHA256:
        raise ValueError('unqualified logf object member bytes')
    dis = subprocess.check_output([str(objdump), '-dr', str(libm)], text=True)
    section = dis.split('libm_a-ef_log.o:', 1)[1].split('\nlibm_a-', 1)[0]
    for word in LG_WORDS:
        if f'{word:08x}' not in section:
            raise ValueError('logf polynomial literal absent')
    for op in ('mul.s', 'madd.s', 'msub.s', '__divsf3'):
        if op not in section:
            raise ValueError('logf target operation absent: '+op)
    result = certificate()
    result.update({
        'libm_sha256': Q.LIBM_SHA256,
        'expf_member_sha256': Q.MEMBER_SHA256,
        'logf_member_sha256': LOG_MEMBER_SHA256,
        'logf_source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'logf_source_url': ('https://github.com/espressif/newlib-esp32/blob/'
                           +Q.NEWLIB_COMMIT+'/newlib/libm/math/ef_log.c'),
        'exact_target_exp_and_log_objects_identified': exp_report['exact_target_object_branch_identified'],
        'source_object_correspondence_method': 'pinned bytes, constants and inspected branch instruction graph',
    })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--libm', type=Path, required=True)
    parser.add_argument('--ar', type=Path, required=True)
    parser.add_argument('--objdump', type=Path, required=True)
    parser.add_argument('--log-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.libm, args.ar, args.objdump, args.log_source)
    args.output.write_text(json.dumps(result, default=Q._json, indent=2)+'\n')
    print('All-input exp relative error:', float(result['exp']['total_relative_exp_error']))
    print('All-input log absolute error:', float(result['log']['total_absolute_log_error']))


if __name__ == '__main__':
    main()
