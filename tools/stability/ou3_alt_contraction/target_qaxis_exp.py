"""All-input error bound for the pinned newlib expf Qaxis branch.

The constant words and operation graph are those in libm_a-ef_exp.o of
esp-14.2.0_20251107's esp32s3 multilib. For z=-x, RN32(.01)<=x<=1/4,
this is the k=0 branch of newlib/libm/math/ef_exp.c at commit
9a0d39153510ec5cbb51eb8c70cecbfeffdbb6ba. No range reduction occurs.

The algebraic certificate uses exact rational polynomial coefficients, Taylor
remainder and forward operation bounds; it does not sample transcendental
values. ADD/SUB/MUL require binary32 round-to-nearest; multiply/add may have
one or two roundings. Division may err by one ulp. Every intermediate is
normal. These scalar instruction/library-helper and final-link premises are
shared deployment obligations; this module does not assert them from hashes.
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

QUALIFICATION = 'OU3_ALT_PINNED_NEWLIB_QAXIS_EXP_ALL_INPUT_V1'
NEWLIB_COMMIT = '9a0d39153510ec5cbb51eb8c70cecbfeffdbb6ba'
SOURCE_URL = ('https://github.com/espressif/newlib-esp32/blob/' + NEWLIB_COMMIT
              + '/newlib/libm/math/ef_exp.c')
SOURCE_SHA256 = 'a7bd3ce3494b7e0099e6f159736fbd1215224eea6600bab650db4654b1313efb'
LIBM_SHA256 = '98e35630b7908f49e8d845e2ba3b1f44920508d4774a41b80eee617a1592a6bb'
MEMBER_SHA256 = 'e983eeb4192c32d1730144257a765dae3ed63c2b428433acc931a77f77222c8c'
COEFFICIENT_WORDS = (0x3e2aaaab, 0xbb360b61, 0x388ab355,
                     0xb5ddea0e, 0x3331bb4c)
P = tuple(F(struct.unpack('>f', struct.pack('>I', w))[0]) for w in COEFFICIENT_WORDS)
X_MIN = B.rn32(F(1, 100))
X_MAX = F(1, 4)
U = F(1, 2**24)


def _polymul(a, b):
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def error_certificate():
    """Prove |target_graph(-x)-exp(-x)|<2^-24 under scalar premises.

    Real graph: t=z², p=P1+t(P2+t(P3+t(P4+tP5))), c=z-tp,
    r=(zc)/(c-2)-z, y=1-r. Thus y=1+2z/(2-c).
    Compare that rational function with exp's degree-12 Taylor polynomial
    after multiplying by 2-c. Bound every residual coefficient exactly.
    """
    a = X_MAX
    t = a*a
    # P(z²)>0 throughout; hence c<0 and 2-c>=2 for z<0.
    p_lower = P[0] - sum(abs(v)*t**i for i, v in enumerate(P[1:], 1))
    assert p_lower > F(1, 7)
    denominator = [F(2), F(-1)]
    for p in P:
        denominator += [p, F(0)]
    denominator.pop()
    numerator = denominator.copy()
    numerator[1] += 2
    taylor = [F(1, factorial(i)) for i in range(13)]
    residual = _polymul(denominator, taylor)
    for i, q in enumerate(numerator):
        residual[i] -= q
    polynomial_error = sum(abs(q)*a**i for i, q in enumerate(residual))/2
    # exp(a)<2, so Lagrange's remainder is at most 2*a^13/13!.
    approximation_error = polynomial_error + 2*a**13/factorial(13)

    # For each Horner stage, m bounds its exact real magnitude, e the
    # computed-minus-real value. Two-rounding multiplication/addition bounds
    # also contain a single rounding MADD. Inputs and outputs are normal.
    t_error = U*t
    m, e = abs(P[-1]), F(0)
    horner_errors = []
    for p in reversed(P[:-1]):
        product_magnitude = (t+t_error)*(m+e)
        local = U*product_magnitude + U*(abs(p)+product_magnitude*(1+U))
        e = (t+t_error)*e + t_error*m + local
        m = abs(p)+t*m
        horner_errors.append(e)
    product_error = ((t+t_error)*e + t_error*m
                     + U*(t+t_error)*(m+e))
    c_magnitude = a+t*m
    c_error = product_error+U*(c_magnitude+product_error)
    n_magnitude = a*c_magnitude
    n_error = a*c_error+U*a*(c_magnitude+c_error)
    # c-2 is between -4 and -2: half its binade spacing is 2^-23.
    d_error = c_error+F(1, 2**23)
    assert c_magnitude+c_error < 1
    assert d_error < 1
    # |quotient|<1/16 (also after operand errors), so one division ulp
    # is at most 2^-28. This permits non-correctly-rounded division.
    quotient_magnitude = (n_magnitude+n_error)/(2-d_error)
    assert quotient_magnitude < F(1, 16)
    quotient_error = (n_error/(2-d_error)
                      + n_magnitude*d_error/(2*(2-d_error))
                      + F(1, 2**28))
    # c, numerator and denominator preserve their signs. The quotient is
    # negative and |quotient|<|z|, with a margin already at minimum |z|.
    assert c_error < X_MIN/2
    assert quotient_error < X_MIN*(1-c_magnitude/2)/2
    # Thus 0<quotient-z<1/4 and 1-(quotient-z) lies strictly in (1/2,1).
    # Last SUB half-ulps are respectively <=2^-27 and <=2^-25.
    rounding_error = quotient_error+F(1, 2**27)+F(1, 2**25)
    error = approximation_error+rounding_error
    assert error < U
    return {
        'qualification': QUALIFICATION,
        'argument_x_interval': (X_MIN, X_MAX),
        'newlib_source_commit': NEWLIB_COMMIT,
        'coefficient_words': COEFFICIENT_WORDS,
        'rational_approximation_error': approximation_error,
        'horner_stage_absolute_errors': tuple(horner_errors),
        'c_absolute_error': c_error,
        'division_absolute_error_allowance': F(1, 2**28),
        'total_absolute_exp_error': error,
        'required_absolute_exp_error': U,
        'positive_error_margin': U-error,
        'all_input_polynomial_and_roundoff_bound_proved': True,
        'range_reduction_or_subnormal_intermediate_used': False,
        'both_fused_and_separate_horner_rounding_covered': True,
        'division_correct_rounding_required': False,
        'scalar_RNE_and_one_ulp_division_premises_discharged_here': False,
        'final_firmware_link_correspondence_discharged_here': False,
    }


def audit(libm: Path, ar: Path, objdump: Path):
    """Bind the analyzed branch to exact verified SDK member bytes.

    Hash identity is supplemented by the coefficient literals, no-reduction
    branch and helper relocations. It establishes the supplying program,
    not processor instruction semantics or final firmware selection.
    """
    if hashlib.sha256(libm.read_bytes()).hexdigest() != LIBM_SHA256:
        raise ValueError('unqualified target libm archive bytes')
    member = subprocess.check_output([str(ar), 'p', str(libm), 'libm_a-ef_exp.o'])
    if hashlib.sha256(member).hexdigest() != MEMBER_SHA256:
        raise ValueError('unqualified expf object member bytes')
    dis = subprocess.check_output([str(objdump), '-dr', str(libm)], text=True)
    section = dis.split('libm_a-ef_exp.o:', 1)[1].split('\nlibm_a-', 1)[0]
    for word in COEFFICIENT_WORDS:
        if f'{word:08x}' not in section:
            raise ValueError('expf polynomial literal absent')
    for op in ('mul.s', 'madd.s', 'msub.s', '__divsf3', '3eb17218'):
        if op not in section:
            raise ValueError('expf target branch correspondence missing: '+op)
    result = error_certificate()
    result.update({
        'libm_sha256': LIBM_SHA256,
        'expf_member_sha256': MEMBER_SHA256,
        'newlib_source_url': SOURCE_URL,
        'newlib_source_sha256': SOURCE_SHA256,
        'exact_target_object_branch_identified': True,
        'expf_no_reduction_branch_entry': '0xe4 (k=0), then 0xe6',
        'expf_horner_instruction_interval': '0xe6..0x124',
        'expf_quotient_and_return_instruction_interval': '0x12f..0x151',
        'target_instruction_semantics_qualified': False,
        'whole_firmware_link_resolution_qualified': False,
    })
    return result


def _json(value):
    if isinstance(value, F):
        return {'numerator': value.numerator, 'denominator': value.denominator,
                'decimal': float(value)}
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--libm', type=Path, required=True)
    parser.add_argument('--ar', type=Path, required=True)
    parser.add_argument('--objdump', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.libm, args.ar, args.objdump)
    args.output.write_text(json.dumps(result, indent=2, default=_json)+'\n')
    print('All-input Qaxis exp graph absolute error:',
          float(result['total_absolute_exp_error']))
    print('Instruction semantics and final link remain separate obligations.')


if __name__ == '__main__':
    main()
