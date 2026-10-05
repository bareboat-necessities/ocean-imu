"""Exact frontend mismatch witnesses; no shipping interval binding is asserted.

These regressions refute use of the old idealized reference box as a literal
float32 frontend. They do not refute the filter, the local theorem, or the
pre-normalization quaternion identity. General IEEE-754 transfer remains open.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import hashlib

ROOT=Path(__file__).resolve().parents[3]

class ShippingBindingFailure(RuntimeError):
    """A reference implementation has not been bound to the shipping program (E)."""


def dyadic(bits: int) -> F:
    sign=-1 if bits>>31 else 1
    exponent=(bits>>23)&255;mantissa=bits&0x7fffff
    if exponent==255:raise ValueError('nonfinite binary32')
    if exponent:mantissa+=1<<23;power=exponent-150
    else:power=-149
    return sign*F(mantissa)*(F(2)**power)


def round32(value: F) -> int:
    """Exact round-to-nearest, ties-to-even binary32, including subnormals."""
    value=F(value)
    if value==0:return 0
    sign=(1<<31) if value<0 else 0;value=abs(value)
    e=value.numerator.bit_length()-value.denominator.bit_length()
    if value<F(2)**e:e-=1
    shift=max(e-23,-149);scaled=value/(F(2)**shift)
    m,rem=divmod(scaled.numerator,scaled.denominator)
    if 2*rem>scaled.denominator or (2*rem==scaled.denominator and m%2):m+=1
    if m==1<<24:m>>=1;shift+=1
    if m<(1<<23):return sign|m
    exponent=shift+150
    if exponent>=255:raise OverflowError('binary32 overflow')
    return sign|(exponent<<23)|(m-(1<<23))


def rn(value: F) -> F:
    return dyadic(round32(value))


def inverse_sqrt_word(bits: int) -> tuple[int,list[str]]:
    """Literal non-fused scalar expression in Mahony_AHRS<float>::invSqrt."""
    n=dyadic(bits)
    if n<=0:raise ValueError('positive normal input required')
    y=dyadic(0x5f375a86-(bits>>1))
    x=rn(n*F(1,2));a=rn(x*y);b=rn(a*y);c=rn(F(3,2)-b);out=rn(y*c)
    return round32(out),[str(v) for v in (n,y,x,a,b,c,out)]


def certificate() -> dict:
    bits,steps=inverse_sqrt_word(0x3f800000);q=dyadic(bits)
    # Initial q=(1,0,0,0), zero gyro and aligned nonzero acceleration have zero
    # Mahony feedback; the first quaternion normalization is invSqrt(1).
    norm_sq=q*q
    assert norm_sq<1
    return {
        'qualification':'OU3_PLANAR_FRONTEND_BINDING_AUDIT_V1',
        'result_type':'PROVED exact arithmetic counterexample to an idealized proof replica',
        'classification':'E_IMPLEMENTATION_DIAGNOSTIC_FAILURE',
        'shipping_counterexample':False,
        'binary32_scope':'round-to-nearest ties-to-even, literal non-fused scalar evaluation',
        'input_bits':'0x3f800000','inverse_sqrt_bits':f'0x{bits:08x}',
        'exact_operation_values':steps,'raw_quaternion_squared_norm':str(norm_sq),
        'raw_quaternion_squared_norm_float':float(norm_sq),
        'unit_norm_after_shipping_float_normalization':False,
        'pre_normalization_quaternion_identity_retained':True,
        'unbound_reference_operations':[
            'MahonyBox substitutes exact reciprocal square roots for the shipping float bit seed/Newton map',
            'MahonyBox returns -down.dot(acc) without the shipping gravity subtraction',
            'initial_adaptation does not execute the literal accelerometer tilt seeding/readiness gates',
            'ClosedCausalAdaptationBox advances period before tuner use instead of the literal staged chronology',
            'tuner profile constants, S cadence and independent AW-sync clock are not source-bound'],
        'required_cell_coordinates':[
            'planar nominal mean','P_even','P_odd','raw private Mahony quaternion and integral',
            'guard detector/conditioning state','frequency/variance/period state',
            'staged and applied joint tuner tuple','S scheduler elapsed',
            'AW synchronization clock and pending target','reference/refinement/gate state'],
        'shipping_frontend_interval_binding_verified':False,
        'all_time_magnetic_service_verified':False,'theorem_closed':False,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/ahrs/Mahony_AHRS.h','src/tuner/VerticalAccelComplementary.h',
            'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h')},
    }


def require_shipping_binding() -> None:
    raise ShippingBindingFailure('the interval frontend is an unbound ideal reference: '+
        '; '.join(certificate()['unbound_reference_operations']))


if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
