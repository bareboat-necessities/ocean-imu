"""Exact finite binary32 graph for the initialized private Mahony observer.

The graph is conditional on PROFILE, not on freely chosen reciprocal witnesses.
Each elementary operation rounds an exact rational to binary32, and each
normalization executes the shipping integer seed and single Newton step.
Rounding defects are retained on their actual operands; they are NOT independently
bounded ISS inputs. Source admission, startup seeding, target compiler profile,
nonfinite intermediates, and whole-word stability remain open.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction as F

PROFILE = 'binary32-rne-gradual-no-fma-eigen3-scalar'
MAGIC = 0x5F375A86
SIGN = 1 << 31
MAX_FINITE = 0x7F7FFFFF


def rational(x):
    if isinstance(x, (float, bool)):
        raise TypeError('use an exact rational or explicit binary32 bit pattern')
    return F(x)


def pow2(e: int) -> F:
    return F(1 << e) if e >= 0 else F(1, 1 << -e)


def value(bits: int) -> F:
    """Decode a finite IEEE binary32 word exactly, including subnormal values."""
    if type(bits) is not int or not 0 <= bits < 1 << 32:
        raise ValueError('one uint32 encoding required')
    e, m = (bits >> 23) & 255, bits & ((1 << 23) - 1)
    if e == 255:
        raise ValueError('nonfinite binary32 is outside this finite graph')
    z = F(m) * pow2(-149) if e == 0 else F((1 << 23) + m) * pow2(e - 150)
    return -z if bits & SIGN else z


def _rne_integer(z: F) -> int:
    q, r = divmod(z.numerator, z.denominator)
    return q + int(2*r > z.denominator or (2*r == z.denominator and q & 1))


def round_bits(x, *, negative_zero=False) -> int:
    """Nearest, ties-to-even; finite result or an explicit overflow exception.

    No host float operation is used. This is rounding to the complete binary32
    lattice, not an interval or decimal approximation to that lattice.
    """
    x = rational(x)
    if not x:
        return SIGN if negative_zero else 0
    sign = SIGN if x < 0 else 0
    a = abs(x)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if a < pow2(e):
        e -= 1
    shift = max(-149, e - 23)
    sig = _rne_integer(a / pow2(shift))
    if sig == 0:
        return sign
    if shift == -149 and sig < 1 << 23:
        return sign | sig
    if sig == 1 << 24:
        sig >>= 1
        shift += 1
    exponent = shift + 150
    if exponent >= 255:
        raise OverflowError('finite binary32 graph encountered overflow')
    if not 1 <= exponent <= 254 or not 1 << 23 <= sig < 1 << 24:
        raise AssertionError('binary32 rounding significand invariant lost')
    return sign | (exponent << 23) | (sig - (1 << 23))


def rounding_cell_contains(x, bits: int) -> bool:
    """Independent nearest-neighbor certificate, without calling round_bits.

    The virtual neighbor of MAX_FINITE is 2**128; its midpoint is the IEEE
    round-to-nearest overflow threshold. Zero's positive neighbor is 2**-149.
    At a midpoint, the low significand bit decides the ties-to-even winner.
    """
    x = rational(x)
    out = value(bits)
    mag = bits & ~SIGN
    if x and bool(bits & SIGN) != (x < 0):
        return False
    a, y = abs(x), abs(out)
    even = not (mag & 1)
    low = (y + value(mag-1))/2 if mag else F(0)
    high_neighbor = value(mag+1) if mag < MAX_FINITE else pow2(128)
    high = (y + high_neighbor)/2
    return (a > low or (a == low and even)) and (a < high or (a == high and even))


def exact_bits(x) -> int:
    """Require a value already on the runtime lattice; never silently quantize."""
    x = rational(x)
    b = round_bits(x)
    if value(b) != x:
        raise ValueError('runtime operand is not an exact binary32 value')
    return b


@dataclass(frozen=True)
class Operation:
    kind: str
    operands: tuple[int, ...]
    result: int
    defect: F

    def validate(self):
        args = tuple(value(b) for b in self.operands)
        if self.kind == 'add' and len(args) == 2:
            real = args[0] + args[1]
            negzero = self.operands[0] == SIGN and self.operands[1] == SIGN
        elif self.kind == 'mul' and len(args) == 2:
            real = args[0] * args[1]
            negzero = bool((self.operands[0] ^ self.operands[1]) & SIGN)
        else:
            raise ValueError('unknown elementary rounding relation')
        if not rounding_cell_contains(real,self.result):
            raise ValueError('rounding result lies outside its exact nearest-even cell')
        if not real and self.result != (SIGN if negzero else 0):
            raise ValueError('zero sign does not match the elementary operation')
        if self.defect != value(self.result) - real:
            raise ValueError('rounding defect detached from its operation')


@dataclass(frozen=True)
class Normalization:
    input_bits: int
    seed_bits: int
    output_bits: int
    operation_start: int


@dataclass
class Arithmetic:
    """A straight-line operation ledger, with exact finite binary32 values."""
    operations: list[Operation] = field(default_factory=list)
    normalizations: list[Normalization] = field(default_factory=list)

    def add(self, a: int, b: int) -> int:
        z = value(a) + value(b)
        out = round_bits(z, negative_zero=(a == SIGN and b == SIGN))
        self.operations.append(Operation('add', (a,b), out, value(out)-z))
        return out

    def neg(self, a: int) -> int:
        value(a)  # reject a nonfinite encoding before flipping its sign
        return a ^ SIGN

    def sub(self, a: int, b: int) -> int:
        return self.add(a, self.neg(b))

    def mul(self, a: int, b: int) -> int:
        z = value(a) * value(b)
        out = round_bits(z, negative_zero=bool((a ^ b) & SIGN))
        self.operations.append(Operation('mul', (a,b), out, value(out)-z))
        return out

    def invsqrt(self, number: int) -> int:
        """Literal Mahony_AHRS<float>::invSqrt for nonnegative finite words.

        Zero and positive subnormal norm sums are retained, not replaced by an
        ideal inverse square root or silently omitted from the graph.
        """
        value(number)
        if number & SIGN:
            raise ValueError('norm sum must have a nonnegative encoding')
        seed = MAGIC - (number >> 1)
        value(seed)
        first = len(self.operations)
        t = self.mul(self.mul(self.mul(number, exact_bits(F(1,2))), seed), seed)
        out = self.mul(seed, self.sub(exact_bits(F(3,2)), t))
        self.normalizations.append(Normalization(number,seed,out,first))
        return out

    def verify(self):
        for op in self.operations:
            op.validate()
        for witness in self.normalizations:
            a = Arithmetic()
            if witness.seed_bits != MAGIC - (witness.input_bits >> 1):
                raise ValueError('inverse-sqrt seed is detached from the input word')
            if a.invsqrt(witness.input_bits) != witness.output_bits:
                raise ValueError('inverse-sqrt result is detached from the Newton program')
            start = witness.operation_start
            if tuple(self.operations[start:start+5]) != tuple(a.operations):
                raise ValueError('normalization is detached from the same operation ledger')


@dataclass(frozen=True)
class StepResult:
    vertical: object
    operations: tuple[Operation, ...]
    normalizations: tuple[Normalization, ...]
    profile: str = PROFILE


def step_initialized(state, cfg, *, dt, gyro, acc, profile=PROFILE) -> StepResult:
    """Closed initialized branch under the named arithmetic profile.

    Returns the existing finite vertical Result, so WPE/band/stillness consumers
    retain one vertical successor. Inputs must already be binary32. Scalar Eigen
    dot uses p0+(p1+p2); this differs from Mahony's left-associated norm sums.
    Signed-zero differences are forgotten only at the existing rational API.
    """
    from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
    if profile != PROFILE:
        raise ValueError('unqualified compiler/evaluation profile')
    if not isinstance(state,V.State) or not isinstance(cfg,V.Config):
        raise TypeError('existing finite private-observer state/config required')
    if not state.initialized:
        raise ValueError('startup FromTwoVectors seed remains a separate obligation')
    if len(gyro) != 3 or len(acc) != 3:
        raise ValueError('three-dimensional body-frame packet required')
    if rational(dt) <= 0:
        raise ValueError('positive runtime dt required')
    q = tuple(exact_bits(x) for x in state.q)
    integ = list(exact_bits(x) for x in state.integral)
    g = list(exact_bits(x) for x in gyro)
    raw = tuple(exact_bits(x) for x in acc)
    h, kp, ki, gravity = map(exact_bits,(dt,cfg.two_kp,cfg.two_ki,cfg.gravity))
    elapsed = exact_bits(state.elapsed)
    exact_bits(state.up)
    a = Arithmetic()
    add, sub, mul, neg = a.add, a.sub, a.mul, a.neg
    half, two = exact_bits(F(1,2)), exact_bits(2)
    x,y,z = tuple(neg(x) for x in raw)
    if any(value(x) != 0 for x in raw):
        norm = add(add(mul(x,x),mul(y,y)),mul(z,z))
        r = a.invsqrt(norm)
        x,y,z = mul(x,r),mul(y,r),mul(z,r)
        q0,q1,q2,q3 = q
        vx = sub(mul(q1,q3),mul(q0,q2))
        vy = add(mul(q0,q1),mul(q2,q3))
        vz = mul(half,add(sub(sub(mul(q0,q0),mul(q1,q1)),mul(q2,q2)),mul(q3,q3)))
        error = (sub(mul(y,vz),mul(z,vy)),sub(mul(z,vx),mul(x,vz)),sub(mul(x,vy),mul(y,vx)))
        if value(ki) > 0:
            integ = [add(i,mul(mul(ki,e),h)) for i,e in zip(integ,error)]
            g = [add(v,i) for v,i in zip(g,integ)]
        else:
            integ = [0,0,0]
        g = [add(v,mul(kp,e)) for v,e in zip(g,error)]
    # Exactly zero accelerometer skips feedback AND clearing/applying integral.
    gx,gy,gz = (mul(v,mul(half,h)) for v in g)
    q0,q1,q2,q3 = q
    un = (add(q0,sub(sub(mul(neg(q1),gx),mul(q2,gy)),mul(q3,gz))),
          add(q1,sub(add(mul(q0,gx),mul(q2,gz)),mul(q3,gy))),
          add(q2,add(sub(mul(q0,gy),mul(q1,gz)),mul(q3,gx))),
          add(q3,sub(add(mul(q0,gz),mul(q1,gy)),mul(q2,gx))))
    norm = add(add(add(mul(un[0],un[0]),mul(un[1],un[1])),mul(un[2],un[2])),mul(un[3],un[3]))
    r = a.invsqrt(norm)
    qn = tuple(mul(x,r) for x in un)
    w,x,y,z = qn
    down = (mul(two,sub(mul(x,z),mul(w,y))),mul(two,add(mul(y,z),mul(w,x))),
            add(sub(sub(mul(w,w),mul(x,x)),mul(y,y)),mul(z,z)))
    products = tuple(mul(x,y) for x,y in zip(down,raw))
    up = neg(add(add(products[0],add(products[1],products[2])),gravity))
    elapsed = add(elapsed,h)
    nxt = V.State(tuple(map(value,qn)),tuple(map(value,integ)),True,value(elapsed),value(up))
    return StepResult(V.Result(nxt,value(up)),tuple(a.operations),tuple(a.normalizations))


def readiness():
    return {
        'initialized_profile_specific_finite_binary32_graph': True,
        'inverse_sqrt_seed_and_Newton_operands_bound': True,
        'exact_operation_rounding_defects_retained': True,
        'independent_reciprocal_witnesses_used': False,
        'actual_target_compiler_profile_qualified': False,
        'startup_seed_qualified': False,
        'nonfinite_intermediate_branches_qualified': False,
        'disturbance_dependent_uniform_defect_bound_proved': False,
        'source_uniform_word_qualified': False,
        'complete_word_finite_identity': False,
        'ALT_LIVE_PASS': False,
        'ALT_STARTUP_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
