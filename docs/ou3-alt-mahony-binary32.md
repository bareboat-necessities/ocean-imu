# Initialized private-Mahony binary32 graph

## Place in the ALT word

The controlling object is the actual finite same-history map `F_word`, before
any search for `V(F_word(z),xi_next) <= rho V(z,xi) + supply`.
The private observer supplies the one vertical acceleration used by WPE, the
sigma band and the tracker-input filter. Its normalization cannot remain an
arbitrary positive reciprocal when qualifying `F_word`: changing that value
changes the nominal quaternion, vertical signal, tuner and later Kalman word.

`finite_binary32_mahony.py` removes that free choice on the initialized branch
under the explicit arithmetic profile below. This is representation work, not
a tighter local stability bound or a claim about the sign of a storage margin.
No meaningful change in rho can be computed until the remaining word closes.

## Arithmetic profile and supplying source

The profile `binary32-rne-gradual-no-fma-eigen3-scalar` means:

- IEEE binary32, round-to-nearest/ties-to-even at each source basic operation;
- gradual underflow, no reassociation, excess precision or fused contraction;
- left-associated scalar sums in `Mahony_AHRS<float>::update`;
- the scalar Eigen three-term dot reduction `p0+(p1+p2)` for the vertical readout;
- no floating-point traps affecting state evolution from the unused Euler outputs.

The supplying code is `src/ahrs/Mahony_AHRS.h::update/invSqrt` and
`src/tuner/VerticalAccelComplementary.h::update`. The ordinary float observer,
not its double specialization or an ideally normalized quaternion, is modeled.
The target ESP32 compiler/library evaluation profile is **not** qualified by a
host compilation. Other reduction/FMA profiles must be attached separately.

Only finite, initialized inputs with positive dt are accepted. An overflow or
nonfinite intermediate raises an explicit exception rather than being replaced
by a finite point, a clipped value or an identity event. This domain is the
conditional lemma's domain, NOT a smaller admitted BRMM source. Proving that the
actual word stays on the domain, or representing the remaining branches, is
still necessary. The first `FromTwoVectors` seed is not supplied here.

## Lemma 1: exact nearest-even rounding

Write a finite positive binary32 value as `m*2^s`, with `m` an integer. Normal
numbers have `2^23 <= m < 2^24`; subnormals use `s=-149` and `0 <= m < 2^23`.
For an exact rational `x != 0`, choose

```
e = floor(log2(abs(x))),
s = max(-149, e-23),
q,r = divmod(numerator(abs(x)/2^s), denominator(abs(x)/2^s)).
```

Round q upward precisely when twice the remainder exceeds the denominator, or
at equality when q is odd. A carry to `2^24` renormalizes by incrementing s.
The exponent reaching 255 is the explicit overflow boundary. The sign is
retained, including negative underflow to zero. Exact cancellation and the
signs of zero multiplication follow the separately stated operation rules.

**Proof.** Inside each binade the binary32 lattice spacing is `2^s`; the two
nearest candidates therefore bracket `abs(x)/2^s` at q and q+1. Comparing twice
the remainder to the denominator selects the nearer integer without a host
floating-point calculation. At a midpoint the low significand bit selects the
even integer. At a binade boundary the carry is the common representable upper
endpoint; at the subnormal boundary the fixed spacing makes the same argument
apply, including the zero/minimum-subnormal midpoint. The maximum finite
number's virtual upper neighbor is `2^128`, giving the usual midpoint overflow
boundary. These cases exhaust the finite lattice. Negation is symmetric.

`round_bits` implements this construction using integers and Fractions only.
`rounding_cell_contains` independently checks membership in the midpoint cell of
the returned binary32 word, with ties admitted exactly for even significands.
It does not call `round_bits`. An elementary operation is accepted by its ledger
checker only when that independent cell relation and its exact defect agree.

## Lemma 2: the literal fast inverse-square-root relation

For a nonnegative finite norm word b, let

```
j = 0x5f375a86 - floor(b/2),
y0 = decode_binary32(j),
t1 = RN32(x * 1/2),
t2 = RN32(t1 * y0),
t3 = RN32(t2 * y0),
t4 = RN32(3/2 - t3),
y  = RN32(y0 * t4).
```

The graph stores b, j, y and the five basic-operation rows together. It never
accepts a caller-supplied y. Zero and positive subnormal norm inputs are kept;
in particular, zero is **not** interpreted using ideal `1/sqrt(0)` semantics.

**Proof.** A nonnegative finite float has sign bit zero, so its `int32_t` bit
representation is nonnegative and the shipping right shift equals floor(b/2).
The seed integer and reinterpretation are therefore exact aliases of the same
input word. The five rounding relations then follow Lemma 1 in the source's
parenthesized order. Induction over them gives the shipping Newton result on
the declared profile. This is an all-input program identity wherever the
finite graph returns; a list of sampled normalizer values is not its proof.

The earlier read-only `ou3_fast_inv_sqrt_interval.py` remains a separate
outward-enclosure primitive. Its existing scalar bounds are not replaced or
promoted here. The new object supplies the exact input/output/defect relation
needed by the finite runtime graph, not another independently selected shell.

## Lemma 3: initialized observer successor

For every initialized predecessor and raw gyro/accelerometer packet in the
stated arithmetic domain, `step_initialized` computes the same decoded
quaternion, integral feedback, elapsed time and vertical output as the supplying
shipping observer under the declared profile.

**Proof.** The accelerometer validity predicate compares the same three packet
coordinates to zero. On a nonzero packet Lemma 2 determines its normalization
from the same rounded squared-norm sum. All gravity half-vector and cross-error
products then have the same operands and evaluation order. Positive Ki advances
and applies the new integral; nonpositive Ki clears it on this nonzero-accel
branch. Exactly zero acceleration skips both integral clearing and application,
as shipping does. The gyro feedback and Euler quaternion increment thus agree
operation by operation. Lemma 2 then determines normalization of that SAME
rounded quaternion increment. The final component multiplies, unnormalized
third-row expression and scalar Eigen dot produce the same up acceleration.
The elapsed-time addition uses binary32 too. Induction over this ordered list
proves the equality; branch predicates are not freely supplied witnesses.

For every elementary operation the ledger retains

`defect_j = decode(result_j) - exact_operation_j(decode(operands_j))`.

It is an exact finite identity with endogenous defects. No finite collection
of these realized defects is a source-uniform disturbance bound, and unknown
motion errors have not been relabeled as bounded inputs. Signed-zero encodings
are kept inside the ledger and are forgotten only by the existing rational
`V.State/V.Result` boundary; the latter promises numerical, not signed-zero
bitwise, equality. Unused Euler outputs are not included in that promise.

## Composition and validation boundary

`finite_vertical_complementary_runtime.step(..., arithmetic_profile=PROFILE)`
returns the existing `V.Result`. The raw-packet bridge already forwards this
argument through `vertical_step_from_raw`, so the private observer still
consumes the raw B-frame packet, not the de-heeled MEKF packet. Existing
WPE/band consumers can use that same result without a second vertical input.
The profile path rejects every supplied seed/reciprocal witness and rejects
uninitialized state; it never falls back silently to the conditional real map.
The default real-arithmetic helper remains available for its conditional lemmas.

The standalone regression compiles the actual headers and only observes private
state while public updates advance from ordinary startup. It does not install
synthetic estimator roots or fit a trajectory/metric. It checks normalizer
boundary cases in every finite exponent bin, zero/subnormal cases, both Ki
branches, a public gain change, and zero-acceleration feedback skipping. These
are correspondence tests; none admits its packets as COMPLETE-BRMM histories.

Still open: target-profile qualification, initial seeding, overflow/nonfinite
branches, every upstream/downstream transcendental/solver binding, a uniform
same-history defect enclosure, source admission, the complete 600-step word,
compatible storage, retention and startup capture. All three ALT gates remain
false. The independent original P2/P3/P4/P5 route is not changed by this lemma.
