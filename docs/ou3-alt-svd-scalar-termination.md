# Scalar Eigen seed: uniform bounds under review

`finite_seed_eigen_svd.py` executes the pinned Eigen 3.4.0 dependency graph.
Its named profile uses scalar binary32 RNE, gradual underflow, correctly rounded
sqrt/division, and nontrapping IEEE special-value arithmetic. It computes both
Householder reflectors, Q, the 2×2 work matrix, every Jacobi rotation and the
actual loop predicate. Native comparisons are regression evidence; they do not
establish the uniform bounds below. Universal termination remains unpromoted
until this argument passes independent review.

Let `u=2^-24`, `eta=2^-150`, and `gamma(n)=nu/(1-nu)`. The two preceding
accelerometer normalizations give

```
v1 = e3,
1-10^-6 < ||v0||² < 1+10^-6,
v0.z < RN(-1+RN(10^-5)).
```

These imply `-1.000001 < v0.z < -0.99998`, transverse norm below `.005`, and
`||v0+e3|| < .006`. Scaling by the largest input component preserves these
bounds with room to use `.007` for the last quantity. Neither pivot is excluded.

## QR totality and a coarse axis bound

The first pivot has norm between `.99` and `1.01`, and its last component has
magnitude above `.99`. Thus its Householder norm root and cancellation-free
denominator have magnitude above `.99`. Bounding the construction's component
errors and the reflector application by `256u` gives

```
work = [[a,0],[b,c]],
.99 < |a| < 1.01,   .98 < |b| < 1.02,   a*b < 0,   |c| < .02.
```

For the second reflector, a tail-square result at or below the smallest normal
selects Eigen's exact identity reflector. Otherwise the tail-square exceeds
`2^-126`. Its exact input norm squared is at least `2^-126(1-2u)`.
Square/sum/root rounding therefore gives
`beta² >= (1-8u)||x||²`; the normal sqrt output itself has no underflow charge.
The signed beta choice makes `|RN(x0-beta)| >= (1-u)|beta|`. Component division
and tau evaluation give

```
||essential||² <= 1+32u,
0 <= tau <= 2+32u.
```

These bounds also cover the first reflector. For the exact matrix represented
by the stored reflector, its operator norm is at most
`(2+32u)²-1`. The explicit scalar application contributes less than `64u` times
the input norm, plus its underflow charges; use `3+192u` for the total factor.
Two applications to e3 give axis norm below 10, hence squared norm below 100.
This deliberately coarse bound suffices to establish finite Mahony seeding;
it does not claim a sharp nullspace or orientation error. All QR scalar
intermediates are finite. The omitted norm downdates have two statically bounded
iterations and cannot affect either the second pivot index or Q.

## Proposed two-sweep bound

For the initial work matrix above, Eigen extracts the reversed block
`[[c,b],[0,a]]`. Its first polar ratio is `(a+c)/(-b)`, between `.95` and `1.06`
and positive. The polar coefficients are therefore bounded away from zero.
The symmetrized block has a large second diagonal and a small first diagonal.
The stable Jacobi tangent formula preserves this order.

The proposed first-sweep ledger uses coefficient errors below `32u`, error
below `128u` in the symmetrization identity, and below `128u` in the composed
left rotation. The two matrix applications and asymmetry defect together fit
within `1024u`. It gives, in the original work-matrix order,

```
|work01|, |work10| <= E = 1024u,
1.3 <= |work00| <= 2,
|work11| <= .03.
```

The first two inequalities also follow from the approximately orthogonal
left/right transformations, the initial norm/determinant, and the small
offdiagonal remainder; the diagonal order must be retained in this argument.

If Eigen requests a second sweep, its threshold is at least `3u`: maxDiag is
already at least `.99` and the precision multiplier is `4u`. With one active
offdiagonal above `3u`, any nonzero rounded difference of the offdiagonals is
at least `2^-48`. Thus the second polar quotient cannot overflow. Its exact
cosine has magnitude near one and its sine is below `2E`.

Use the symmetrization identity and the *lower* offdiagonal to bound the
evaluated upper offdiagonal by `2E`. A separate interval bound on two canceling
products would be too large. The second symmetrized diagonal gap exceeds one.
For `|y| >= Y=2^-40`, the Jacobi tangent graph has no overflow and the following
weighted roundoff ledger is proposed:

| Quantity | Absolute error bound |
| --- | ---: |
| Polar cosine | `32u` |
| Polar sine | `128uE` |
| Jacobi cosine | `32u` |
| Jacobi sine | `128uE` |
| Composed left cosine | `128u` |
| Composed left sine | `1024uE` |
| Left-applied diagonal | `512u` |
| Left-applied offdiagonal | `4096uE` |
| Final offdiagonal, including asymmetry | `16384uE` |

The small sine bounds are essential: multiplying a generic coefficient error
by a full diagonal would lose termination. The elementary errors must be
charged on the actual scalar operands and signs in the implemented graph.

For `|y| < Y`, the calculation may overflow in `tau²`, despite finite input.
When this gives infinite w, the cancellation-free tangent denominator gives
signed zero t. The remaining cosine/sine are finite. The finite-tau case has
`|t| <= 4|y|`; either case is covered by an additional `256Y` bound. No
nonfinite branch is deleted from the implemented graph.

If every row of this ledger is verified, the second-sweep residual satisfies

```
16384u(1024u) + 256*2^-40 < 3u,
```

so the next loop predicate is false. The bound would prove at most two sweeps
for every admitted normalized source input, including rank-one inputs. At
present the exact operation graph is complete, but this uniform error ledger
is still under review; it must not be replaced by the native sample result.
