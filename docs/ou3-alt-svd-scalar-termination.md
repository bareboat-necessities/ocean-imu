# Eigen seed: uniform termination and scalar producer

`finite_seed_eigen_svd.py` executes the pinned Arduino Eigen 0.3.2 / Eigen
3.4.0 dependency graph for the actual real 2×3 `JacobiSVD` seed. It computes
the scale, both QR pivots, two Householder reflectors, Q, work matrix, Jacobi
rotations and retained loop predicate. The requested axis is always computed.
A legacy axis witness is accepted only when equal to that result.

`finite_seed_svd_roundoff.py` proves that every normalized near-antiparallel
source input returns after at most two Jacobi sweeps. The exact rational
argument has undergone independent review. Native bit comparisons check the
implementation; they do not establish the quantified result.

The proof uses binary32 round-to-nearest-even, gradual underflow, correctly
rounded division/sqrt and nontrapping IEEE arithmetic. It also covers local
multiply-add/subtract contractions within the same expression tree: deleting
an internal rounding selects zero from its allowed error interval. It does
not permit reassociation or cross-statement fusion. The exact Python producer
computes the unfused trace. Coverage of the local-contraction family does not
assert equality between its trace bits and a compiled target trace; target
compiler/ABI qualification is a separate obligation.

## Source and QR domain

Put `u=2^-24`, `eta=2^-150`. A raw accelerometer component cap of 160 and the
actual first norm guard `norm>.001f` make both source normalizations finite.
For exact `S=sum(x_i²)`, the nonnegative dot product lies between
`(1-u)^3*S-5eta` and `(1+u)^3*S+5eta`. Normal sqrt and component-division
bounds then establish

```
v1 = e3,
1-10^-6 < ||v0||² < 1+10^-6,
v0.z < RN(-1+RN(10^-5)).
```

After scaling, both column norms lie in `.99..1.01`, their last components
have magnitude above `.99`, their transverse norms are below `.005`, and
the norm of their sum is below `.007`. Both pivot outcomes are retained.

The first reflector uses a cancellation-free denominator. Operation bounds
for beta, its denominator, essential vector and tau give an application norm
error below `256u`, establishing

```
work = [[a,0],[b,c]],
.99 < |a| < 1.01,  .98 < |b| < 1.02,  a*b < 0,  |c| < .02.
```

For either reflector, a tail-square at or below `MIN_NORMAL` selects literal
identity. Otherwise the exact input squared norm exceeds
`MIN_NORMAL*(1-8u)` and the squared beta ratio lies in `1±16u`.
The signed beta prevents cancellation in `x0-beta`. Thus
`||essential||²<1+32u` and `0<=tau<2+32u`, including tiny nonzero tails.
The stored reflector operator norm is at most `(2+32u)²-1`.

`Bound` propagates the actual multiply/add node errors. Each applied output
component has error below `32u*X+32eta` for input norm `X<=4`; the vector error
is below `64u*X+64eta`. The error expression is affine in X, so the exact
endpoint checks establish the entire interval. The effective reflector
factor is below `3+256u`. Two applications to e3 give squared axis norm below
`81.002<100`. This bound suffices for finite seeding without a sharp nullspace
or orientation estimate.

`finite_seed_svd_axis_reduction.py` proves the write-index invariant: all
post-QR Jacobi rotations and sorting swaps touch V columns 0 and 1. V column 2
is the unchanged QR column on every returning rank-one or rank-two run.

## Two-sweep induction

Eigen extracts the reversed block `[[c,b],[0,a]]`. Its first polar ratio is
positive and bounded away from zero; both exact coefficients exceed `.68`.
The symmetrized diagonal gap exceeds `1.3`, so the stable Jacobi tangent is
below `.05`. The exact eigenvalue shifts and propagated diagonal errors give

```
|work01|, |work10| < E = 1024u,
1.3 < |work00| < 2,
|work11| < .03.
```

All coefficient error bounds are derived from the stable formula's individual
rounding factors. Matrix multiplication uses exact magnitude identities and
retains the same stored-work reference. In particular, the ideal polar output
is symmetric; its lower offdiagonal supplies the small upper-offdiagonal bound
without independently enclosing two canceling large terms.

If a second sweep is active, its retained `maxDiag` is at least the initial
`.99`, so its threshold exceeds `3u`. A nonzero rounded difference of active
offdiagonals has magnitude at least `2^-48`; the polar quotient and its square
cannot overflow. A zero difference selects literal identity, and the stored
matrix is already symmetric. The ideal polar sine is below `2E`.

The second symmetrized diagonal gap exceeds one. Its evaluated offdiagonal
is below `2E`. For `|y|>=Y=2^-40`, all stable Jacobi coefficient operations are
normal and finite. Their relative bounds yield absolute sine errors of order
`uE`. `Bound` propagates these weighted errors through the actual composition
and matrix-application nodes. The discrepancy between ideal polar symmetry
and the symmetric matrix used by the exact Jacobi rotation is charged as
`2*offdiag_error+8E*diag_error`, again of order `uE`.

For `|y|<Y`, the guard can select identity. Otherwise an overflowing `tau²`
produces infinite w and signed-zero tangent, also identity. In the remaining
finite branch the stable same-sign denominator gives
`|t_machine| <= (1+u)*2Y/(1-u)^2+eta <3Y`; the ideal tangent is below Y.
The resulting absolute coefficient error is below `8Y`. Reapplying the same
node graph proves an additional residual below `256Y`.

The exact rational second-sweep residual is below `0.045402u`, while the
retained next threshold exceeds `3.959999u`. Therefore the next predicate is
false. All discarded special-value branches remain in the executable graph;
no fixed iteration cutoff replaces the shipping loop.
