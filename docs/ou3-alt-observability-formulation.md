# OU-III ALT corrected observability formulation

## Selected architecture

The old requirement of strict coercive joint24 contraction on every admitted 3 s word is retired. PR #533 proves that an admitted quiet ungauged word contains the exact unipotent pair `(theta_z,bg_z)` with spectral radius one. This falsifies that theorem formulation, not filter stability.

ALT now uses two regimes. During ungauged intervals the proof works modulo the unobservable heading/axial-bias centre dynamics and must bound finite-prefix growth without claiming absolute-heading contraction. A yaw-only quotient is insufficient because `bg_z` remains neutral. The quotient is legal only when the actual map descends to it: `Q_after^T A N_before = 0`; deleting coordinates is not a proof of equivariance.

After magnetic gauge is established, the full-state theorem is conditional on an informative magnetic-service class. Attempted calls at gaps <=40 ms are retained runtime facts but are not themselves accepted informative observations. The conditional Normal-Live theorem now explicitly assumes a bound on gaps between accepted informative events and a positive transported heading/bias information floor over every overlapping superword. Membership is an additional hypothesis for that regime, not a consequence of attempted calls and not a restriction on the unconditional source theorem. See `ou3-alt-informative-service-theorem.md`.

## Controlling measurement

Before interval enclosure, measure on the actual same-history carried superword

`rho_W = lambda_max(Z^T A_W^T M_after A_W Z, Z^T M_before Z)`.

`Z` injects the canonical unsupplied coordinates; all output rows of `A_W Z` remain, including A21 bias corrections. A common metric must be the same `M` on every word. Compatible endpoint metrics must come from one declared storage law with uniform coercivity; choosing a separate favorable metric for each word is invalid.

The gauged motion-block spectra near 0.9964 reported by PR #533 are feasibility evidence only. They are not an upper bound on this storage ratio and their roughly 3.6e-3 distance from one is not yet an interval-enclosure budget.

## Magnetic service

For accepted gauged magnetic events in a superword, transport the actual whitened residual sensitivity to the same initial heading/axial-bias coordinates and form `G_B = sum B_i^T B_i`. The selected conditional service class is parameterized by explicit constants `T_B` and `alpha_B>0` such that informative-event gaps are <=`T_B` and `G_B >= alpha_B I_2` on the declared window. A deployment must declare these constants and assume or prove the resulting accepted-service property. The proof does not infer constants from callback cadence or fit them to diagnostic extrema. Rejected calls, ungauged calls, and heading-degenerate rows provide no such service.

## Carried native storage feasibility

The prescribed endpoint law is `M(P)=diag(P21^-1,I3)`. It uses the **full**
carried 21-state covariance in both modes, including held-bias uncertainty;
it never inverts an H18 marginal or drops active-bias output energy. No metric
is fitted to a word. The source is the admitted stationary COMPLETE-BRMM member
`p=v=a=S=0`, identity physical attitude, zero BIAS0, and `B_W=(32,0,0) uT`.
Its field norm and horizontal component meet MAG-BMM150-DET-v1. Raw acceleration
is the binary32 representation of gravity; its rounding residual is within the
commissioned sensor bounds. Zero physical motion is retained as a necessary
member of the source class, not substituted for the full wave theorem.

`carried_storage_rho_diagnostic.py` runs the actual wrapper from construction,
with default tuning and 25 Hz attempted magnetic calls. The public external
bias hold selects H18; A21 uses actual magnetic unlocking; a third history
releases the external hold after the first measured window. The observed root
is sample 20,048 (about 100.24 s), after actual Live/north/unlock ancestry.
No filter state, covariance, gain, clock, or tuner snapshot is installed.

Three consecutive 600-sample windows per history retain all native memory:

| History | First window | Second window | Third window |
| --- | ---: | ---: | ---: |
| H18 held | 0.9760520113 | 0.9765933549 | **0.9771585509** |
| A21 active | 0.9760555562 | 0.9765928606 | 0.9771579208 |
| H18 then release to A21 | 0.9760520113 | 0.9765963647 | 0.9771580464 |

These are projected **local tangent storage ratios**, not finite-error or
source-uniform rho. The probe requires exactly zero nominal state and
innovations: only on this branch are prediction factors and `I-KH` sufficient
without nonzero-residual gain/reset terms. It retains actual binary32 operands,
full native covariance updates/floors, the actual R_S and R_acc, and the release
covariance operation. The unchanged uninstrumented build produces bit-identical
sample states. Derivatives of discontinuous floating-point rounding are not
claimed; arithmetic defects remain a separate theorem obligation.

Every window contains 75 accepted native magnetic innovations. Their
information rows are `R_mag^-1/2 H_mag Phi_before E_(theta_z,bg_z)`, using the
actual H, R and preceding event product. The minimum pair-information eigenvalue
is above 26,324 in these windows; it is a measured value, not an assumed uniform
floor. The earlier `[1,t]` proxy supplies no runtime qualification.

The worst storage margin is `1-rho=0.0228414491`; the corresponding amplitude
margin `1-sqrt(rho)=0.0114866966` is a **point** sensitivity budget. The limiting
direction combines roll-axis gyro bias with lateral motion/primitive errors.
The same words have identity-metric ratios around 20.4–20.8. Thus failure of
that one metric does not falsify compatible storage.

A separate extended-precision composition and 60-decimal-digit terminal
inverse/eigensolve agree with binary64 within 3e-15. This is a conditioning
check, not an outward enclosure of the word. The report includes the full
maximizing direction and an energy-change decomposition by actual operation.
Its largest decrement is the accelerometer event; magnetic service constrains
the heading pair but is not the limiting direction in this quiet history.

Results and source hashes are in
`reports/results/ou3_alt_storage/carried-quiet-diagnostic.json`. Reproduce with:

```sh
PYTHONPATH=.:tools/stability OPENBLAS_NUM_THREADS=1 python3 -m \
  tools.stability.ou3_alt_contraction.carried_storage_rho_diagnostic \
  --work /tmp/ou3-alt-carried-storage --output /tmp/ou3-alt-carried-storage.json
```

The remaining controlling question is whether this same prescribed law can
satisfy a uniform finite-word inequality over the actual carried source family.
It requires uniform bounds `0<p_min I <= P21 <= p_max I`, all informative-service
histories, nonzero finite errors, nonlinear resets/projection, BIAS1/BIAS2 and
machine supplies, with the complete finite-master prerequisites retained.
Even identity dynamics can have a subunit inverse-covariance ratio if covariance
grows without bound; an executable negative control protects this distinction.
The nine quiet windows justify keeping this storage law as a candidate. They
do **not** authorize interval refinement or establish eventual startup for
other admitted histories. Do not spend effort fitting metrics to these traces
or relabeling their extrema as universal bounds.

## Non-promotion

This selection does not certify source-uniform magnetic acceptance, quotient equivariance, a common/compatible storage, a source-uniform rho, interval enclosure, startup, Live, or end-to-end stability. `storage_search_allowed`, `ALT_STARTUP_PASS`, `ALT_LIVE_PASS`, and `ALT_END_TO_END_PASS` remain false. The finite-master inventory retains six open and five component-closed qualifications; the corrected formulation does not silently close any of them.

The independent P2/P3/P4/P5 proof track is unchanged, including frozen P3 `1e-18`.
