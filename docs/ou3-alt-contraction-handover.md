# OU-III parallel ALT proof handover

## Resume point and independent scope

Read `AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover,
`docs/ou3-alt-finite-measurement-proof.md`, `docs/ou3-alt-finite-core-composition.md`,
`docs/ou3-alt-runtime-primitives.md`, `docs/ou3-alt-contraction.md`, the ALT section
of `docs/ou3-proof-research-state.md`, and `docs/ou3-brmm-main-handover.md`.
Continue the open ALT PR on its branch; after it is merged, start a new PR from
latest main.

The original P2/P3/P4/P5 route remains independently continuable and unchanged.
Its PASS labels do not discharge ALT obligations. The shipping implementation,
P3 threshold `1e-18`, physical source and original gates are not modified here.

## Target and immutable physical contracts

Retain `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)` in R^24, all motion/bias cross
terms, and coercive storage satisfying

`V_next <= rho V + w'Gamma w + c_supply'Beta c_supply`, `0 < rho < 1`.

Held H18 bias and nonrelaxing BIAS2 truth may use independently justified
bounded supply. Unknown motion errors may not be renamed as bounded inputs.
Keep physical forcing and correlated state/source ports inside the finite graph.

Corrected COMPLETE-BRMM retains the padded acceleration bound 8.8 m/s^2,
all-time centered primitive bound D_S <= 1100 m*s, one-time Live S origin,
zero lever arm, actual anisotropic R_S, same-signal WPE/bandpass/sigma/tau/T_S
ancestry, staged tuner commits, scheduler guards, and separate BIAS0/1/2
physical histories. Fresh centered e_S=0 is not an arbitrary 300 m*s entry ball.
The remaining hard entry coordinates must come from actual startup and truth,
not filter zero initialization or covariance consistency.

## Proved finite algebra

The finite measurement note supplies these conditional real-arithmetic identities:

- finite inverse-free measurement, residual secants, both quaternion injection
  branches and same-beta radial projection on full joint24;
- physical prediction with the actual physical rotation increment, nonzero
  angular/model defect, correlated q15 translation moments, gyro-bias drift,
  and one driver shared by e_ba and beta;
- masked-branch Joseph cancellation using `K Sigma=N`, without assuming `N=P H'`.

`finite_core.py` carries finite mean/error, nominal quaternion, current full P,
accepted residual/gain/Joseph/reset/projection and the persistent physical
reference through consecutive represented prefixes. Prediction reuses the same
`PhysicalSegment`; no Jacobian cocycle, fresh covariance or new Live origin is
inserted.

`finite_prediction_covariance.py` proves the ordinary pre-floor/pre-S-service
shipping covariance update as the full block congruence

`F=diag(F_AA,F_LL,phi_b I3)`, `Q=diag(Q_AA,Q_LL,Q_BB)`, `P+=F P F'+Q`,

retaining every AA/linear/BA cross covariance. `F_LL` uses the same coefficient
tuple as the finite mean and the correlated/independent Q_LL assembly follows
the shipping block layout.

`finite_ou_runtime_primitives.py` removes free mean coefficients, free active
`phi_hat`, and free active `Q_BB` at the higher-level paired runtime entry. One
`(h,tau,alpha)` root generates `phi_va`, the literal small/general
`safe_phi_A_coeffs` branch, and `F_LL`. One active `(tau_b,phi_b,Q_bacc)` root
generates both the BA mean factor and
`Q_BB=Q_bacc*tau_b*(1-phi_b^2)/2`; H18 remains exactly `phi_b=1,Q_BB=0`.
The source relation still must prove the deployed exponential/tuner ancestry of
those scalar roots. This is not source admission or finite-precision closure.

H18 retains the full latent BA covariance in accelerometer innovation. The
error/storage state is never reduced to 18 coordinates. A nonzero numerical
solve defect is rejected rather than dropped.

## Implementation correspondence, not source qualification

`shipping_finite_identity.py` builds the actual wrapper twice from startup.
One temporary include overlay adds passive observations at selected core
boundaries; the other build is uninstrumented. Observation calls must erase to
the original source and recorded sample states must match bit-for-bit. No runtime
state, covariance, gain, schedule or mode is forced by the harness.

Three 600-step windows exercise H18, A21 and an actual H18->A21 release. The
audit checks prediction, covariance/floor, physical S residuals, latent H18 BA,
masked numerator/full innovation, Joseph, reset, same-beta projection, due/not-
due S scheduling and the observed release. This is implementation correspondence,
not a universal source cover, branch proof, contraction result or roundoff bound.

Run the isolated finite checks with:

```sh
PYTHONPATH="$PWD:$PWD/tools/stability:$PWD/tests/ou3_alt_contraction" \
  python3 -m unittest test_finite_physical_prediction \
  test_finite_measurement_graph test_finite_covariance_rank3 test_finite_core \
  test_finite_prediction_covariance test_finite_ou_runtime_primitives \
  test_shipping_finite_identity test_core test_proof_plan test_bias_families -v
```

The inherited broad ALT suite remains separate and may fail on the unchanged
shared Mahony prerequisite. Do not weaken that premise to make ALT green.

## Decisive missing object and next work

**The complete source-uniform finite 600-step word is still NOT materialized.**
The controlling gap is now narrower but remains C/E: finite representation and
same-history runtime/source attachment.

Next bind the actual attitude/noise prediction primitives (`F_AA,Q_AA`) including
rotation/B construction, fast/structured Q branch, Simpson evaluations and PSD
hygiene. Bind analytic Qaxis generation/regularization and the scalar exponential
roots to the same tuner/runtime predecessor. Then compose pending a_w covariance
positive-part inflation, final covariance hygiene, periodic S service, applied
R values, safe-LDLT accept/retry/reject branches, frontend/tuner staged commits,
asynchronous magnetometer state and every H18/A21 edge. Expose every literal
prefix rather than selected observation points.

Only after that finite master closes may the high-precision feasibility
diagnostic and first common joint24 storage search run. Then come rigorous source
cover, endpoint dissipativity, every-prefix retention, the physics-compatible
Live basin, hybrid landing, aggregate fresh entry, startup capture and finite
precision.

## Gate status and failure classification

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Original `P4_PASS` / `P5_MAY_START` remain untouched by ALT.

The current gap is C/E, not evidence of filter instability or a failed common
metric. There is no certified rho, ultimate bound, retained basin or capture
time. Do not substitute more seeds, thinner boxes, a frozen observer, wordwise
S reset, position reanchor or a smaller entry set for the missing finite source
relation.
