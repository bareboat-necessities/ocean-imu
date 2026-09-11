# OU-III parallel ALT proof handover

## Resume point and independent scope

Read `AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover,
`docs/ou3-alt-finite-measurement-proof.md`, `docs/ou3-alt-contraction.md`, the
ALT section of `docs/ou3-proof-research-state.md`, and
`docs/ou3-brmm-main-handover.md`. Continue the open ALT PR on its branch; after
it is merged, start a new PR from latest main.

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

The finite measurement note is the supplying derivation for these conditional
real-arithmetic identities:

- finite inverse-free measurement, residual secants, both quaternion injection
  branches and same-beta radial projection on full joint24;
- physical prediction with the actual physical rotation increment, nonzero
  angular/model defect, correlated q15 translation moments, gyro-bias drift,
  and one driver shared by e_ba and beta;
- masked-branch Joseph cancellation using K Sigma=N, without assuming N=P H'.

`finite_physical_prediction.py` is the physical prediction implementation.
Its attitude and translation identities have all-coefficient polynomial checks.
The older `finite_prediction_graph.py` / `finite_prediction_deployed_step.py`
attitude routines describe a sampled shadow; they do not replace a continuous
physical increment or prove its defect vanishes. Their retained algebra callers
are not a physical-word supplying lemma.

H18 still needs the full latent BA covariance in accelerometer innovation.
On the held zero-cross-block invariant the reduced covariance uses R_acc+B0,
not bare R_acc. The error/storage state is never reduced to 18 coordinates.
The checked rank-three Joseph evaluator rejects a nonzero solve defect rather
than dropping its E K' contribution. No deployment roundoff bound is supplied.

## Implementation correspondence, not source qualification

`shipping_finite_identity.py` builds the actual wrapper twice from startup.
One temporary include overlay adds passive observations at selected core
boundaries; the other build is uninstrumented. Observation calls must erase to
the original source and the recorded sample states must match bit-for-bit.
No runtime state, covariance, gain, schedule or mode is forced by the harness.

The regression checks three 600-step windows: H18, A21 and an actual H18->A21
release. It checks prediction, covariance/floor, physical S residuals, full
innovation/numerator, Joseph, finite reset, same-beta projection, due/not-due S
scheduling, and consecutive core ancestry. It retains one Live origin and
fresh centered e_S=0. The audit records selected frontend fields, not a complete
frontend transition relation. Unseen repair/rejection, active radial projection,
watchdog and other guards remain explicitly unqualified. A test trajectory is
not a universal physical source cover or startup/capture proof.

Run the isolated checks without the old frozen-word observer:

```sh
PYTHONPATH="$PWD:$PWD/tools/stability:$PWD/tests/ou3_alt_contraction" \
  python3 -m unittest test_finite_physical_prediction \
  test_finite_measurement_graph test_finite_covariance_rank3 \
  test_shipping_finite_identity test_core test_proof_plan test_bias_families -v
python3 tools/stability/ou3_alt_contraction/shipping_finite_identity.py \
  --output /tmp/ou3-alt-finite-identity.json
python3 tools/stability/ou3_alt_contraction/benchmark_finite_rank3.py \
  --output /tmp/ou3-alt-rank3.json
```

The rank-three benchmark reports runtime, peak memory and exact rational
entrywise equality for H18 and A21. It does not report source-uniform enclosure
quality, subdivision savings or a dissipativity margin. The independent
`finite-physical-identities` CI job retains these reports. The inherited broad
ALT suite is separate; its shared Mahony prerequisite can fail as recorded in
the research ledger. A green isolated job must not hide that failure.

## Decisive missing object and next work

**The requested source-uniform finite 600-step word is NOT materialized.**
Local finite identities plus observed same-history operands do not supply the
analytic frontend/tuner/covariance/guard graph. The present `physical_word.py`
Jacobian cocycle remains rejected by `assert_finite_storage_master`.
No high-precision finite-word feasibility or common-storage search is authorized
by the regression reports, and none is claimed here.

Bind the finite predictor and measurement descriptors to explicit actual runtime
successors over the analytic source family. Preserve physical angular defects,
all correlated moments and bias histories, every coefficient product, all
accepted/rejected/not-due branches, floors, conditioning, asynchronous events,
and the actual H18->A21 release guard. Expose every completed literal prefix;
selected observation points are not that complete prefix graph.

Only after that object closes, run its high-precision feasibility diagnostic
and first common joint24 storage search with bounded neutral/source supply.
Then proceed to rigorous source cover, endpoint dissipativity, every-prefix
retention, the physics-compatible Live basin, hybrid landing, aggregate fresh
entry, both startup paths, and finite precision. The two-strike architecture
review rule still applies; no common-storage attempt or strike is implied by
an incomplete graph or an inherited startup-validation failure.

## Gate status and failure classification

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Original `P4_PASS` / `P5_MAY_START` remain untouched by ALT.

The current gap is C/E (finite representation / missing source attachment), not
A (filter instability) or a proved infeasible common metric. There is no certified
rho, ultimate bound, retained basin or finite capture time yet. Do not substitute
more seeds, thinner boxes, a frozen observer, wordwise S reset, position reanchor,
or a smaller entry set for the missing finite source-uniform relation.
