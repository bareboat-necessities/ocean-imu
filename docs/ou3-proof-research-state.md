# OU-III proof research state

## Current hypothesis and controlling inequality

The finite-window source definition is insufficient for the requested indefinite
finite bound on all 18 motion errors. Before enclosing a full source cover, test
whether a finite uniformly coercive working tube can contain every admitted
continuation. The controlling necessary inequality is `W <= L < infinity`, not a
frozen-word spectral radius. Canonical P3 remains `delta=1e-18`; P4/P5 remain false.

## Executed result and exact failure analysis

`ou3_brmm_infinite_continuation.py` materializes the exact quiet-source ambiguity
column `B(h)=E_p+h E_S`, with one closed radial position parameter retained for
all time. Its polynomial event identities cover H18, A21 and the joint 24-state
coordinate. Prediction advances physical time; S Joseph residuals cancel against
the SAME physical S column; vector updates, finite common resets/projections,
release, reinitialization and identity branches preserve the ambiguity. No
independent gain, covariance, tuner coefficient or residual box is used.

The published primitive qualification freezes V_m and P_m but explicitly supplies
only `S_dot=p` and short-window Delta S, not an indefinite S bound. It admits
`p=+d` and `p=-d`, `v=a=omega=0`, fixed attitude and identical sensors. Zero bias is
separately admitted from each BIAS0/1/2 module. Fresh entry still has centered
`e_S,L=0`. Thereafter the same shipping estimate gives

`e_S,L^+ - e_S,L^- = 2 h d`,

so `max(||e_S,L^+||,||e_S,L^-||) >= h ||d||`. At least one fixed-sign history has
unbounded limsup. Uniform coercivity gives `max(W_+,W_-) >= m_- h^2 ||d||^2` in
fixed units, even with different source-dependent metrics and compatible metric
memory. For a 1/8 m witness, any 300 m*s WORKING radius fails for at least one of
the pair after 2400 s post-Live. This is not a fresh-entry radius counterexample.
The continuous physical-S energy on a 3 s word is
`||d||^2 (3 h^2+9 h+9)`, despite every Delta S being at most `3 ||d||`.

Classification **B for the finite-window-only indefinite target**: no finite
ultimate bound/working tube follows from that source definition. If the broader
phrase “bounded primitives or explicit forcing budgets” was intended to impose
an additional uniform S bound, its missing qualification is **E**; the witness
is not claimed to satisfy such an added premise. The source declaration's
Normal-Live admission flag remains false. This is not an established nominal
filter instability, a P3 counterexample, or a refutation of conditional ISS with
an unbounded physical-S input. The exact common-S-origin lemma remains valid.

## Native experiment failure and replan

The first paired-wrapper regression assumed H18->A21 release within 3 s after
Live. It failed only that assertion; states, covariance and source-error
identities stayed equal. The shipping default holds accelerometer bias until
magnetic-reference refinement completes its **30 s** window. The invalidated
hypothesis is the test's release deadline, not the symbolic event identity or
a runtime capture theorem. The revised regression follows that unchanged guard
for 33 s, without forcing release or modifying filter constants, and passes
6,600 post-Live IMU samples and 3,300 magnetic callbacks with an actual A21
release. All eight new exact/mutation tests and fourteen existing entry tests
pass; the theorem-input LaTeX smoke build also passes. These finite tests
regress the implementation binding, not the infinite-time conclusion.

## Validation compatibility

The source/runtime audit suites initially failed on unchanged main fixtures:
`unfrozen_physical_constants` no longer exists in the V4 declaration, and a
5 m/s^2 test signal no longer exceeds its 8 m/s^2 acceleration cap. These are
stale test assumptions, not proof or runtime failures. The source test now
checks V_m/P_m/A_m against the qualified primitive object, keeps the remaining
constants including S_m unresolved, and rejects mutations/promotion. The runtime
test places positive and negative controls on opposite sides of the actual
acceleration/body-rate caps. Both suites pass without changing source constants,
shipping code or audit thresholds.

## Retained facts and frozen dead ends

Fresh v/p/S/aw/bg/ba means are held zero; one common S-origin is removed only once.
BIAS0/1/2 keep separate driver lifts and bounded-bias projection proofs; BIAS2
separation is not mandatory. Qualified runtime Live/H18 timeout is 150 s, not a
P4 capture time. Conditional binary32 ISS remains separate from deployment.

The joint 24-state compatible-storage candidate retains useful point rates
H18=0.9996524356, A21=0.9959531012. The A21 marginal-motion covariance route
(rho about 10.119) and full-product nonlinear 3 s route (rho about 1.103 at
scale 5.5) remain rejected. Those failures do not justify weakening P3 or losing
motion/bias cross-information. Nor can longer words or tighter interval arithmetic
remove the exact position/S ambiguity lower bound.

## Alternatives and independent critic pass

The strongest objection to the negative result is that a *stronger* intended
source already bounds S or its forcing uniformly. The certificate explicitly
limits its negative conclusion to the published finite-window definition and
must not be presented as a counterexample to that stronger theorem.

Three qualitatively different resolutions are: qualify an actual same-history
uniform centered-S/forcing bound; change the performance objective to a specified
observable/regulated-reference output; or add physical position information to
the estimator. The latter two change the requested theorem or shipping filter
and are not made here. Shrinking the fresh-entry set or re-zeroing S every word
is not a resolution. A bounded S increment alone does not qualify the first option.

## Next falsifiable experiment

Establish the intended source's indefinite primitive/forcing premise from actual
physical assumptions, including Q, rather than supplying an arbitrary S cap.
At minimum every observation-equivalence class must have eventually bounded
centered-S diameter; a uniform working tube requires a uniform diameter bound.
Only after that scope is justified can full continuation coverage feed the joint
24-state compatible endpoint/prefix masters, first exit and H18 capture.
No complete source cover, numerical maximum P4 basin, endpoint/prefix margin,
post-Live capture time or end-to-end theorem is certified by the obstruction.
