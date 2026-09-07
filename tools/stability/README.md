# OU-III stability tooling

This package contains the retained OU-III stability theorem, interval arithmetic, SEA3 source, and P1-P4/P5 certificate tooling.

Production/validation utilities and unrelated OU-III engineering studies remain in `tools/`. Retired P2 history/source-node routes are not compatibility-shimmed here: canonical P3/P4 use the complete SEA3 moving-Riccati route.

## Current proof plan

The primary target is **bounded accelerometer-bias error plus regional practical
ISS of the other 18 errors**, in both held and active modes. The active filter
still executes all 21 states and every covariance/gain cross term. The target
is stated in `thm:sea3-bounded-bias-motion` in
`doc/kalman_ou_iii/w3d-sea3-stability-theorem.tex-part` and built by
`ou3_p4_bounded_bias_motion.py`. The earlier full-state contraction theorem is
a stronger unclosed extension, not a prerequisite for this weaker objective.
The physical bias package is stated and proved conditionally in
`doc/kalman_ou_iii/w3d-mems-bias-preconditions.tex-part` and recorded by
`ou3_mems_bias_contract.py`. Its open physical constants are not inferred from
the filter settings, simulation extrema, or the desired P4 margin.

| Step | Required use of the bias premises | Evidence needed |
| --- | --- | --- |
| BIAS0 / SEA0 | Qualify the assembled sensor, actual calibration, true residual root, temperature/strain mismatch and GM parameter range. Compose with the existing complete SEA3 source. | Stationary/thermal/restart/unit data, fit uncertainty and residual bounds; currently missing. |
| P1 entry and H18 hold | Retain the true physical bias while the estimate is held. Bound the held error and actual uncompensated offset as source inputs. | Reachable entry and continuing physical history; an estimate clamp is insufficient. |
| P2 / source cells | Carry one bias root, common parameters and driver history through the same window lineage. Split source coordinates and replay descendants. | BIAS1 dependence, same source ID and joint SEA3 membership; neither fresh event boxes nor ID attachment alone proves it. |
| P3 / detectability | Retain the full-state process matrices, finite tau, and actual R_S at delta=1e-18. | Rebuild the existing conditional H/A gate; separately establish coverage where the weaker target admits the projection boundary. Held H18 does not certify active-mode motion gains. |
| H18 to A21 | Transport the physical bias without reset; keep the rectangular error lift, actual held error, release covariance and all subsequent corrections. | Same-history reachable release, not an appended A21 fixture block. |
| Bias compactness | Projection preserves the estimate ball; BIAS0/1 bounds the same-history true bias. | Conditional error bound B_e=B_true+0.4; no bias convergence or hard Gaussian-OU cap is inferred. |
| P4-motion nonlinear graph | Keep the full bias prediction/correction/projection recurrence and actual active P/H/R/K. | The new domain includes the closed 0.4 ball; pre-projection auxiliaries may leave it. Existing 0.35-interior coverage cannot be reused silently. |
| P4-motion endpoint | Certify W_H,N <= rho_H W_H,0 + gamma_b D_b + gamma_s D_s + gamma_n D_n. | Strict motion decay factor, quantified coupling gains and exact source/error/forcing attachment. D_b bounds internal corrected bias error, not exogenous noise. |
| BIAS2 / gain sharpening | If used, certify the sector on the actual same-history graph. | Optional for the weaker objective; positive separation is not needed for compactness itself. |
| P4-motion every prefix | Certify finite Gamma and channel gains at every completed shipping event. | Endpoint-level invariance and strict chart retention, including all actual R_S events; prefix contraction is unnecessary. |
| P5-motion / capture | Carry bias/source/forcing through startup, H-to-A and all allowed hybrid events. | Start only after P4-motion, including a useful residual bound and retention, is certified; compactness alone supplies neither capture nor motion accuracy. |

The new master is C_N^T M_N C_N-rho_H C_0^T M_0 C_0 minus the same-graph
bias/source/noise energy forms. The full graph and its cross terms remain;
only the performance storage changes. Before enclosure work, evaluate this
inequality on an exactly attached complete word with candidate channel gains.
Its limiting performance quantity is the ultimate bound, not a raw VN/V0
ratio with nonzero forcing. No metric grid, longer-window search or source
reduction is authorized. The immutable old witness remains evidence about
the stronger storage only, subject to its documented attachment defects.

With B_e=B_true+0.4, the endpoint forcing budget is
C=gamma_b*T*B_e^2+gamma_s*Dbar_s+gamma_n*Dbar_n. The endpoint error floor is
C/(1-rho_H), and the every-prefix floor includes Gamma and the prefix input
budget. Model mismatch remains when sensor noise is off. The observed 2--3%
residual is not asserted as a uniform certified number; its output metric,
normalization, source range and quantitative gain still require proof.

CI rebuilds bias/source/P3 and emits separate `P4_MOTION_PASS` and
`P5_MOTION_MAY_START` flags alongside the unchanged stronger full-state flags.
The compactness and composition lemmas are conditional results, not numerical
motion-gain certificates. Both motion flags remain false until the source,
endpoint, prefix, retention and accuracy obligations close. BIAS0 device
qualification remains a separate deployment requirement; conditional work
does not wait for measurements. Literature does not supply numerical gains.

The retained-witness experiment is
`tests/kalman_ou_iii/ou3_p4_bias_witness_admissibility.py`, run by the existing
physical finite-map CI workflow. It pins the original H18/A21 payloads and
directions, checks the initial and every-event bias state, and evaluates BIAS2
with the full stacked nonlinear/corrected-error cross term and Xi=V0. It also
audits the required F/Q/H/covariance transports and exact finite reset energy
identity using `ou3_p4_retained_word_attachment.py`. The reset-deleted
diagnostic fails attachment as a gauge representative; its expansion is not
a canonical counterexample. Source membership remains undetermined because
the capture lacks its hard driver, phase, joint response and nominal history.
Hardware
qualification is not an execution gate for this conditional experiment.
The current result retains interior A21 expansion even with positive point
separation; see the research ledger. A positive sampled BIAS2 ratio is neither
a uniform sector nor proof of endpoint dissipation.

The 4.6 MiB fixture `tests/kalman_ou_iii/fixtures/ou3_p4_retained_payloads.tar.xz`
preserves the four original payloads from run `34102009588`, artifact
`10011024466`. The manifest pins their individual SHA256 values and the
archive SHA256. CI extracts these bytes and uses the frozen directions/scales;
it does not regenerate a simulation, select another window or depend on an
expiring Actions artifact. Run the audit with `--unpack-retained`,
`--payload-prefix /tmp/ou3_p4_physical_payload` and an `--output` JSON path.

Physical admission is checked separately by `ou3_physical_sea3_membership.py`
and the public-API C++ probe `ou3-sea-generator-membership.cpp`, against the
pinned v1.1.3 generator. `ou3_sea3_response_union.py` defines the authorized
SEA3+ union: the entire unchanged linear-vessel branch plus a correlated
third-order Stokes wave-following branch. The branch is fixed at the source
root; neither eventwise branch switching nor independent harmonic phases or
coefficients are allowed. All root directions/phases are included, not one
seed, and the finite Stokes model does not replace the legacy continuum.

The exact linear-branch squared-gain test at 2.4 Hz remains 2<=1 (rejected).
CI additionally checks all actual Stokes amplitudes, steepness and dependent
coefficients; model admission there does not assert complete word admission.
The linear moment/period lemmas stay linear-only; Stokes uses its full output
and same-history discrete frontend. The conditional P3 matrix implication is
rebuilt for both branches at delta=1e-18 under unchanged common Normal-Live
premises. The retained stronger P4 gate reports endpoint, finite prefix gain
and retention separately for both branches and both modes; missing coverage
prevents promotion of that stronger P4/P5 pair. The motion gates have their own
gain, projected-domain and retention obligations above.

The next motion proof-plan gate is common-root source/bias plus exact finite
error/forcing/reset/storage attachment, followed by the motion/channel-gain
master. Valid BIAS2 sectors may sharpen it but are not mandatory.
The simulator's turn-on bias/random walk differs from the audit's unforced
zero true-bias premise, and its finite-difference gyro needs a separate defect
bound. Neither is silently repaired by the response-family enlargement.

`tools/stability/ou3-source-endpoint.cpp` captures a new conditional
Stokes point from power-on through the two fixed observer-clock windows. One
shipping observer owns every nominal update, covariance reset, tuner commit
and actual R_S. Root phases/directions, raw warmup inputs, true zero residual
bias and actual sample/magnetic prefixes are retained. No initial state or
covariance is overwritten, and this is not a replacement for the immutable
retained witness. The Python endpoint audit independently reconstructs the
source and evaluates full-state endpoint storage at 80 decimal digits from
the float shipping states/covariances. It explicitly measures the physical
truth's defect relative to the committed F_LL prediction. A forced trajectory
ratio is not the homogeneous P4 ratio; no result of this capture alone can
authorize a sector search or falsify homogeneous contraction. Intermediate
S/accelerometer/reset graph capture and homogeneous attachment remain open.

The connected subevent test is
`tests/kalman_ou_iii/ou3_p4_connected_motion.py`. Its read-only trace overlay
must strip back to the shipping header and reproduce the uninstrumented
root/input/endpoint files byte-for-byte. It checks the full actual gains,
Joseph updates, finite resets/projection and source/error residuals at every
shipping event. In particular, physical S=0 has residual
`delta_S-S_true`: the `-S_true` source input is retained alongside the existing
actual-R_S stabilizing correction. The output includes full-information
18-error storage and eventwise signed energy changes, not a fitted gain or
an outward source-uniform certificate. Its stopping gate is point attachment;
the nonlinear graph sectors and augmented motion/gain master must follow only
on the attached graph. No earlier full-state ratio is promoted by this test.
