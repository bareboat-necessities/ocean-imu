# OU-III stability tooling

This package contains the retained OU-III stability theorem, interval arithmetic, SEA3 source, and P1-P4/P5 certificate tooling.

Production/validation utilities and unrelated OU-III engineering studies remain in `tools/`. Retired P2 history/source-node routes are not compatibility-shimmed here: canonical P3/P4 use the complete SEA3 moving-Riccati route.

## Current proof plan

The controlling target is strict endpoint dissipation and finite every-prefix
gain plus chart/domain retention on one complete correlated SEA3 word, for both
H18 and A21. The physical bias package is stated and proved conditionally in
`doc/kalman_ou_iii/w3d-mems-bias-preconditions.tex-part` and recorded by
`ou3_mems_bias_contract.py`. Its open physical constants are not inferred from
the filter settings, simulation extrema, or the desired P4 margin.

| Step | Required use of the bias premises | Evidence needed |
| --- | --- | --- |
| BIAS0 / SEA0 | Qualify the assembled sensor, actual calibration, true residual root, temperature/strain mismatch and GM parameter range. Compose with the existing complete SEA3 source. | Stationary/thermal/restart/unit data, fit uncertainty and residual bounds; currently missing. |
| P1 entry and H18 hold | Retain the true physical bias while the estimate is held. Bound the held error and actual uncompensated offset as source inputs. | Reachable entry and continuing physical history; an estimate clamp is insufficient. |
| P2 / source cells | Carry one bias root, common parameters and driver history through the same window lineage. Split source coordinates and replay descendants. | BIAS1 dependence, same source ID and joint SEA3 membership; neither fresh event boxes nor ID attachment alone proves it. |
| P3 / detectability | Consume the matched finite-tau prediction premise; retain all full-state process matrices and actual R_S. | Rebuild the existing conditional H/A gate at delta=1e-18. Device qualification remains a separate deployment blocker. If actual tau differs, retain mismatch forcing or explicitly re-certify the changed homogeneous family. |
| H18 to A21 | Transport the physical bias without reset; keep the rectangular error lift, actual held error, release covariance and all subsequent corrections. | Same-history reachable release, not an appended A21 fixture block. |
| P4 nonlinear graph | Use BIAS1 true-bias evolution for projection, and the full prediction/correction/projection recurrence for bias error. | Root-derived projection inputs; retain deterministic Kalman corrections inside the homogeneous graph. |
| BIAS2 / P4 endpoint | Certify positive separation on the actual corrected-error graph before using Pi_sep in the dense augmented master. | Uniform mu_sep and source identity are currently missing. A SEA3 spectrum and an OU PSD do not prove separation. |
| P4 every prefix | Use restricted graph sectors and finite Gamma_k, then prove chart, bias and source-domain retention. | Strict outward augmented LDLT for finite gain; prefix contraction is unnecessary. |
| P5 / practical ISS | Carry bias source/forcing through capture and all hybrid events, then consume certified P4. | Finite capture and disturbance propagation; stationary spectral variance is only a stochastic corollary. |

Before a new sector/enclosure search, replay the retained expanding A21 witness
through the exact admitted source/error map with BIAS1 attached. True bias and
bias estimation error are different: the latter receives every Kalman update.
If that witness remains admitted with rho>1, the same storage is falsified;
valid graph sectors cannot repair it. Conversely, failure of an S-procedure
relaxation alone does not falsify the storage. The qualitative alternatives are
dense source-structured storage, path-dependent storage, or a direct theorem
falsification. No metric grid, longer-window search or source reduction is
authorized by the new bias model alone.

CI rebuilds the bias contract, source, frozen P3 and canonical P4 chain and
publishes their separate closure flags. BIAS0 qualification, BIAS2 separation,
the complete source cover, P4 and P5 remain open until their actual evidence is
available. Literature establishes the model rationale, not those certificates.
