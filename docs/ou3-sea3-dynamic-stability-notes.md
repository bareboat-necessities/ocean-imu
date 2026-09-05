# OU-III SEA3 P3 continuation after #490

This branch corrects the canonical P3 source architecture before continuing the numerical closure.

## Immutable source semantics

P3 proves the complete Normal-Live H18/A21 word for the compact, phase-continuous SEA3 sea-state source

\[
\zeta_k=(x^s_k,\lambda_k,z^t_k,q_k),
\]

where \(x^s\) is the SEA3 oscillator/shaping state, \(\lambda\) is the compact three-partition directional-sea parameter state, \(z^t\) is the exact measurement-only front-end/tuner state, and \(q\) carries scheduler/hybrid source memory.

There is no alternate promotable source. In particular P3 may not be generated or pruned by replay, a Gaussian finite-horizon event, spectral moments alone, arbitrary bounded inputs, a tuner rectangle, independent \(\tau/\sigma/R_S/T_S\) extrema, the retired P2 predecessor graph, or a selected four-S word.

Stochastic concentration remains only a later forcing/corollary calculation. Configured measurement covariance remains in every Riccati update.

## Full shipping word

The same SEA3 realization must generate every committed \(\tau_k,\sigma_{aw,k},R_{S,k},T_{S,k}\), every \(F_k,Q_k\), every due S=0 correction using the actual deployed SpectralMSE per-axis \(R_S\), every valid accelerometer Joseph correction, asynchronous PE events, covariance-floor events, and both H18/A21 recursions.

The useful P3 gate remains

\[
\Omega_W-10^{-18}P_W\succeq0
\]

by validated full-matrix LDLT in both modes. P4 remains blocked until this is actually closed.

## Current execution status

The #490 scaffold never executed the source-reachable full word. This branch now:

- restores SEA3 compactness as a theorem-domain property;
- removes the Gaussian good-event source shortcut from canonical P3;
- converts shipping `0.2f`/`0.3f` measurement literals to binary32 before outward variance enclosure;
- adds validated shipping prediction-matrix primitives that consume only SEA3-derived per-sample coordinates.

Still open: interval-materialize the complete phase-continuous SEA3 family through its compact transition relation and the exact front-end/tuner/scheduler state, then execute all samples through the literal H18/A21 word and close both endpoint LDLTs. No diagnostic subset may set the canonical pass flag.
