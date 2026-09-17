# OU-III ALT corrected observability formulation

## Selected architecture

The old requirement of strict coercive joint24 contraction on every admitted 3 s word is retired. PR #533 proves that an admitted quiet ungauged word contains the exact unipotent pair `(theta_z,bg_z)` with spectral radius one. This falsifies that theorem formulation, not filter stability.

ALT now uses two regimes. During ungauged intervals the proof works modulo the unobservable heading/axial-bias centre dynamics and must bound finite-prefix growth without claiming absolute-heading contraction. A yaw-only quotient is insufficient because `bg_z` remains neutral. The quotient is legal only when the actual map descends to it: `Q_after^T A N_before = 0`; deleting coordinates is not a proof of equivariance.

After magnetic gauge is established, the full-state theorem is conditional on an informative magnetic-service class. Attempted calls at gaps <=40 ms are retained runtime facts but are not themselves accepted informative observations. The theorem must establish a bound on gaps between accepted informative events and a positive transported heading/bias information floor over each superword.

## Controlling measurement

Before interval enclosure, measure on the actual same-history carried superword

`rho_W = lambda_max(Z^T A_W^T M_after A_W Z, Z^T M_before Z)`.

`Z` injects the canonical unsupplied coordinates; all output rows of `A_W Z` remain, including A21 bias corrections. A common metric must be the same `M` on every word. Compatible endpoint metrics must come from one declared storage law with uniform coercivity; choosing a separate favorable metric for each word is invalid.

The gauged motion-block spectra near 0.9964 reported by PR #533 are feasibility evidence only. They are not an upper bound on this storage ratio and their roughly 3.6e-3 distance from one is not yet an interval-enclosure budget.

## Magnetic service

For accepted gauged magnetic events in a superword, transport the actual whitened residual sensitivity to the same initial heading/axial-bias coordinates and form `G_B = sum B_i^T B_i`. The selected service class eventually needs explicit constants `T_B` and `alpha_B>0` such that informative-event gaps are <=`T_B` and `G_B >= alpha_B I_2` on the declared window. Those constants are deliberately unset until derived from runtime/source semantics. Rejected calls, ungauged calls, and heading-degenerate rows provide no such service.

## Non-promotion

This selection does not certify source-uniform magnetic acceptance, quotient equivariance, a common/compatible storage, corrected rho, interval enclosure, startup, Live, or end-to-end stability. `storage_search_allowed`, `ALT_STARTUP_PASS`, `ALT_LIVE_PASS`, and `ALT_END_TO_END_PASS` remain false. The eleven finite-master qualifications remain open but subordinate until the corrected rho measurement is feasible.

The independent P2/P3/P4/P5 proof track is unchanged, including frozen P3 `1e-18`.
