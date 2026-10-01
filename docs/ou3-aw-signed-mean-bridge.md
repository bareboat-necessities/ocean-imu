# Source-uniform signed AW-mean bridge

The physical side is no longer an independent open gauge assumption for the
candidate qualification. The existing bounded-velocity/jerk theorem plus
H_a=60 s,C_a=.05 m/s and the slow amplitude bound leave

    M_AW = .89553839501604595 m/s^2

inside Corollary A*'s g/5 field-separation threshold on every complete 16-s
window.

What remains is one full-state causal-reader inequality, not AW-only passivity.

Let alpha_k be the normalized nonnegative accelerometer-row weights of one G0
window and mu_hat=sum alpha_k a_hat_w,k. The exact literal AW loop gives

    mu_hat = mu_phys + mu_slow + mu_fast + L_AW b_full,        (AWB1)

where b_full is the chronological stacked source containing the SAME
accelerometer innovation, S corrections, OU prediction/source terms,
attitude/BA nuisance, reset/projection and arithmetic defects. L_AW is formed
by backward propagation through the complete 21-state literal factors.

The decomposition is defined by algebraic pairing: physical process and
accelerometer columns sharing a physical p/v/S/a primitive are added before a
norm; internal slab endpoints telescope; the S source is retained in the same
reader. No marginal AW Joseph inequality is used.

The candidate closes Corollary A* whenever

    sup_{same reachable 16-s history} |P_B L_AW b_full|
       < .89553839501604595 m/s^2.                            (AWB2)

This is a much larger and dimensionally correct target than the historical
.003056118-rad candidate bridge. If AWB2 holds, then

    |P_B mu_hat| < g/5,

so Corollary A* supplies a positive attitude Gram floor. MAGNETIC SERVICE plus
Lemma T and Theorem G0 then supply the injection-free aggregate six-column
floor. Literal injection transport and finite-error/prefix retention remain
their existing downstream obligations; AWB2 does not silently close them.

AWB2 must be evaluated with the complete paired two-Abel reader. Independent
bounds on process, acceleration, S, or gain extrema are not a proof because
they destroy the cancellation. A carried replay below the target is feasibility
evidence only; theorem closure requires a factor/interval enclosure over the
compact reachable 16-s coefficient/history class.
