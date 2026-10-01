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

## Carried falsification result

The committed six-profile native AW audit contains the adversarial
sync-locked-rectification profile. Its largest literal 16-s signed AW-loop
error is

    0.371143 m/s^2.

Against AWB2's 0.89553839501604595 m/s^2 allowance this leaves

    0.52439539501604595 m/s^2

of carried margin, i.e. the allowance is about 2.41 times the observed worst
stress value. Therefore the literal mechanism passes the required first
falsifiable calculation. This number is not promoted to a source-uniform
constant.

## Source-uniform enclosure target

Do not intervalize K, tau, sigma_aw, R_S, T_S or covariance independently.
That relaxation is already known to destroy the same-history cancellation.
For a complete 16-s word h define the composed scalar/vector functional

    Phi(h) = P_B L_AW(h) b_full(h).

The theorem target is the single enclosure

    sup_{h in H_16} |Phi(h)| < 0.89553839501604595,

where H_16 is the compact reachable class carrying the one SLOW+FAST history,
literal scheduler, tuner state, covariance, reference and operation branches.
A valid interval/factor proof must propagate common dependency symbols through
the chronological factors and pair physical process/acc/S columns before each
wrapping step. The current generic independent interval Riccati boxes are not
a valid substitute.

The carried margin 0.52439539501604595 m/s^2 is the available wrapping,
nonlinear and float32 headroom for this enclosure.

## Uniformity theorem and remaining quantitative diameter

On any fixed regular event stratum, the shipping prediction, Joseph correction,
S scheduler, AW synchronization, reset and projection maps are continuous in
the carried state/covariance/tuner/reference variables as long as their already
declared positive innovation-noise floors and solve guards hold. The candidate
physical history class is equibounded/equicontinuous on 16 s: p,v,a and jerk
are bounded; SLOW bias is bounded/Lipschitz; FAST held histories are bounded
and have the all-placed-window primitive cap; attitude lies in compact SO(3);
BG/BA estimates are projected; tuner variables are clamped and coupled.
Held-H18 LIN BIBO and the captured release result provide the needed finite
LIN root set on the qualified branch. Closed scheduler/gate strata and the
finite number of operations on 16 s therefore make the reachable regular
16-s word class compact. Phi is continuous on each stratum and has one-sided
continuous limits at hard-event boundaries retained by the literal branch
semantics. Hence

    B_AW := sup_{h in H_16} |Phi(h)|

exists and is finite.

This proves source-uniform FINITENESS, but not the strict numerical inequality
B_AW < .89553839501604595. For strict closure choose the committed carried
stress family as centers h_j and a dependency-preserving factor metric d_F on
the COMPOSED chronological factors. A sufficient finite-cover certificate is

    max_j |Phi(h_j)| + L_Phi * delta_F + E_nl + E_f32
      < .89553839501604595,                                (AWB3)

where every reachable h lies within delta_F of some center in the same event
stratum and L_Phi is a rigorous Lipschitz bound of the composed Phi map on
that stratum. Since max_j|Phi(h_j)|=.371143, the available combined radius is

    L_Phi delta_F + E_nl + E_f32 < .52439539501604595.

Neither L_Phi nor delta_F is presently certified in the repository. This is
the sole quantitative source-uniform gap in AWB2; independent gain/tuner boxes
are not substituted.

## Why a six-center numerical cover is not yet a proof

The candidate MARINE/SLOW+FAST assumptions exclude an all-time exact
field-axis collinearity execution: persistent nominal field alignment would
require the gravity-sized transverse physical acceleration identified by the
sampled force/field theorem, and integrating that acceleration contradicts
the all-time |v|<=5.5 bound. Together with the compact same-history class this
gives a positive qualitative distance from persistent exact collinearity.

It does NOT imply that every reachable 16-s composed factor word lies within
a known radius of one of the six committed stress words. The Riccati analysis
already shows no proved positive distance from the transverse AW-gain
cancellation manifold at a single finite word: S/magnetic corrections,
prediction and the changing accelerometer Jacobian can move the AW gain
numerator additively, and no current invariant bounds those linked changes by
less than the incoming singular distance.

Therefore assigning a numerical delta_F from the six carried words would be
an invented qualification. The finite-cover inequality AWB3 is a valid
sufficient condition, but its delta_F is presently unproved. The correct
source-uniform consequence currently available is qualitative: exact
persistent collinearity is excluded, so a finite superword strictness margin
exists by compactness. A numerical .895538395 AW-reader ceiling needs either
(a) a rigorous reachable-factor cover, or (b) a direct contradiction proving
that |Phi|>=.895538395 itself forces the forbidden persistent acceleration
pattern. Neither implication is currently established.
