# Bounded internal bias: the motion cocycle and its remaining P4 obligations

## Current result

The corrected sample-entrance motion operator is Schur on both attached
600-sample coefficient words. This establishes a **point, frozen-coefficient
quadratic-storage construction**, not uniform nonlinear P4. All P4 and P5
promotion flags remain false. The full shipping 21-state factors, corrected
bias, actual anisotropic R_S, finite reset/projection factors and physical
forcing are retained. No filter, P3 threshold, physical domain or quality gate
is changed.

## Exact reduction with feedback

Write the completed-event error as e=(x,b), with x in R^18 and accelerometer
bias error b in R^3. For each actual nonlinear history, use the retained exact
finite factorization e_next=A(e,zeta)e+B(e,zeta)u. Compose **all 21 rows within
one IMU sample**, including each actual S, accelerometer and asynchronous
magnetometer correction and its immediate reset/projection. At its completed
prefix ell this gives

    e_(i,ell) = F_(i,ell) e_i + f_(i,ell).

Both coefficients and forcing depend on the same physical/error history.
The forcing includes the latent truth defect, gyro discretization defect,
measurement mismatch and -S_true. Neither zero sensor noise nor a motion-only
performance storage removes these terms.

Let E_H and E_b select motion and bias, respectively, and define

    A_i = E_H F_i E_H^T,  G_i = E_H F_i E_b^T,  d_i = E_H f_i.

Then the following identity holds on the **actual** joint history:

    x_(i+1) = A_i x_i + G_i b_i + d_i,
    x_n = Phi_(n,0) x_0
          + sum_(i<n) Phi_(n,i+1) (G_i b_i + d_i),
    Phi_(n,j) = A_(n-1) ... A_j.

Proof: select the motion rows of the full sample identity and substitute
recursively. The same substitution using F_(i,ell) proves every-prefix
identities. No independence of b_i, no vanishing motion-to-bias feedback, and
no convergence of b_i is needed. The active gains still contain the full
covariance cross blocks. This is an identity for the executed estimator,
not a reduced estimator or a reset of its bias.

Crucially, Phi_(n,0) is generally **not** the motion principal block of the
full 21-state word product. That latter block contains feedback through bias
across sample boundaries. Conversely, multiplying motion blocks at every
subevent would incorrectly discard within-sample feedback and change the
sample-entrance bias-energy convention. Full-sample composition before
selection avoids both mistakes.

## Bias compactness is independent of the cascade condition

For a projection radius R, Pi_R(z)=z min(1,R/||z||), with Pi_R(0)=0, satisfies
||Pi_R(z)||<=R for every z. Holds, norm-preserving frame changes, predictions
with scalar factor in [0,1], and projection after injections therefore
preserve the estimate ball by induction. Under a qualified same-history
true-bias bound ||b_true||<=B_true,

    ||b_error|| <= B_true + R.

This is the existing conditional compactness lemma in the BRMM theorem.
It does not require true bias to lie in the estimate ball, and it does not
require zero motion-to-bias feedback in active mode. The unprojected injection
is an auxiliary coordinate and may leave the ball.

A linear map frozen at an interior trajectory uses that trajectory's
projection multiplier, often one. Extrapolating this multiplier to other
initial states is not the radial projection at those states. Thus the earlier
A21 frozen-map value 0.400126 does not falsify nonlinear estimate-ball
invariance or introduce a new exact-real compactness premise. It remains
useful evidence that frozen coefficient extrapolation is not a uniform
nonlinear enclosure. The captured binary32 radius is 0.4000000059604645;
shipping floating-point projection still needs its own rounding enclosure.
Neither a reference projection test nor this real-arithmetic proof supplies
that enclosure.

## Constructive point storage, without an eigenbasis assumption

For a fixed motion word T with spectral radius r<1, choose q=(1+r)/2 and a
positive physical scaling Q=diag(1/radii^2). The convergent Stein series

    M = sum_(k>=0) ((T/q)^k)^T Q (T/q)^k

is positive definite and satisfies

    M - (T/q)^T M (T/q) = Q,
    T^T M T = q^2 (M-Q) < q^2 M.

Convergence follows from Schur stability, including nontrivial Jordan blocks:
a polynomial times (r/q)^k is summable. The k=0 term gives positive
definiteness. Subtract the shifted series to obtain the equation. Therefore
this construction does not depend on diagonalizability or invertible
numerical eigenvectors. The producer solves the Stein equation numerically
and reports conditioning and residuals; it does not claim outward proof of
its floating-point candidate.

The statement `for every individual word there exists a contracting metric`
does not provide compatible metrics for consecutive words. For example,

    U = [[1/2, 2], [0, 1/2]],  V = [[1/2, 0], [2, 1/2]]

each has spectral radius 1/2, but VU has spectral radius greater than one.
This is a counterexample to that logical implication, **not** a claimed
admissible OU-III source or a counterexample to the OU-III theorem. A
source-dependent metric must contract from M(zeta_before) to
M(zeta_after), have uniform coercivity, and respect the actual continuation;
a separately selected metric for each captured word proves none of these.

## Separate bias and common-template supplies

For a fixed coefficient sequence let the motion lift at a prefix be

    x = T x_0 + B beta + f alpha,
    D_b = sum_i h_i ||b_i||^2 = ||beta||^2,

where beta stacks sqrt(h_i)b_i at **sample entrances only**. Keep one scalar
alpha on the entire captured physical forcing template; alpha is not an
independent input at each event. The Gramian Q_b=B B^T can be propagated
without storing all columns. At a sample prefix with motion blocks A,G,

    T_prefix = A T_previous,
    Q_b,prefix = A Q_b,previous A^T + G G^T/h_i,
    f_prefix = A f_previous + E_H f_(i,ell).

For an endpoint metric M and rho exceeding ||T||_M^2 but below one, put

    X = M^-1 - T (rho M)^-1 T^T > 0.

Then

    Q_b/gamma_b + f f^T/gamma_s <= X

implies x^T M x <= rho x_0^T M x_0 + gamma_b D_b + gamma_s alpha^2.
Proof: concatenate [T,B,f], apply the inverse block-diagonal input cost,
and compare its output Gramian with M^-1. Equivalently, the weighted map
has operator norm at most one. This bounds all coefficient input choices
and hence their actual correlated subset. A zero channel is omitted rather
than divided by zero. The producer chooses two separate gains using a
fixed 2.01 split of the whitened channel norms, not a gain grid or observed
trajectory-energy fit. The result is a sufficient outer test, not a claim
that independently varied bias slots describe physical histories.

For the retained nonlinear theorem, the same inequality must hold uniformly
for coefficients generated by its actual graph. Prefix supplies, admissible
source budgets, endpoint level invariance and strict physical/chart retention
are still required. BIAS2 may sharpen the actual joint-graph coupling, but
no unproved BIAS2 sector is used in this construction.

## Executed point diagnostic

The attached source is the existing PM 1.5 m 28 ft RAO capture, with its fixed
H18 and A21 clock windows. No new source, direction-dependent replay, or
window search is used. Maximizing coefficient directions are diagnostics,
not asserted physical states. The output records capture hashes, every
prefix, 80-digit evaluations of the worst endpoint/prefix directions, and
signed operation costs. Decimal arithmetic starts from binary64 factors,
not from a high-precision nonlinear observer.

| Quantity | H18 | A21 |
| --- | ---: | ---: |
| Spectral radius of sample motion cocycle | 0.997663960671 | 0.960675313400 |
| Constructive complete-word quadratic ratio | 0.995550332926 | 0.953993263418 |
| Worst completed-prefix ratio in that metric | 11.2178082213 | 1.63544556380 |
| Endpoint supply rho | 0.997775166463 | 0.976996631709 |
| Bias-energy gain | 4847.73017438 | 117.303616700 |
| Common-template gain | 0.237245019053 | 0.025801527123 |
| Whitened supply ratio (must be <1) | 0.756418187803 | 0.781065858240 |

The active sample cocycle differs from the principal block of the complete
21-state product by operator norm 0.516227538638. It reconstructs motion
when driven by the actual full-factor bias sequence to below 1e-15 in both
modes. The respective point metric condition numbers are approximately
3.33e9 and 5.11e9; small residuals and 80-digit directional agreement do
not turn those solves into uniform rigorous certificates.

## Failure analysis and critic decision

The endpoint coefficient supply is feasible, but compactness alone plus the
relaxed L2 bias budget does **not** certify useful coordinate retention.
With zero initial motion deviation, B_true=0, bias radius .4 and unit
captured-template amplitude, the sufficient every-prefix bounds reach
1.3212 chart radii, 1.1134 velocity radii and 3.5159 latent-acceleration
radii in H18. A21 reaches 3.3272 latent-acceleration radii. These are upper
bounds over a relaxed coefficient input class, not attained excursions of
the nonlinear source graph. Initial motion and repeated-word accumulation
have not even been added; no retention claim is justified.

Classification: loss-of-correlation / supply-enclosure failure, not a
counterexample to the theorem, projection invariance or point contraction.
The failed hypothesis is that compactness with an unconstrained entrance
energy budget alone supplies a useful retained motion tube. H18's actual
held bias, for example, is not a fresh independent vector each sample.

The strongest reason to abandon this bound as the final architecture is that
it discards exactly the bias temporal structure that limits physical forcing.
Refining scalar norms or increasing interval subdivisions does not address
that defect. Alternatives are: retain the true-bias/held-estimate recurrence
and active projection in a joint graph supply; use a source-dependent metric
with compatible consecutive-word contraction; or propagate correlated
source-centered coordinate tubes with a qualified word-entry set. The next
falsifiable experiment is a consecutive-word joint-graph test retaining that
bias recurrence and source continuation, with explicit coordinate budgets,
before any uniform enclosure work. No further independent-energy refinement
is justified by this diagnostic.

The remaining closure is uniform nonlinear/source-dependent contraction,
useful endpoint and every-prefix retention with admitted budgets, and the
separate canonical P3 execution-premise/source checks. Device qualification
is not a reason to stop conditional mathematical work. P4-motion and the
stronger P4 remain unclosed; neither P5 may start.
