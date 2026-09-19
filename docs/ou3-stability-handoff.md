# OU-III stability handoff

## Architecture

There is one theorem under simultaneous MARINE MOTION, IMU BIAS, and MAGNETIC
SERVICE assumptions on one persistent physical execution. Shipping code is
authoritative.

The proof route is
`construction -> capture -> finite H18 bridge -> finite reference refinement
and bias release -> recurring magnetically informed A21 -> regional practical
stability`.

## Key correction in PR #558

H18 is no longer asked to provide asymptotic contraction. PR #557 showed why a
particular full-state held-bias contraction formulation cannot work. That
negative result is retained in prose, while its experimental export and
finite-difference scripts have been removed.

The shipping wrapper itself supplies the reason to move on: accelerometer-bias
learning is externally held while the magnetic reference is provisional, and
the hold is released when refinement completes. The internal gate additionally
requires its accepted-magnetometer-update threshold and a one-second guard.
Under recurring MAGNETIC SERVICE the count/guard part is finite once reference
refinement is finite.

## Productive proof obligations

1. Qualify the assembled sensor/bias limits and all-time marine-motion
   membership.
2. Prove finite history-dependent startup/capture.
3. Prove finite magnetic-reference refinement and hence finite H18 release.
4. Bound the finite H18 bridge and prove the literal release retains the A21
   tail domain.
5. Establish a source-uniform complete A21 information floor in normalized
   covariance coordinates. MAGNETIC SERVICE supplies the heading/axial
   gyro-bias component; gravity/accelerometer and integral pseudo-updates must
   close the remaining directions. A full floor mu gives the linear comparison
   rho0<=1/(1+mu).
6. Carry the literal dissipative accelerometer-bias projection sector; do not require projection-inactive release. Bound the remaining nonlinear reset/tuner remainder and close finite-error contraction with (sqrt(rho0)+eta)^2 < 1. Arithmetic remains a separate additive supply.
7. Close every-prefix retention, recurring service, and finite-precision
   arithmetic.

## Reproduction

`cd tests/validation && python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias test_ou3_magnetic_service test_ou3_marine_motion test_ou3_no_mag_obstruction test_ou3_same_execution test_ou3_theorem_status test_ou3_tail_stability`

`python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json`

`make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test && ./tests/kalman_ou_iii/shipping_contract-test && ./tests/kalman_ou_iii/shipping_transition-test`

No H18 numerical diagnostic is a proof gate. The next implementation work
should target finite reference refinement/release and the analytic A21
dissipativity inequality.


## Latest analytic closure

The A21 operating point is now handled as genuinely time varying. The
tau-scaled, progress-preserving S scheduler plus compact dt/tau bounds yields a
uniform positive translation observability minor through an extended-Chebyshev
argument. Gravity supplies a uniform 1.00665 m/s^2 tilt sensitivity floor;
MAGNETIC SERVICE supplies the missing heading/axial-gyro-bias information.
After quotienting the strictly stable active accelerometer-bias OU mode, the
neutral information Gramian is pointwise positive on strict compact hybrid
cells, hence has an existential uniform floor mu_N>0.

Positive process-noise densities give uniform complete controllability, and the
shipping covariance sync/release operations preserve compactness. Therefore
the linear A21 LTV Kalman error word has a source-uniform rho0<1. On the strict
6-degree / 0.15 m/s^2 inner domain, projection is inactive and the same-history
nonlinear remainder is smooth with eta(r)->0, so a positive local radius exists
with sqrt(rho0)+eta<1. Arithmetic is additive supply.

Next work is to make these existential margins constructive: enclose numerical
mu_N/rho0, derive an explicit nonlinear radius, bound arithmetic supply, and
prove capture/release retention into that radius.


## Validation note

The proof branch was rebuilt once on current main to remove stale generated
evidence from the PR. Subsequent validation runs should compare ordinary
descendant commits; the one synchronization run immediately following that
history rewrite can fail its shallow before/after classifier because the old
pre-rewrite object is intentionally no longer in branch history.


## Marine-forced attitude diversity

The attitude information floor follows from the existing physical contracts.
For `f=a-g` and true field `B`,

`integral_0^T f x B dt = Delta v x B - T g x B`.

Thus
`T^-1 ||integral f x B dt|| >= g B_h,min - 2 V_max B_max/T`.
With the declared limits this becomes positive after 5.60844 s and is
43.97475 (m/s^2) uT for T=8 s. Every 8 s marine history therefore contains
non-collinear gravity/specific-force and magnetic sensitivity. Recurring
MAGNETIC SERVICE transports an applied magnetic sensitivity into that interval.
Two 8 s subwindows form the 16 s A21 proof word and expose gyro bias through
attitude propagation. No additional attitude-excitation assumption is used.


## Quantitative enclosure status

A first fully analytic normalized translation enclosure is now constructive.
On the 16 s word, selecting S updates after times 0, 8 and 16 s with the
literal maximum scheduler delay 0.156 s, state scales
(V,p,S)=(5.5,8.1,1100), and worst declared S-noise standard deviation 100,
the determinant/Frobenius certificate gives

`mu_trans >= 2.04734e-3`.

This is a real lower bound, not a sampled singular value.

A deliberately sparse two-epoch attitude/gyro calculation gives a much weaker
candidate scale (~6.18e-7) and therefore is **not promoted as the final
shipping mu_N certificate**. The reason is important: MAGNETIC SERVICE is an
already-transported two-coordinate heading/axial-bias Gramian over a window,
not an instantaneous pure-heading row. A tight full certificate must compose
that actual 2-D service Gramian directly with the many accelerometer rows in
the same 16 s word; replacing it by a fictitious instantaneous attitude
measurement would be an invalid shortcut. The next quantitative proof step is
therefore the aggregate Schur/Gramian bound on the literal partition, using all
recurring accelerometer and S information rather than two sparse rows.

The float32 arithmetic path is likewise separated correctly. A straight-line
kernel with n rounded operations has the standard gamma_n bound
`gamma_n=n*u/(1-n*u)`, u=2^-24. This is implemented as a certificate
primitive, but no whole-word arithmetic supply is claimed until literal kernel
operation counts and magnitude envelopes are composed. Roundoff remains
additive supply and is not allowed to consume the nonlinear derivative margin.


## Constructive long-word neutral information

The weak gyro-bias scale is handled by accumulation, not by strengthening the
physical assumptions. The proof word is now 2048 s. At every disjoint one-second
MAGNETIC SERVICE root, the service inequality J_hb>=I permits extraction of one
unit of **pure root-heading information** (subtract diag(1,0); the remainder is
PSD). The simultaneous root accelerometer sample has no preceding gyro-bias
transport. With the commissioned detector-band residual bound 0.3 m/s^2,
shipping vibration gain 0.75, and nominal accelerometer std 0.2 m/s^2, the
literal effective accelerometer std is bounded by

`sigma_acc,eff <= hypot(0.2,0.75*0.3)=0.3010399 m/s^2`.

Together with `|f_z|>=g-A_max=1.00665`, `|f_h|<=8.8`, and one unit of
heading service, each service root supplies a pure attitude information floor
greater than 0.0129 in the declared proof coordinates.

Let `y_k=theta_0+C_k b_g` be those extracted root-attitude observations.
Bounded angular rate gives

`sigma_min(C_(k+1)-C_k) >=
  s_bg * 2 sin(Omega_max/2)/Omega_max`

with `s_bg=0.02`, so the normalized increment is at least about 0.01969.
The path-graph inequality
`sum |y_(k+1)-y_k|^2 <= 4 sum |y_k|^2`, combined with the root observation
`y_0=theta_0`, gives the joint attitude/gyro-bias floor

`mu_ag = j*c/(1+c),  c=(N-1) beta^2/4`.

For N=2048 this exceeds 2.1e-3. The independent 16-s S certificate remains
`mu_trans>=2.04734e-3`; therefore the complete **fixed-coordinate neutral
Gramian** satisfies

`mu_N >= 2.04e-3`.

This construction never integrates attitude over 2048 s: it uses only
one-second transport increments, so bounded rotations cannot cancel the
gyro-bias information.

### Important normalization correction

The number above is not yet a covariance-metric contraction factor. The
information-form identity `rho0<=1/(1+mu)` requires the Gramian to be whitened
by the actual root covariance. If, in the same proof coordinates,
`P_root>=p_min I`, then

`mu_cov >= p_min mu_N`,
`rho0 <= 1/(1+p_min mu_N)`.

Using the fixed-coordinate 2.04e-3 directly in the latter formula would be a
coordinate-dependent and invalid shortcut. The next quantitative certificate
is therefore a recurring lower enclosure `p_min>0` for the shipping A21 root
covariance. Only after that enclosure is established will rho0, the explicit
L2/r* radius, and the arithmetic/practical-radius bounds be promoted.


## Verified interval Riccati certificate implementation

The covariance-normalization blocker now has executable fail-closed certificate
machinery. A candidate covariance trajectory may propose a spectral box
`p_min I <= P <= p_max I`, but it cannot promote it.

For every accelerometer, integral, and magnetic innovation interval, write
`S=S_0+E`. The certificate requires
`lambda_min(S_0)>=a`, `||E||_2<=r<a`. Then the midpoint inverse residual is
bounded by `r/a<1`, and Neumann/Krawczyk inclusion gives
`||S^{-1}||<=1/(a-r)`. If this strict inequality fails, the certificate fails
closed.

The complete interval Riccati image must then satisfy, with outward-rounding
slack included,

`[R(P,parameters)] subseteq [p_min I,p_max I]`

for every covariance in the proposed box and every admitted shipping parameter
box. The certificate also requires explicit inclusion of the periodic
`a_w` covariance synchronization, H18-to-A21 bias-release covariance floor,
the complete bounded tuner/scheduler envelope, and float32 rounding. Only after
all flags and inclusions verify does the machinery compute

`mu_cov=p_min*mu_N`, `rho0=1/(1+mu_cov)`.

A committed status artifact currently remains OPEN. This is deliberate: the
residual/inclusion checker is now present, but the full outward-rounded shipping
Riccati image has not yet been generated. No midpoint-only p_min is accepted.


## Current covariance-normalization result

The literal entrywise midpoint-radius 21-state Riccati enclosure was executed after its one permitted near-identity Rodrigues refinement. It failed as an interval-conditioning mechanism: the predicted covariance spectral enclosure was [-18.7907040, 2502.41442]. The S innovation still verified with r/a=0.666663864<1, but accelerometer and magnetometer innovation boxes did not admit strict inverse certificates. Pointwise shipping covariances remain PSD; the negative lower endpoint is interval dependency, not a physical covariance counterexample.

Per the research protocol, do not subdivide this entrywise mechanism again. The next covariance-normalization construction must preserve PSD structure (spectral/square-root factor enclosure) and must itself be run before p_min is promoted. Until then constructive_root_covariance_floor, mu_cov/rho0, explicit r_*, arithmetic practical radius, capture/release-to-tail, and prefix retention remain open.


## PSD-preserving covariance review result

A one-sample source-range feasibility check of the replacement PSD-preserving representation is positive: sigma_min(F) >= 0.6440871, with covariance floor 4.14848e-7 after prediction and 1.49403e-7 after the maximal S/accelerometer/magnetometer correction sequence. These are non-promoting feasibility margins.

The remaining recurring-covariance task is now the exact process-noise factor: certify a full-rank lower factor/floor for shipping Q_AA, integrated-OU Q_LL, and active-bias Q_BA, then iterate the PSD representation over the actual cadence. A one-step scalar Q eigenvalue is expected to be extremely conservative for the triple-integrated chain, so it must not be used to manufacture a useless rho0; the factor/block structure must retain multi-step controllability.
