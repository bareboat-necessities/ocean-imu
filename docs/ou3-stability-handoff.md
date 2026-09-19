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
6. Retain A21 inside the inactive accelerometer-bias projection region
   (bias-error radius <0.1748334 m/s^2), then bound the remaining nonlinear
   reset/tuner/arithmetic remainder. Close finite-error contraction with
   (sqrt(rho0)+eta)^2 < 1.
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
