# Declared-domain coordinate retention

`domain-retention.json` records this experiment on the same passively attached
physical capture the three storage routes use: H18 at source time
32.995–35.995 s and A21 at 1142.995–1145.995 s, 600 predictions and 600
accepted accelerometer corrections each, with the full 21-state finite map,
actual covariance and anisotropic R_S, finite quaternion resets and bias
projection. Coefficients are frozen at that one history. No source search, no
coefficient fit, no filter change, no reduced domain.

Route 1 of `three-storage-routes.md` bounds a scalar information storage and
then converts that scalar back into an attitude excursion through the worst
direction of the metric. This experiment removes both conversions: the same
reachable set is propagated in the declared physical coordinates and each
coordinate group is compared with its own declared bound from
`tools/stability/ou3_proof_operating_domain.json`.

Every prefix reports two numbers. The certified upper bound is subadditive
over the independent declared balls and never understates the reachable
excursion. The attained lower bound is a maximizing unit functional and is
always achieved by an admissible initial state. A group is retained only when
its upper bound stays inside its own ball, and is definitely violated only
when its attained bound leaves it by more than the 1e-9 relative margin the
enclosure check itself allows.

The two bounds are maximized over the prefixes independently, because the
prefix whose subadditive sum is loosest need not be the prefix actually
reached furthest. Each is reported with the prefix index that attains it.

## The bias ball is not the limiter

With the closed .4 m/s² bias-estimate ball and the same common forcing
template as the only initial deviation, every motion group is retained with a
large margin.

| Reached / declared bound | H18 | A21 |
|---|---:|---:|
| 30° attitude chart | .0703 | .0041 |
| Gyro-bias ball | .0101 | .0045 |
| Velocity ball | .1713 | .2382 |
| Position ball | .0444 | .0844 |
| Integral-displacement ball | .0020 | .0046 |
| Latent-acceleration ball | .1431 | .1397 |
| Accelerometer-bias ball | 1.0000 | 1.0004 |

The rows are certified upper bounds. On the four that carry the conclusion —
velocity, position, the integral state and the chart — the attained bound is
within 1.3% of them, so those are tight numbers rather than a loose enclosure.
The gyro-bias row is the widest, 19.3% in H18, and both of its bounds are two
orders of magnitude inside the declared ball either way.

Route 1 reports the same initial set as a sufficient storage 6.4111 (H18) and
2.2693 (A21) times the 30-degree chart level. Storage is quadratic, so the
comparable linear excursions are 2.532 and 1.506 chart radii, against .0703
and .0041 here — factors of 36.0 and 364. The 6.41/2.27 failure is therefore
manufactured by charging a bias-driven velocity and displacement excursion to
attitude through the worst direction of the information metric. It is a
property of that scalarization, not of the reachable set.

This retires the separated bias budget as the source of route 1's chart
failure. It does not prove P4, admit a physical source, or establish anything
uniform over the nonlinear source family.

## The closed bias ball is not itself invariant in A21

The accelerometer-bias row above is the one that is not retained. H18 maps the
closed .4 m/s² ball into itself exactly. A21 does not: its attained bound is
.400126 m/s², exceeding the declared radius by 1.259e-4 m/s², a relative
3.146e-4. That is roughly 2e12 ULPs and is not roundoff, so A21's complete
word grows the worst bias direction rather than contracting it, and the ball
is invariant for one mode only.

The growth is small but it is a real open obligation: a bounded-bias statement
over repeated words has to carry it rather than assume the ball is invariant.
It is separate from the chart question above — every motion row is reached to
at most .2382 of its declared bound, and none of them depends on the bias
ball's own behaviour — so it qualifies the bias premise without reinstating
the separated budget as the chart limiter.

## The declared product box is not invariant

With the complete declared initial set the same word leaves the domain by a
wide margin, and the limiting group is velocity in both modes.

| Initial set | H18 worst | A21 worst |
|---|---:|---:|
| Full declared box | 35.358 (velocity) | 9.357 (velocity) |
| Without the 300 m·s integral ball | 4.033 (velocity) | 4.921 (velocity) |
| Bias ball and template only | 1.0000 (bias) | 1.0004 (bias) |

Every bound is positively homogeneous in the declared radii and the template
amplitude, so the largest box of the declared shape that this word retains is
2.83% (H18) and 10.69% (A21) of the declared one.

The single dominant source is the declared 300 m·s integral-displacement ball:
on its own it drives velocity to 159.7 m/s against a 5 m/s bound and position
to 379.3 m against 20 m in H18. That ball is not independently reachable. The
integral state is the running integral of the position state, so a 300 m·s
integral error together with a position error of at most 20 m requires a
sustained 20 m position error for 15 s. The declared bounds are a product of
independent balls and therefore contain states the filter's own kinematics
cannot produce.

Removing that ball still leaves velocity at 4.03 (H18) and 4.92 (A21), with
position and attitude as the next dominant sources, so a correlated
integral/position fact is necessary but not sufficient on its own.

## The correlated initial set is very nearly invariant

The product box is not the only correlated alternative, and the filter already
carries one: the covariance of the word's initial point, which holds exactly
the position/integral/attitude cross terms the box discards. It is a single
convex set rather than a product, so its image needs no subadditive step and
the excursion `sqrt(c)*||Pi_G T_j L||_2 + |alpha|*|Pi_G r_j|`, with
`P = L L^T`, is exact.

Reading the result in initial standard deviations makes the two sets
comparable. For each group the critical level is the largest `sqrt(c)` this
word keeps inside the declared radius, and the declared radius itself sits at
some number of initial sigma. Their ratio is what the complete word actually
costs.

| Group | H18 critical / radius, in sigma | growth | A21 critical / radius, in sigma | growth |
|---|---:|---:|---:|---:|
| Attitude | 9.679 / 9.679 | .001% | 113.48 / 113.49 | .005% |
| Gyro bias | 9.956 / 10.000 | .443% | 123.75 / 123.75 | .000% |
| Velocity | 4.974 / 4.979 | .117% | 174.92 / 176.01 | .622% |
| Position | 19.37 / 19.57 | 1.021% | 231.30 / 233.44 | .922% |
| Integral displacement | 550.3 / 567.4 | 3.111% | 1593.9 / 1607.7 | .869% |
| Latent acceleration | 11.93 / 12.79 | 7.200% | 112.87 / 112.88 | .008% |
| Accelerometer bias | 100.0 / 100.0 | .000% | 27.45 / 27.46 | .049% |

Every critical level is within 7.2% of the level at which the initial set
already touches its own declared radius, and within 1.1% on five of the seven
H18 rows. The complete word therefore barely expands this set at all: the
binding quantity is how many initial standard deviations the declared radius
is, not the word's dynamics. Against the same word the declared product box
overshoots velocity by 3436%.

The overall critical levels are 4.974 sigma (H18, limited by velocity) and
27.452 sigma (A21, limited by the bias ball). Both are what a covariance
consistency requirement looks like once it is stated quantitatively.

This exchanges one unproved premise for another and does not close P4. The box
was never a reachable set; this ellipsoid is the covariance the filter
believes, not a qualified bound on its actual error, and the runtime audit
records that actual covariances differ from the frozen P3 premises. The level
above is the requirement, not the certificate.

## The forcing template is not the obstruction

The critic's second architecture asks whether a particular response restarted
every word is what fails. With every initial ball at zero, the same common
template reaches at most .0470 of a declared bound in H18 and .00472 in A21,
both in velocity. It consumes under 5% of the smallest budget in play, so a
globally bounded particular solution would not change any conclusion here.
That architecture is answered at this capture rather than left open.

## Every route executed on this capture

| Route | Verdict | Status |
|---|---|---|
| Independent-port common-gain storage | bounds at least 1.4037e15 / 2.2313e13 | dead end, frozen |
| Two-occurrence scalar PE transport | 1 − omega_max·3 s/2 = −5.544984695 | dead end, sign unrepairable |
| Common-template common-gain frozen coefficients | H18 every-prefix fails at 2.70996 | dead end pending a new fact |
| Fixed SI diagonal metric | rho 1072.993 / 6659.686 | dead end, congruence-invariant |
| Fixed gravity/3-second diagonal metric | rho 11.183777 / 6.576082 | dead end, congruence-invariant |
| Separately scalarized information and remainder | −.380703 / −.257014 | rejected |
| BIAS2 conditional multiplier grid | no nonstrict master becomes strict | coefficient-feasibility failure |
| Source-centered motion storage | endpoint .999572659 / .932948382, prefix max 1.003245112 / 1.000001212 | open; its 6.41 / 2.27 chart budget is a scalarization artifact |
| Full signed vector information | margins .00042734 / .06705162, eta6 minima 972.647 / 862.261 | open; blocked on source admission |
| Declared-domain product box | velocity 35.358 / 9.357 | not invariant; the box is not a reachable set |
| Bias ball and template only | every motion group at most .2382 | retires the separated bias budget |
| Correlated covariance ellipsoid | word growth at most 7.2%; critical 4.974 / 27.452 sigma | open; needs covariance consistency |
| Forcing template alone | at most .0470 / .00472 | not the obstruction |

## Scope

These are point diagnostics on one frozen capture. The complete-word endpoint
ratios (.999572659 and .932948382) and prefix maxima (1.003245112 and
1.000001212) are unchanged by this experiment. Neither is the A21 bias-ball
growth a counterexample to anything: it is one more unmet premise. Nothing here is a uniform
nonlinear certificate, a legal counterexample, physical source admission, or a
promotion: P4 is unproved and P5 may not start.
