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
velocity, position, the integral state and the chart — the attained bound sits
within .37% of them in H18 and within 4.65% in A21, the widest of those being
A21's chart row. The gyro-bias row is wider still, 19.27% in H18. Every one of
these rows is at most .2382 of its declared ball on either bound, so the
enclosure width changes no conclusion; it is reported so the reader can see
which rows are tight and which are not.

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

## The best quadratic storage, in closed form

Every route above picks a metric and measures its ratio. The systematic
question is what the best metric could do, and it has a closed-form answer.
For `V(x) = x^T P x` the complete-word ratio is the squared `P`-weighted
operator norm of the word transition, and the infimum of that norm over all
`P > 0` is the spectral radius, so

    inf over P > 0 of rho_w(P) = rho(T)^2.

The infimum is attained when the eigenvalues of largest modulus are semisimple
and is otherwise approached, which the reported eigenbasis condition number
exposes. One metric is used before and after, which is what a uniform
statement needs; route 1 compares two different actual covariances and is a
different quantity, so its .999572659 / .932948382 are not comparable with the
numbers below.

| Word transition | rho(T) | best achievable rho_w | contracting storage |
|---|---:|---:|---|
| H18, 18 motion errors | .997663961 | .995333378 | exists |
| A21, 18 motion errors | .960813214 | .923162032 | exists |
| H18, full 21 states | 1.000000000 | 1.000000000 | none |
| A21, full 21 states | .999399433 | .998799226 | exists |

On the motion block the metric was never the obstruction: a contracting
quadratic storage exists for both modes, and the eigenbasis witness
`P = (S^-1)^H (S^-1)` attains the infimum. The two fixed diagonal metrics
frozen above at rho 1072.993 and 6659.686 were poor choices rather than
evidence against the architecture.

H18's full-state unit eigenvalue is structural. Its accelerometer bias is
unobserved over this capture, so it contributes an eigenvalue of exactly one
and no quadratic storage contracts the full word. That predicts both measured
bias-ball rows above -- H18 invariant to 2.2e-16, A21 growing 3.146e-4 -- and
it says the full state wants a bounded-bias statement rather than strict
contraction, which is what the retained BRMM hypothesis assumes.

These are the frozen capture's own transitions. Existence of a contracting
storage here is not a uniform certificate over words, sources or the nonlinear
coefficient dependence, and it promotes nothing.

## The bounded-bias theorem form, analytically

The unit eigenvalue above says strict contraction of the full state is the
wrong target. The right one is bounded bias, and it has a standard analytic
form once the word is split as

    x+ = A x + G b + r,    b+ = C x + Phi b.

H18 is an **exact cascade**: its bias-from-motion block `C` is `0` in every one
of the 2593 per-step factors, not merely in their product, because that mode
never updates the bias. Cascade ISS then applies directly. Choosing `P` to
attain `||A||_P = rho(A) < 1` and iterating over words gives, for any bias
ball `||b|| <= beta`,

    limsup ||x||_P <= (||G||_P * beta + ||r||_P) / (1 - ||A||_P).

| Constant | H18 | A21 |
|---|---:|---:|
| Per-step bias-from-motion defect | 0 exactly | 2.171e-3 |
| `rho(A) = ||A||_P` | .997663961 | .960813214 |
| `||A||_2` (non-normality) | 32.76 (32.8x) | 81.61 (84.9x) |
| `||Phi||` | 1.000000000 | .998728897 |
| `||G||_P` | 12.05 | 125.4 |
| ISS gain `1/(1-||A||_P)` | 428 | 25.5 |
| ISS limit at `beta = .4` | 3020.5 | 1317.5 |

Three things follow, and the third is the reason this is recorded rather than
promoted.

The theorem form is settled: bounded bias, not strict contraction, and A21
needs the perturbed form since its `C` does not vanish.

The severe non-normality, `||A||_2` about 33x and 85x `rho(A)`, is the single
explanation for the failed metrics and for every prefix ratio above one. A map
can contract asymptotically while amplifying transiently by that factor, which
is exactly what the prefix maxima 1.003245112 / 1.000001212 record.

The constants are loose. The same bias ball and the same word give a direct
reachable-set excursion of at most .2382 of any declared bound, against an ISS
limit of 3020.5, because ISS worst-cases the bias as adversarial and persistent
in the worst direction at every step and leans on submultiplicativity. The
analytic route supplies the architecture; the direct computation supplies the
numbers.

### What the lemma reduces P4 to

Its hypotheses are uniform, and this capture supplies them only pointwise. A
uniform certificate needs, over every admissible word `w`:

1. `rho(A_w) <= alpha < 1` **in a common metric** `P`, not a per-word one;
2. `||Phi_w|| <= 1`;
3. `||G_w|| <= G_max`;
4. `||r_w|| <= R_max`.

Only the pointwise versions are measured here, in per-word metrics of
condition 4.2e3 and 2.9e4. Pointwise spectral radius below one does not imply
uniform stability for a time-varying family, which is the classical gap and is
precisely where P4 sits. The lemma therefore localises the open problem to
these four bounds rather than closing it.

## Every route executed on this capture

The list below is the routes executed, not a claim that the space is
exhausted. The last row is the one that bounds the rest.

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
| Optimal quadratic storage, closed form | motion rho(T) .997663961 / .960813214; full 1.0 / .999399433 | metric was never the motion-block obstruction |
| Cascade ISS, analytic | H18 exact cascade, ISS limit 3020.5 vs direct .2382 | theorem form settled; constants loose; four uniform bounds open |

## Scope

These are point diagnostics on one frozen capture. The complete-word endpoint
ratios (.999572659 and .932948382) and prefix maxima (1.003245112 and
1.000001212) are unchanged by this experiment. Neither is the A21 bias-ball
growth a counterexample to anything: it is one more unmet premise. Nothing here is a uniform
nonlinear certificate, a legal counterexample, physical source admission, or a
promotion: P4 is unproved and P5 may not start.
