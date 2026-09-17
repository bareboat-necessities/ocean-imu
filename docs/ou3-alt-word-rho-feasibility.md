# ALT complete-word rho feasibility

The controlling object of the ALT theorem is the complete-word ratio

```
rho_w = V(F_w(x)) / V(x)
```

for the coercive joint24 storage `V` with the admitted bounded neutral supply on
`[e_ba, beta_true] = 18:24`. This note records its first measurement.

The result is a **theorem failure**: an admitted legal word has `rho_w >= 1` for
every common storage. No sharper enclosure, multiplier or metric can change it.
`tools/stability/ou3_alt_contraction/finite_word_rho_diagnostic.py` is
authoritative and never promotes a gate.

## Why the floor is metric independent

`common_storage_master` reduces the existence of a finite bounded-coordinate
supply completion to the projected Finsler restriction

```
Z' (A_w' M A_w - rho M) Z < 0,      Z injects ker(C).
```

The floor that holds in **both modes** comes from an invariant subspace. If
`V` is contained in `ker(C)` and `A_w V` is contained in `V`, then for every
`x` in `V` the restriction reads `(A_w x)' M (A_w x) < rho x' M x` with
`A_w x` again in `V`, so

```
rho > spectral_radius(A_w restricted to V)^2
```

for every storage `M`. No multiplier, metric or sharper enclosure can go below
it, and the infimum is not even attained when the restriction is not
diagonalizable.

**H18 additionally has an invariant kernel.** Its accelerometer bias is held, so
no literal event writes a motion coordinate from a neutral one and the word map
is block lower triangular,

```
A_w = [[A_mm, A_mn], [0, A_nn]].
```

The diagnostic checks this entrywise rather than assuming it; the worst
neutral-to-motion entry is exactly `0`. Then `ker(C)` is itself invariant, the
restriction collapses to `A_mm' M_mm A_mm < rho M_mm` on the 18-dimensional
motion subspace, and `spectral_radius(A_mm)^2` is a floor.

**A21 does not.** Its released accelerometer bias is corrected from motion
coordinates, so the neutral-to-motion block is nonzero, `ker(C)` is not
invariant, and the A21 motion-block spectrum is reported as an indicative ratio
rather than a floor. This is why the falsification below is stated through the
invariant pair, which is checked at full joint24 width and holds in both modes.

## The ungauged heading obstruction

For an ungauged word -- one containing no magnetic event -- the pair
`(theta_z, bg_z)` is an **exactly invariant** subspace of the full composed
joint24 map, with every entry of the surrounding rows and columns identically
zero, in H18 and A21 alike. Both coordinates are motion coordinates, so the pair
lies in `ker(C)` and the lemma above applies. It carries

```
[[1, T],
 [0, 1]]
```

for the elapsed word horizon `T`. The algebra is elementary and does not depend
on word length:

- the accelerometer residual sensitivity `-skew(f)` annihilates the specific
  force direction, which is the vertical in the quiet member;
- the S=0 rows touch only `S`;
- the attitude/gyro-bias prediction advances `theta_z += h * bg_z` and leaves
  `bg_z` constant;
- the seeded covariance keeps zero cross-covariance between that pair and every
  observed coordinate, so the Joseph gain rows for both are exactly zero and no
  correction ever leaks in.

A product of such unipotent blocks is unipotent, so a word of any length has
restricted spectral radius exactly `1` and `rho_floor = 1`. The diagnostic
exercises both modes and both literal sample shapes, with and without the
scheduled S=0 event, and certifies every prefix.

This word is admitted, not contrived. The ungauged timeout path leaves
`mag_ref_set_` false, so `updateMag()` is never reached before north lock, and
`MAG-CALL-SCHEDULE-v1` imposes post-gauge service rather than a pre-gauge
acquisition deadline. The magnetometer is asynchronous with no declared ODR, so
zero magnetic events inside one 3 s word is legal. The Live covariance seed
itself records the situation: the ungauged attitude seed "belongs to the gravity
quotient until magnetic regauge".

## Measured floors on the canonical 3 s word

600 IMU samples, shipping Live-entry covariance seed, committed tuner schedule,
S=0 at the committed pseudo-period, accelerometer on every valid sample:

| Source phase | rho floor | Limiting direction |
| --- | --- | --- |
| quiet ungauged | `1` exactly | `theta_z`, `bg_z` |
| wave ungauged | `0.99954648` | `v_y`, `p_y` |
| quiet gauged (25 Hz magnetic) | `0.99597048` | translation |
| wave gauged (25 Hz magnetic) | `0.99604593` | translation |

The composed word is a product of hundreds of binary64 Cayley/Joseph operations.
The deviation of the certified unipotent diagonal from its exact value `1` is a
direct drift proxy and is `2.3e-8` at 600 samples, so the measured floors are
meaningful only well outside that. The falsification above rests on the exact
invariant-subspace algebra, never on a measured radius near one.

## What this does and does not invalidate

Falsified: strict `rho<1` for a common coercive joint24 storage whose only
bounded supply is `[e_ba, beta_true]`, over the admitted Live word language.
Adding subordinate lemmas cannot recover it.

Not invalidated:

- the shipping filter. An unobservable heading under gravity alone is the
  expected physical situation, not an instability;
- motion ISS with a supply, or a theorem stated on the gravity quotient;
- the gauged word, whose floors above leave real room: about `4.0e-3` distance
  to one per 3 s word;
- the independent P2/P3/P4/P5 route, which is untouched;
- the eleven finite-master qualifications, which remain open on their own terms.
  They are subordinate to `rho_w` and finishing them would not reach the
  declared contraction.

## Current limiting quantity

For the ungauged word, `spectral_radius(A_mm) = 1` on `(theta_z, bg_z)`: an
exact obstruction with no margin to improve.

For the gauged word, the limiter moves to weakly observable surge/sway (`v_y`,
`p_y`), at a floor of about `0.9960` per 3 s word. That margin is the one a
rigorous enclosure would have to survive.

## Next falsifiable experiments

The formulation must change before any further enclosure work. Three
qualitatively different routes, none a refinement of the failed one:

1. **Quotient the gauge.** State the theorem on the gravity quotient, with
   storage coercive only transverse to the ungauged heading direction, and prove
   the quotient map is well defined across the literal reset/projection
   branches. Falsifiable: exhibit a legal word whose quotient map still has
   spectral radius one.
2. **Make gauging a hypothesis of the Live word.** Declare a magnetic service
   class for Normal Live, as the non-ALT route already does for
   `MAG-CALL-SCHEDULE-v1`, and re-measure. Falsifiable directly: the gauged
   floors above are already below one, so the question becomes whether the
   `0.9960` margin survives interval enclosure over 600 steps -- measure the
   enclosure width before building it.
3. **Supply the heading pair.** Move `(theta_z, bg_z)` into an independently
   bounded neutral port, as `[e_ba, beta_true]` already is. This needs an
   independent bound on heading error that is not the desired conclusion;
   charging unknown motion error as bounded supply is forbidden and would be
   circular. Falsifiable: produce the independent bound, or show none exists in
   the admitted language.

Route 2 is the cheapest to falsify and route 1 is the only one that keeps the
current ungauged source language. Route 3 is listed for completeness and is the
most likely to be circular.

## Running it

```sh
PYTHONPATH="$PWD:$PWD/tools/stability" python3 \
  tools/stability/ou3_alt_contraction/finite_word_rho_diagnostic.py \
  --samples 600 --output /tmp/ou3-alt-word-rho.json
```

`finite_master_guard` carries a short memoised probe of the same measurement in
`complete_word_rho_feasibility` and lists the falsification under
`falsified_prerequisites`. All ALT gates stay false.
