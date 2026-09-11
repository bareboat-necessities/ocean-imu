# OU-III proof research state

## Current hypothesis and controlling inequality

COMPLETE-BRMM describes wave displacement relative to a local equilibrium, not
absolute vessel translation. Its physical generator must supply one bounded
potential `phi` with `phi_dot=p`. The spectral alternative has uniformly finite
hard vector-amplitude moments, including the inverse-frequency moment; the
shaping alternative has an exact output identity and a hard state invariant.
`ou3_brmm_physical_wave_source.py` checks these alternatives in exact rational
arithmetic. They imply `S_L=phi-phi_L` and a derived all-time primitive diameter.
Neither a power spectrum, Hs, tuner cutoff nor the P4 working radius sets that
physical diameter. The full family's numerical generator envelope is still **E**.

The controlling P4 inequality remains the complete-word, same-history **joint24**
error/true-bias inequality, with compatible consecutive metrics, endpoint
`rho<1`, every literal prefix and first-exit retention. Generator qualification
alone does not prove it. P3 delta remains `1e-18`; P4/P5 remain false.

## Evidence and retained facts

- Exact harmonic output identities and shaping `LA=C`, `LB=D`, positive-definite
  metric and dissipativity checks derive `D_S=2 M_-1` or `D_S<=2 r ell`.
  Certificate fields are recomputed, not trusted. A derived 400 m*s example is
  deliberately not clipped to the 300 m*s working scale.
- The old exact indistinguishable `p=+/-d` proof remains a regression: **B under
  the old finite-window-only source**, **E as the intended source-specification
  omission**. For any finite derived D and nonzero component magnitude a,
  `h=(D+a)/a` gives `a*h-D=a>0`, excluding an indefinite DC displacement under
  the corrected source. Finite constant-position segments and zero waves are
  not excluded by a pointwise flag. No filter instability is asserted.
- One Live potential is carried across all words. The exact primitive transition
  and typed sample/event payload retain potential, p/v/S, source and generator
  identities. Endpoint/interval continuity is only a necessary outer relation:
  the continuous generator-to-p/v/a/IMU output attachment remains open.
- The trusted event attachment now consumes the existing same-signal JOINT
  frontend at the literal post-IMU/pre-magnetometer transition. Its H18 and A21
  prior-frequency regression succeeds without calling the obsolete
  measured-period-only frontend. It retains the physical-potential payload in
  every event cell; dropping it fails source binding. This is not a 601-sample
  cover or a finite-capture proof. The codec retains absent prior WPE states as null, not an invented period;
  the legacy canonical typed executor is still only a measured-period slice.
- Preserve the shipping initialization, one-time common S-origin reduction,
  BIAS0/1/2 separate drivers and true-bias ancestry, full P/H/R/K, Joseph/reset/
  projection splitting, candidate/active tuner and scheduler state. No `src/`
  change is part of this continuation. The qualified runtime Live/H18 handoff
  remains distinct from P4 capture. Conditional binary32 mathematics remains
  separate from deployment qualification.

## Current failure analysis / limiting quantities

**E — numerical physical source and output qualification.** Existing hard
position/velocity/acceleration caps do not determine a uniform inverse-frequency
amplitude budget. The analytical family `p_z=cos(t/n)`, n>=10, has uniform p/v/a
bounds but exact primitive diameter `2n`. It demonstrates missing uniformity of
those old premises, not a counterexample to a corrected fixed-budget source.
The limiter is the full family's hard generator envelope and its continuous
same-history output attachment. The current 300 m*s tube cannot qualify it.

**C — source-indexed measurement-Phi representation.** Before the sparse-identity
repair, both H18/A21 rebase smoke checks give `xi_exact=false`, while the identity
residual contains zero and there is no interior event. For the known exact
physical map `C=L=I`, `rho=0`, generic interval matrix operations unnecessarily
widen structural zeros in `xi=E_aw*(epsilon_new-epsilon_old)`. This invalidates
the generic representation, not the identity or the shipping filter. The
limiting quantity is equality of the exact zero rows, not an endpoint margin.
The sparse experiment passes in H18 and A21: all other rows stay exactly zero,
and uncertain epsilon differences remain outward (even for equal interval
boxes). Two regression tests pass. This fixes that representation failure but
does not close a numerical source-uniform rebase-defect or endpoint bound.

**C — coupled H18/A21 innovation inversion.** The refreshed two-sample
construction was invoked separately for BIAS0, BIAS1 and BIAS2. All stop on the
second sample's **H18 accelerometer Joseph update**, before A21 completion:
Gauss-Jordan pivot 1 is `[-6.741051843328836, 7.659816852449009]`. The entrywise
innovation hull admits the exact matrix `(1/10)*ones(3,3)`, which has identical
columns, rank one and zero determinant. Thus no universal inverse of that
independent hull exists. This is **C**, not an admissible singular innovation or
filter instability: the same-history Joseph/PSD and positive-R ancestry have
been lost in the rectangular representation. A midpoint-preconditioned residual
infinity-norm diagnostic is about 4.12, not below one. Another unstructured
inverse algorithm cannot remove an actual singular member of the hull.

**D/E — complete prior-root window attachment.** The codec now preserves absent
raw/log periods as explicit nulls and rejects a forged usable latch. The
same-signal event attachment supports prior-frequency evolution. However, the
retained legacy canonical typed executor still calls the measured-only frontend;
a 601-sample source-uniform prior-root window is not certified by these tests.

**F/G and open P4 obligations.** No source-uniform worst endpoint margin,
every-prefix gain, maximum retained basin, finite H18 capture time or H18/A21
basin landing is certified. Joint24 generator forcing still needs to enter the
same-source augmented master, and every consecutive metric/hybrid transition
needs its actual compatibility bound.

## Failed approaches / DEAD_ENDS

The old finite-window-only all-18 theorem and the search that chooses physical
D_S from `[0,300]` are not valid routes. Sign symmetry does not make arbitrary
oscillatory `p` and `-p` sensor-indistinguishable, so it supplies no general
D_S<=300 necessity theorem. A21 marginal-motion 18-state storage discards the
needed motion/bias cross information. Independent coefficient rectangles,
packetwise S resets, covariance-as-hard-membership and replay-fitted source
constants remain disallowed.

## Independent critic / alternatives

The strongest reason to abandon a numerical P4 search now is the absence of a
qualified family-wide physical generator envelope: a margin for an invented
D_S would answer the wrong theorem. Distinct source alternatives are a hard
RAO amplitude-measure envelope, a bounded driven shaping-state invariant, or an
explicitly separated moving equilibrium with a qualified wave potential.
Choose using physical provenance, not the most favorable P4 radius.

For the repeated structural-zero rebase issue, alternatives are (1) exact sparse
identity embedding, (2) symbolic DAG normalization throughout the finite map,
and (3) unshifted joint coordinates avoiding this rebase. The first is justified
by the newly isolated exact identity, not by a subdivision guess: it removes
all zero-product widening algebraically. It makes no prediction about rho.

For this inverse limiter, the critic rules out repeating raw LDLT, generic
Gauss-Jordan or midpoint preconditioning of the same box. Distinct alternatives
are retaining the Joseph Gram/square-root factors with the same P/H/R graph,
retaining a joint information/Schur covariance set, or proving a source-family
subdivision actually removes the spurious singular member. The first two retain
PSD ancestry; merely relabeling the same intervals does not. A new experiment
must certify a complete useful word, not only replace the inverse exception by
an enormous norm bound.

## Next falsifiable experiments

The sparse H18/A21 rebase and prior-frequency/event ancestry regressions pass.
Regenerate the corrected source/final gate for the current commit, then qualify a physical family envelope without consulting the P4 domain, bind
its continuous generator to the existing source relation, and report the first
complete H18/A21 same-history joint24 endpoint/prefix failure. Only a useful
complete-word margin warrants outward enclosure refinement. Numerical source
qualification, literal prefix retention, capture and hybrid transport must all
close before the indefinite theorem can promote.
