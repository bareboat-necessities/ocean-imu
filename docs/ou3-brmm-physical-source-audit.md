# Physical BRMM source audit

## Coordinates and the missing premise

The wave coordinate is `x_CoG = x_equilibrium + p_wave` in fixed inertial axes,
with `v=dot(p_wave)` and `a=dot(v)`. Global origin, current, propulsion, leeway,
secular wave drift and slow equilibrium translation belong outside `p_wave`.
Baseline acceleration must be zero, compensated, or retained in a declared
same-history disturbance. This is not a software position reset or removal of
an actual acceleration from the IMU.

The audit starts from main `f7123bd874f0adeb508b9308b1586e7a39264da6` (PR #515).
The root Makefile uses the 28 ft vessel RAO dataset from
`oceanography-waves-lib` tag `v1.2.1`. Its `VesselRao.h` constructs displacement,
velocity and acceleration from the same fixed complex harmonic responses. Its
mean-pose convention excludes drift; its general API nevertheless accepts zero
frequency. API acceptance and physical theorem admission are different tests.

| Intended physical fact | Existing executable premise | What was missing / correction |
| --- | --- | --- |
| No acceleration DC | A bounded velocity primitive and `a=dot(v)` | This eliminates acceleration DC, not displacement DC. |
| Oscillatory displacement | Position cap, derivative chain and short-window moments | A constant integration offset still satisfies them. Require a physical output/potential realization. |
| Zero displacement DC | Not encoded as a hard all-time consequence | Derive it from the bounded potential; do not use asymptotic zero mean as the proof. |
| Nonzero minimum frequency | Tuner limits and RAO shaping cutoff existed | Neither is a hard physical low-frequency cutoff. Qualify physical support or a finite inverse amplitude moment. |
| Finite amplitudes | Hard p/v/a caps and statistical Hs/spectrum descriptors | They do not determine a hard spectral total-variation budget. PSD is not a pathwise certificate. |
| Bounded shaping states | Finite-window output/moment representation | Add an exact potential output identity and all-time hard state invariant; the full numerical family remains unqualified. |
| Stationarity / parameter changes | Same-history lineage and finite-window variation/phase constraints | Independent window resets or local rate bounds can accumulate indefinitely. Carry one potential and prove invariant/jet transport at switches. |
| Translation/reference removal | Vessel-response mean-pose convention | Encode the baseline split in the source contract and retain any baseline acceleration as a disturbance. |

For the old source, `p=d, v=a=0` has bounded p and, on each fixed short window,
`Delta S = h*d`. No premise limited its sum over all windows. This explains the
oversight without asserting instability of the physical estimator.

## Hard source correction and derived bound

The paper's physical-generator section and lemma are authoritative. The exact
backend `ou3_brmm_physical_wave_source.py` provides two sufficient realizations:

- Hard vector amplitude measures with fixed uniform moments `M_j`, j=-1,0,1,2.
  The same measures generate p/v/a and
  `phi=integral (sin(wt)dA-cos(wt)dB)/w`; hence `phi_dot=p` and `D_S=2 M_-1`.
  With a physically qualified lower frequency and hard amplitude mass,
  `D_S<=2 M_0/w_min`. Finite support is sufficient but is not required when the
  inverse moment is already finite.
- A bounded physical shaping state `xdot=A*x+B*u`, `p=C*x+D*u`, `phi=L*x`.
  Exact `LA=C`, `LB=D`, a positive-definite metric, a dissipativity inequality
  and hard initial/input bounds yield an invariant radius r and
  `D_S<=2*r*sqrt(trace(L*P^-1*L'))`. Exact rational LDL and outward rational
  square roots check this certificate. The current checker is constant-matrix;
  more general modulation/switching needs its own transported invariant proof.

Both imply `S_L(t)=phi(t)-phi(t_L)`. The actual Live origin is taken once.
Neither position nor S is reset at word boundaries. A nonzero constant p would
force its integral past the derived finite D, so the old witness fails physical
admission by theorem. Zero waves and finite flat segments are not excluded by
an ad hoc pointwise rule.

The exact old paired-history regression remains in
`ou3_brmm_infinite_continuation.py`: B for the old finite-window source, E for
the intended physical omission, and excluded under the corrected definition.
The new regression explicitly checks
`constant_nonzero_position_zero_velocity_history_admitted == False`.

## What this does not qualify

There is not yet a hard numerical generator envelope for the entire declared
family. In particular, `p_z=cos(t/n)` has uniform old jet bounds but primitive
diameter `2*n`. Individually positive frequencies do not give uniformity across
n. This is an E gap in the old premises, not a corrected-source instability.
The 300 m*s working radius cannot supply the missing physics. The backend accepts
an independently derived 400 m*s example and rejects rewriting its certificate
to 300. No physical maximum or working-tube comparison is claimed for the full
family until its hard constants are qualified.

The hard provider now requires a recomputable generator certificate and one
potential/Live-potential ancestry across words. The exact boundary bridge and
typed event payload retain that relation alongside the existing p/v/S moments.
Those necessary outer constraints do not prove that every interval sample is
an output of the continuous generator. Code-owned numerical-envelope and
continuous-output qualification gates therefore remain false. The source-cover
contract also leaves physical-generator forcing consumption by the joint24
augmented master false. Flags supplied by an artifact cannot override them.

## Preserved execution and tested continuation

No shipping header or tuning constant is changed. P3 delta is `1e-18`. The real
same-signal JOINT frontend now owns the event attachment even before WPE has a
usable period, at the literal post-IMU/pre-magnetometer stage. Candidate/active
commit order and current P/H/R/K remain checked. Absent prior WPE states are
serialized as null, not a manufactured period. The complete prior-root window
cover is still open.

The source-indexed measurement-linearizing Phi rebase is distinct from the
physical wave potential. Its known identity map now uses the exact sparse
`E_aw*(epsilon_new-epsilon_old)` representation, preserving uncertain epsilon
and exact zero rows. This repairs a representation failure, not a filter.

The joint24 target retains motion/bias cross terms, BIAS0/1/2 separately,
true-bias ancestry, Joseph/reset/projection splitting, finite precision and
consecutive metric compatibility. P4 endpoint, every-prefix, first-exit,
maximum-basin, finite-capture and release obligations are not closed.
See the research ledger and generated fail-closed gate for the current limiter.

## Reproduction and evidence

Run the physical-wave-source, old infinite-continuation, centered-S, hard-provider,
codec, source-indexed rebase and source-uniform event-attachment tests under
`tests/validation`. The physical-source CI workflow regenerates the source and
P3/final-gate JSONs from the checked-out commit. The source-cover-construction
workflow preserves the original negative witness and native Live-entry audit.
`docs/ou3-brmm-physical-source-evidence.json` records local results and exact
input hashes; it is not a source-uniform stability certificate.

Local `make all` was attempted and stopped at missing Eigen headers. The full
article reached an absent generated PGF chart; the edited six-page BRMM section
was separately compiled and rendered with its IEEE preamble. Neither limitation
is described as a passing full build. Remote CI remains authoritative for its
own checked-out commit and installed dependencies.
