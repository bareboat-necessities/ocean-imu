# OU-III proof research state

## Current hypothesis

BRMM is the primary physical source: bounded accelerometer-bias error and
practical ISS of the other 18 errors, using the full 21-state estimator.
Keep the common physical history, all actual anisotropic R_S corrections,
full Q/floors/cross covariance, finite reset/projection, physical -S_true,
and the closed .4 bias-estimate ball. BRMM-0 requires a_m=dv_m/dt and bounded
same-history velocity; fixed physical acceleration DC is excluded. Constant
sensor offset remains BIAS0. Bounded velocity alone does not bound p or S.

## Current dataset and evidence

All current replays use oceanography-waves-lib v1.2.1,
`sim-data-files-vessel-rao-28ft.zip`, SHA256
`6d6eb92db78e97f1d2456c6b92387993bf23be14a9a6cf26943f0a22fc6fee60`.
Its generator commit is e442150682f560384be427df4cc7815956a091c5.
It applies an estimated stationary 28 ft fin-keel sailboat RAO to incident
components and emits exact CG kinematics and body-rate derivatives. Incident
spectra and filenames retain their incident-sea meaning. This engineering
preset is not hull-specific measurement or a continuum BRMM certificate.

Source/runtime audits and the conditional proof are being regenerated.
Until replaced, checked-in source/runtime JSON describes its own recorded
archive identity and must not be relabeled as vessel evidence. All phases
and all violating samples must remain in the audits.

The source-connected PM 1.5 m probe now uses this RAO and one common incident
root. Independent Python reconstruction agrees with the C++ source to
5.22e-15 (H18) and 1.42e-13 (A21). Both fixed 600-sample windows satisfy the
checked point caps/branch premises. The actual forced endpoint ratios are
7.2745201904 (H18) and .4378846487 (A21); the A21 prefix maximum is 1.297748031.
These ratios include nonzero physical forcing and are not homogeneous P4
contraction tests or counterexamples. No source or direction search was used.

## Conditional P3 scope

The established matrix implication uses delta=1e-18 for H18 and A21 under
explicit configured Normal-Live premises: acceleration 4 m/s^2, body rate
30 deg/s, Racc=.04 I, Rmag=.09 I and accepted-vector geometry/recurrence.
The retained H18 outward LDLT pivot is 4.987499868870966e-14 and first active
A21 bias margin 1.2499987189052501e-9. The new rerun must report its own revision.
These are conditional matrix facts, not full physical admission, nonlinear
projection stability, P4, or P5. BRMM recurrence alone does not imply PE.
Actual replay Racc=.000866618473 I and Rmag=3.68640018 I remain different.
The RAO migration does not silently change those premises or the filter.

## Failure analysis and current limiter

The old claim that startup exclusions alone reconcile surface-reference
histories with the configured caps failed source admission, not stability.
The new RAO source audit still observes cap/impulse violations in larger
incident seas. Qualification must use complete runtime phases, actual
accepted vectors and actual measurement matrices before widening P3.

The forced H18 ratio above one invalidates interpreting forced endpoint
energy as homogeneous contraction. It does not invalidate ISS with a supply
term, P3, or the filter. The limiter remains a useful same-history practical
supply/storage bound with every-prefix retention and the actual nonlinear
projection. A favorable forced A21 endpoint cannot discharge that obligation.

CI's first full-evidence attempt failed at `--repo: command not found`
(exit 127) after the archive checksum had passed. Class: workflow migration
implementation defect; leftover command continuations were removed. This
invalidated execution of that attempt, not archive identity or proof algebra.
The next falsifiable check is the unchanged full replay pipeline at the fixed
workflow revision, with genuine provenance rather than fingerprint restamping.

## Retained facts / DEAD_ENDS

- Independent-port common-gain P4 storage was quantitatively unusable:
  point factors .9999569976489486/.9788191291615017 and gains 2^32/2^36 gave
  bounds at least 1.4037e15/2.2313e13. Do not refine that discarded coupling.
- Two-occurrence scalar PE transport cannot cover the proposed 250 deg/s,
  one-second recurrence: 1-omega_max*(3 s)/2=-5.544984695. Interval tightening
  cannot repair that sign. Those surface-derived cap proposals are unfrozen.
- Bias qualification plus the .4 estimate ball bounds bias error; BIAS0/1
  need qualification and optional BIAS2 needs its actual nonlinear sector.
- Passive runtime instrumentation must compare byte-for-byte to its baseline;
  matched -ffp-contract=off resolved the prior native FMA instrumentation defect.

## Critic pass, alternatives and next experiment

The strongest reason to abandon the old architecture is its loss of physical
source correlation and uninformative common gain. Retain full matrix/group
structure and compare: (1) common latent-increment/p/S primitives with explicit
channel budgets; (2) full transported accepted-vector information instead of
two-occurrence scalar transport; (3) a separately qualified recurrence/response
premise from deployment. None may be fitted to replay extrema for PASS.

Finish the RAO phase/source audits, rerun canonical P3 without changing its
configuration, and report conditional proof versus physical admission
separately. P4 needs useful practical bounds and every-prefix retention before
rigorous source covering. P5 remains blocked by the unclosed P4 obligation.
