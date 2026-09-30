# OU-III proof research state

## Current hypothesis

One persistent MARINE MOTION / IMU BIAS / MAGNETIC SERVICE execution follows
construction -> capture -> informed H18 -> refinement/release -> recurring
A21 -> regional practical stability. No estimator, physical assumption,
quality gate, tuner law or carried diagnostic is changed. Contraction words
are 100 s; 17 s remains nuisance/root warm-up. All history coordinates,
covariances, frontend state, actual event ordering and committed
(tau,sigma_aw,R_S,T_S) values are carried, not restarted.

The active analytical calculation is [the linked full-covariance soft
return](ou3-linked-soft-return.md). It retains the complete same-word Pi,
quotient contribution, real next-word information and cross-covariances.
The earlier detailed research record is preserved verbatim in
[the superseded ledger snapshot](ou3-proof-research-state-before-linked-return.md).
That record is supporting history, not the current theorem status. In
particular, its claims of automatic kernel invariance from next precision,
and LP10/LP14--LP22 of ou3-aw-adjoint-cancellation.md, are superseded below.
The LK/AL/DC/DP exact scalar reductions in ou3-corrected-word-proof.md remain.

The complete nonlinear stability theorem is OPEN. This calculation is an
analytical reduction of O2, not a closure of capture, covariance invariance,
source-uniform contraction, finite-error retention or float32 transfer.

## Evidence

- Exact scalar return: G_c=J+nu nu'/c, w=Phi_tilde' n and
  D_W(c)=n'Pi n+w'G_c^-1 w. Its dual is a rank-one generalized eigenvalue;
  D_W(c)<=c_next is exactly the block PSD test LS3 in the linked note.
- General Schur coordinates give D_W(c)=d_perp+c ell^2/(s+c j),
  s=||nu||^2, with d_perp retaining BOTH Pi and diffuse quotient uncertainty.
  The exact fixed-point polynomial is
  j c^2+(s-ell^2-j d_perp)c-s d_perp. For j=0 and d_perp>0, a finite
  ceiling requires ell^2/s<1. These reductions are consistent with the
  existing DP/DC argument; finite gain alone is not a substitute.
- New full-baseline identity: conditioning B+a yy' on REAL information J
  gives C+a(Ly)(Ly)'/(1+a y'Hy), where C=(B^-1+J)^-1,
  L=I-CJ and H=J-JCJ. No next proof precision is inserted as an observation.
- New sharp linked maximum on u in Range(E), ||u||=1:
  n'Cn+p v'(E'E+p E'HE)^-1 v, v=E'L'n. It keeps the same baseline and
  action together. For E=I it equals n'[(B+pI)^-1+J]^-1 n exactly for
  this one functional. An enlarged direction set is only a relaxation.
- The correlated kernel/quotient plane has sharp total maximum
  b_nn+p-b_nq^2/(b_qq+p+1/j). The worst angle is generally intermediate,
  not coincident. The covariance decrease and apparent excess amplification
  must be combined before a bound is taken.
- Twelve new exact-rational tests pass. The existing word-diameter suite
  (6 tests) and corrected-word suite (12 tests) pass. No numerical example
  is promoted to a shipping-reachable or source-uniform certificate.

## Current limiter

The linked return is not yet bounded by a common finite c on every admitted
same-history adjacent 100-s pair. O1's compactness/boundary-nullity and
ordered-gap premises must retain their established scope; they do not supply
O2. On a fixed compact class the established large-c argument still requires
strict transfer ell^2/s<1 on the exact j=0 subset, or an actual blockwise
replacement. Finite deterministic gain is insufficient.

Capture, compatible-class stationary practical stability, H18/release
retention, constructive O1 constants, coercivity where required, nonlinear
radius, every-prefix retention and implementation/arithmetic totality remain
separate open obligations. No physical counterexample to the estimator is
claimed by the matrix examples.

## Failed approaches / DEAD_ENDS

1. **Treating rank-one information as rank-one covariance.** Failed quantity:
   claimed next scalar variance from a bare rank-one prior. Classification:
   mathematical model mismatch. Invalidated hypothesis: the known-root Pi
   and diffuse quotient may be discarded. Exact example Pi=I, Phi=I,
   J=diag(0,1), nu=n=e1 gives D(c)=1+c for every c. Retain Corollary K and
   the full Schur covariance expression; next test is LS3, not that shortcut.
2. **Using next ceiling precision as a measurement.** Failed quantity:
   claimed automatic kernel ceiling from 1/c_next. Classification: circular
   invariance argument. In the preceding example the synthetic update gives
   c/2 although the real return is c+1. Retain the scalar update only as
   synthetic algebra; actual observation information must be kept separate.
3. **Compressing Pi before inversion.** Failed quantity: LP10's compressed
   condition-number bound. Classification: lost cross-covariance. For
   Pi=[[1,4/5],[4/5,1]], compression to e1 has condition number 1 but the
   linked product is 25/9. Retain the Schur-shortened inverse, not A^-1.
4. **Transferring coincident-line maximality to a full covariance.** Failed
   quantity: angle monotonicity with background correlations. Classification:
   omitted background. B=[[1,-.9],[-.9,1]], j=p=1 has worst total 1.73,
   versus aligned 1.595. Retain LS8 and the linked LS12 total cancellation.

Earlier failed mechanisms and their retained derivations are preserved in
the superseded ledger snapshot; none is reactivated by this calculation.

## Retained facts

The source-faithful Riccati, Joseph, reader, physical-history and nuisance
identities retain their previously qualified scope. The scalar prior is an
upper covariance premise implemented as a lower INFORMATION bound. Actual
process, accelerometer, S and magnetic sources remain jointly represented.
No additional physical excitation, event schedule or covariance ceiling is
assumed. O2 and the end-to-end theorem remain false/open in status outputs.

## Alternatives

Use LS8 for the literal carried soft direction, LS9 for its allowed propagated
subspace, or LS3 directly. At positive eigenvalue multiplicity retain the
whole actual least eigenspace; nullity<=1 does not limit its multiplicity to
two. Full-space maximization is a rigorous relaxation, not reachability.

## Next falsifiable experiment

Analytically substitute the literal homogeneous compatibility transport and
its full background covariance into LS3/LS8. Prove a uniform linked ceiling
or identify the exact same-history failure, with a block return if necessary.
Do not multiply independent dbar/Hbar extrema or add next-root precision.
The immediate algebra audit is the exact-rational linked-return test class.

## Validation limitations

The inherited PR contained literal escaped-newline syntax errors in
`theorem_status.py` and `magnetic_strata_certificate_cli.py`; repairing only
those textual errors allows theorem tooling to compile. No theorem flags or
quality thresholds were changed.

The complete rank-loss test module runs 37 tests: 35 pass, including all 12
new tests; two inherited interval tests error. The exact errors are
`ValueError: rectangular midpoint/radius matrices required` for an empty
source factor, and `ArithmeticError: singular midpoint` for the expected
fail-closed magnetic inverse. The theorem-status module runs 8 tests: 6 pass;
committed JSON differs from status_report(), and an older test still demands
true for a retracted S-chain claim. Those are existing implementation/status
contract failures, not failures of LS1--LS12. They remain to be reconciled
without re-promoting retracted mathematics. Full native CI was not run here.


### Literal graph alpha / block refinement

The exact-kernel slope is now tied directly to the homogeneous compatibility
graph, not to an arbitrary soft eigenvector. For
`r_W=(a,0,-A_Wa)`, fixed physical metric M and same-history successor
`r_+`,

`alpha_W=|r_+' M T_W r_W|^2/[(r_W'Mr_W)(r_+'Mr_+)]`.

This follows because `J_Wr_W=0` annihilates both the Schur diagonal and
cross block, so the LS4 coupling is exactly the endpoint graph pairing.
`A_W,T_W,r_+` retain the literal physical history, covariance-generated
gains, S/magnetic/accelerometer chronology and committed coupled tuner values.

No current lemma proves the one-pair supremum is <1. This is not a proof
failure by itself. For m words, compose the FULL Riccati word before
scalarization and use LS1 on that superword. On an exact block kernel the
large-c slope is the same endpoint graph formula with
`T_[j,m]=T_(j+m-1)...T_j`; if an intermediate word charges the carried
mode, the block has positive Schur information and zero large-c slope for
that mode. Along a wholly exact persistent chain the block alpha is the
product of the boundary alphas. Therefore one-word alpha=1 is harmless if
unit transfer cannot persist through a finite block.

The existing equality-chain compactness theorem now applies to this literal
quantity: absence of an infinite unit-transfer same-history execution implies
existence of finite m,delta with block alpha<=1-delta. The unresolved question
is global continuation/escape of the coupled compatibility zero dynamics. No
estimator assumption or tuner law was changed.


### Literal zero-dynamics substitution result

Substituting r=(a,0,-Aa) into the actual Live chronology gives the graph
cocycle `A_+ T_att a=T_ba A a` plus the literal magnetic/accelerometer
compatibility rows.  Prediction, due S=0 update, accelerometer correction,
measurement-only tuner update, one-sample-later tau/sigma/R_S commit, AW
covariance sync and asynchronous magnetic callbacks are retained in order.

On a regular stratum the two transverse compatibility equations solve the two
transverse accelerometer-innovation components by IFT, leaving one
longitudinal physical input.  The tuner is a bounded lagged coefficient
sequence, the S residual is endogenous, and AW covariance sync is mean-neutral.
No term diverges or has a sign forcing escape as BA q->0; q=0 is a regular
candidate compatibility manifold.  Bounded physical primitives likewise do
not force escape because the required physical forcing can be zero-mean.

Conclusion: existing tuner/S/BA/physical inequalities do NOT imply finite
escape.  Local compatible continuation is proved conditionally wherever the
transverse compatibility Jacobian and physical/service/gate margins are
strict.  Global forward completeness remains OPEN because no source-uniform
lower bound on that Jacobian or invariant magnetic/gate margin is proved.
The next decisive obligation is a compact invariant strict-margin ZG patch or
a theorem that every constrained ZG trajectory loses one of those margins.
