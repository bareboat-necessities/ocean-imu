# OU-III stability handoff

Continue only in PR #560, branch `ou3-explicit-regional-practical-stability`.
Read AGENTS.md and the current research state before editing. Shipping source
is authoritative; do not change behavior, physical assumptions, or gates.

The only path is construction, finite history-dependent capture, finite H18,
finite magnetic-reference refinement/release, recurring informed A21, and
regional practical stability. Carry the same physical execution and estimator
state through every boundary. No new common startup deadline is assumed.

## Current result

Read `ou3-literature-applicability.md` for the primary theorem/hypothesis map
and the realized-gain auxiliary-observer contract. Read `ou3-factor-loss.md`
for the factor derivation, exact singular nuisance elimination and source
export reproduction. Neither an IEKF theorem nor the 2025 modified SO(3)
filter theorem is applied to unchanged OU-III.

The complete root-whitened factor engine now agrees with the dense rational
reference on coupled full21 and adversarial words. A carried 16-s quiet-water
shipping-core export has diagnostic rho=.8327094089, with source-double parity
below 5e-14 relative and a matching 70-digit prefix check. This establishes
feasibility of the representation on that trace, not a uniform theorem margin.
The rational residual bound on a stored floating QR matrix deliberately leaves
whole-word factor error, source uniformity and float32 transfer unverified.

The complete rational 4x4 LIN endpoint-action matrix is now certified for the
bound default profile. The explicit comparison is
`P_end >= E_LIN A^-1 E_LIN'`. Its singular full-state form is intentional.
`lin_matrix_certificate.py` preserves all within-LIN factors; raw neutral
S-information normalization is at least 3.59e-11 over the event-time boxes.
This number is not a full-state rho.

The constructive joint lower covariance is now closed at regular A21
post-prediction roots after the 16-s window. Combine P>=diag(q_AG I,0,q_BA I)
from fresh process injection with P>=E_LIN A^-1 E_LIN^T using equal convex
weights. The resulting full 21-state matrix is strictly positive without
discarding cross covariance. AG is 4.9991e-10; BA is at least 4.99999e-10;
the LIN block remains the matrix A^-1/2. Uniform upper bounds, prefix
retention and float32 transfer remain open.

The earlier reported mu_N and rho0 were overpromoted: a restricted magnetic
heading/bias Gramian cannot be lifted to independent full-state heading
information. The exact counterexample in `word_energy.py` has restricted
information I_2 while a nuisance-cancellation direction preserves energy.
This falsifies the implication, not the shipping filter or its physical domain.

The correct complete-word identity is
`M_end' P_end^-1 M_end + D_word = P_root^-1`, with process and actual-innovation
loss terms. To obtain rho0<1, certify `D_word >= delta P_root^-1` source-uniformly
in all 21 coordinates, with recurring covariance coercivity. The root transport
includes every preceding correction/reset and cannot be replaced by a bare
attitude/gyro kinematic map.

## Required continuation

1. Use the literature hypotheses and factor/range construction to bound full
   loss and its cross blocks, or an equivalent joint full-state path-action
   comparison. Include stable a_w/BA directions and complete varying histories.
2. Use the new joint root lower bound and certify uniform upper covariance
   and every-prefix bounds through literal correction/reset cadence. A
   smallest LDL pivot is not an eigenvalue floor; the new lower bound uses
   two full Loewner comparisons and their convex combination.
3. Close the nonlinear projection/reset/tuner remainder and explicit radius.
4. Compose float32 supply, covariance coercivity, every-prefix retention, and
   finite capture/H18/release entry into that radius.
5. Complete physical qualification and recurring applied-service certificates.

Do not revive entrywise Riccati subdivision, the failed endpoint-batch lemma,
restricted-service lifting, or obsolete parallel theorem architectures.

## Reproduction

```
python3 -m py_compile tools/stability/ou3_theorem/*.py
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json
python3 tools/stability/ou3_theorem/block_factor_numeric.py
make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test
./tests/kalman_ou_iii/shipping_contract-test
./tests/kalman_ou_iii/shipping_transition-test
```

The block assembly command reports partial evidence successfully while keeping
`verified=false`, `rho0=null`, and `theorem_closed=false`. Passing reproduction
must never be interpreted as proof completion. Generated matrix/audit evidence
is recomputed and compared by the evidence builder and bound by provenance.
