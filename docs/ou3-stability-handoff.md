# OU-III stability handoff

Continue only in PR #560, branch `ou3-explicit-regional-practical-stability`.
Read AGENTS.md and the current research state before editing. Shipping source
is authoritative; do not change behavior, physical assumptions, or gates.

The only path is construction, finite history-dependent capture, finite H18,
finite magnetic-reference refinement/release, recurring informed A21, and
regional practical stability. Carry the same physical execution and estimator
state through every boundary. No new common startup deadline is assumed.

## Current result

The complete rational 4x4 LIN endpoint-action matrix is now certified for the
bound default profile. The explicit comparison is
`P_end >= E_LIN A^-1 E_LIN'`. Its singular full-state form is intentional.
`lin_matrix_certificate.py` preserves all within-LIN factors; raw neutral
S-information normalization is at least 3.59e-11 over the event-time boxes.
This number is not a full-state rho.

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

1. Bound full loss and its cross blocks, or an equivalent joint full-state
   path-action comparison. Include stable a_w/BA directions.
2. Derive recurring AG/BA factors through the literal correction/reset cadence.
   A smallest LDL pivot is not an eigenvalue floor; isolated marginal floors
   cannot be added into a block-diagonal full-state lower bound.
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
