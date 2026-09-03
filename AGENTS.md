# AGENTS.md

## Repository expectations

- Keep changes minimal and targeted.
- Do not change public APIs, filenames, or build targets unless the task requires it.

## Build and validation

- Build requires data fetched from release of https://github.com/bareboat-necessities/oceanography-waves-lib (Check .github/workflows/build.yml for instructions how to fetch simulation data file sim-data-files.zip)
- Primary validation command: `make all`
- After changing any `.c`, `.cc`, `.cpp`, `.h`, `.hpp`, `.mk`, or `Makefile`, run `make all`.
- If the build fails, report the exact failing command and error.

## Eigen handling

- This repo may run in environments where `Eigen/Dense` is not on the default compiler include path.
- Prefer an existing vendored Eigen copy if the repo already contains one, such as `third_party/eigen`.
- Otherwise preserve or use an include path equivalent to:
    - `CPPFLAGS += -I/usr/include/eigen3`
    - or `make all CPPFLAGS+='-I/usr/include/eigen3'`
- Do not remove Eigen usage just to make the build pass.
- If the environment itself is missing Eigen, note that Codex cloud setup must install `libeigen3-dev`.

## Makefile policy

- Prefer small Makefile fixes that keep both local builds and Codex builds working.
- Prefer a fallback pattern like:
    - use vendored Eigen if present
    - otherwise use `/usr/include/eigen3`

## Change review checklist

Before finishing:
- run `make all`
- keep diffs small
- avoid unrelated formatting churn
- summarize exactly what changed and why

## Making PR

While preparing PRs make sure you do not leave dangling and obsolete
tests, .md files, tools, .py files and data files. Everything that PR
makes obsolete should be clean up before making PR ready, obsolete .yml GitHub workflows too.
Do not keep all historical context in .md files. It creates unnecessary bloat.
Documentation should state current behavior. It doesn't need to explain what
PR changed. What PR changed belongs to PR meta, but not to the committed files.

## OU-III proof research protocol

The objective is to establish or falsify the declared OU-III theorem stages. Existing proof constructions are hypotheses, not requirements. Do not preserve a proof method merely because earlier commits used it.

The mandatory research state machine is:

`EXECUTE -> FAILURE ANALYSIS -> REPLAN -> EXECUTE`

Never go directly from a failed mathematical attempt to another refinement of the same idea.

### Failure analysis is mandatory

After every mathematical, enclosure, conditioning, or numerical failure, stop implementation and record in `docs/ou3-proof-research-state.md`:

1. the exact failed inequality, quantity, source class, state direction, or numerical operation;
2. whether the failure is a theorem failure, proof-method failure, enclosure failure, conditioning failure, implementation defect, CI failure, or infrastructure failure;
3. which hypothesis or strategy the evidence invalidates;
4. what the failure does not invalidate;
5. the current limiting quantity;
6. the next falsifiable experiment.

If essentially the same mechanism has already failed once, list at least three qualitatively different alternatives before choosing the next method. The alternatives must not all be variants of subdivision, tighter scalar norms, or increased search depth.

### Two-strike rule and quantitative permission

A proof tactic gets one initial implementation and at most one mathematically motivated refinement. After the second failure of substantially the same mechanism, freeze it as a dead end and perform an architecture review.

Do not repeat subdivision, interval refinement, tighter scalar norm bounds, additional remainder lemmas, or increased search depth unless a scaling argument predicts quantitatively that the change can cross the required theorem threshold. "Try it and see" is not sufficient.

A rejected route may be revisited only after documenting the new mathematical fact that makes the previous rejection inapplicable.

### Master inequality before implementation

Before implementing a new proof lemma, identify where it enters the controlling theorem inequality and estimate whether improving it can materially change the final margin.

For P4, the controlling object is the complete-word ratio

`rho_w = V_after(F_w(x)) / V_before(x)`.

Local reset, Joseph, sector, measurement, or remainder bounds are subordinate lemmas. Do not create a new P4 micro-certificate unless it is required by a complete-word formulation already shown to be feasible and its expected contribution can materially move `rho_w` through 1.

### Distinguish a false theorem from a bad enclosure

Before rigorous certification of a new P4 formulation, run a non-promoting, high-precision complete-word feasibility diagnostic that preserves the exact matrix/group structure and introduces essentially no interval pessimism.

Interpret the result as follows:

- diagnostic `rho` clearly below 1: rigorous enclosure work is justified;
- diagnostic `rho` near 1: identify the limiting state direction and reconsider the theorem/metric before interval implementation;
- diagnostic `rho` above 1: abandon that P4 formulation rather than sharpening unrelated bounds.

The diagnostic must report worst H/A ratios, the limiting legal word/source/phase, a maximizing state direction or generalized eigenvector, operation-by-operation margin consumption, and distance to `rho=1`.

### Independent critic pass

After every important proof failure, perform a critic pass before further implementation:

- assume the current proof architecture is wrong;
- state the strongest reason to abandon it;
- propose an alternative that does not share the failed mechanism;
- compare the alternatives against the actual limiting quantity.

Only after this pass may implementation resume.

### Research ledger

Maintain `docs/ou3-proof-research-state.md` as a short current-state ledger with:

- Current hypothesis
- Evidence
- Current limiter
- Failed approaches / DEAD_ENDS
- Retained facts
- Alternatives
- Next falsifiable experiment

Keep it concise and current. It is research state, not a chronological PR history.

### One research question per PR

Do not let an OU-III proof PR absorb successive unrelated theorem stages. A feasibility result that justifies a new rigorous proof stage starts a new PR unless the current PR explicitly defines that stage as its research question.

### Immutable OU-III constraints

Do not move these constraints for proof convenience:

- canonical P3 usefulness threshold = `1e-18`;
- no replay fitting;
- no artificial operating-domain reduction to obtain PASS;
- no deployed-filter change to make the proof easier;
- same-history source correlation must be preserved;
- both H=18 and A=21 are required;
- lever arm remains zero/disabled in the current proof scope;
- dormant vibration-guard branch remains transparent;
- P4 cannot promote before canonical P3;
- P5 cannot promote before strict canonical P4 contraction.

A failed proof attempt is evidence. Prefer abandoning a bad proof architecture over accumulating increasingly elaborate bounds around it. The goal is not to make the current algorithm pass; the goal is to determine the mathematically correct route to the theorem, or determine that the current theorem formulation cannot be certified.
