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

## OU-III Proof Research Protocol

This applies to `tools/ou3_*.py`, `tests/validation/test_ou3_*.py`,
`doc/kalman_ou_iii/**` and `.github/workflows/ou3-proof.yml`.

Read `docs/ou3-proof-research-state.md` before modifying any proof code.

The primary objective is to establish or falsify the declared OU-III theorem
stages. Existing proof constructions are hypotheses, not requirements. Do not
preserve a proof method merely because previous commits used it.

After every mathematical or numerical failure, STOP implementation and perform
a failure analysis before making another change. Record:

1. the exact failed inequality, quantity, source class, state direction, or
   numerical operation;
2. whether the failure is a theorem failure, proof-method failure, enclosure
   failure, conditioning failure, implementation defect, or infrastructure
   failure;
3. what previous assumption or strategy the failure invalidates;
4. what the failure does NOT invalidate;
5. at least two qualitatively different alternative approaches when the same
   mechanism has already failed once;
6. a quantitative reason why the selected next approach can improve the
   limiting quantity.

Do not repeat the same remedy more than once. Repeated subdivision, interval
refinement, tighter scalar norm bounds, additional remainder lemmas, and
increased search depth are prohibited after a second failure unless a scaling
argument demonstrates that the refinement can cross the required threshold.

Two-strike rule: if essentially the same mechanism fails twice it is frozen as
a dead end. Do not attempt a third refinement until at least three
qualitatively different alternatives have been stated and compared.

The workflow states are EXECUTE, FAILURE ANALYSIS, REPLAN. Never go from
EXECUTE straight back to EXECUTE on the same idea. After an important failure,
take an adversarial pass first: assume the current proof architecture is wrong,
state the strongest reason to abandon it, and find an alternative that does not
share the failed mechanism.

Iteration budget per idea: one initial implementation, at most one
mathematically motivated refinement, then a mandatory architecture review.

Before implementing a new proof lemma, identify exactly where it enters the
master theorem inequality and estimate whether tightening it can materially
affect the final margin. A lemma worth 1e-5 against a failure of +0.4 is
irrelevant and must not be implemented.

For P4 the controlling object is the complete-word ratio

    rho_w = V_after(F_w(x)) / V_before(x).

Local bounds and certificates are subordinate lemmas only. Do not create new P4
micro-certificates unless they are required by a complete-word formulation
already shown to be feasible.

Before rigorous certification of a new P4 formulation, run a non-promoting
high-precision feasibility diagnostic preserving the exact matrix/group
structure, and interpret it as:

* diagnostic rho clearly below 1 -- rigorous enclosure work is justified;
* diagnostic rho near 1 -- identify the limiting state direction and reconsider
  the theorem/metric before interval implementation;
* diagnostic rho above 1 -- abandon that P4 formulation. Do not attempt to
  rescue it by unrelated bound sharpening.

Always distinguish "the theorem may be false" from "our bound is bad".

Immutable constraints:

* canonical P3 usefulness threshold = 1e-18;
* no replay fitting;
* no artificial operating-domain reduction to obtain PASS;
* no deployed-filter change for proof convenience;
* same-history source correlation must be preserved;
* H=18 and A=21 are both required;
* lever arm remains zero/disabled in the current proof scope;
* dormant vibration-guard branch remains transparent;
* P4 cannot promote before canonical P3;
* P5 cannot promote before strict canonical P4 contraction.

One research question per PR. A feasibility result starts a new PR; it does not
grow the current one.

A failed proof attempt is useful information. Prefer abandoning a bad proof
architecture over accumulating increasingly elaborate bounds around it. The
goal is not to make the current algorithm pass. The goal is to determine the
mathematically correct route to the theorem, or to determine that the current
theorem formulation cannot be certified.

After a failure, permission to continue the same idea must be earned by
predicting quantitatively why the next modification will fix it. Otherwise
replan.

## Making PR

While preparing PRs make sure you do not leave dangling and obsolete
tests, .md files, tools, .py files and data files. Everything that PR
makes obsolete should be clean up before making PR ready, obsolete .yml GitHub workflows too. 
Do not keep all historical context in .md files. It creates unnecessary bloat.
Documentation should state current behavior. It doesn't need to explain what
PR changed. What PR changed belongs to PR meta, but not to the committed files.
