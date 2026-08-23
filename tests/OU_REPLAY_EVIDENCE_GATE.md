# OU full-study replay gate

The OU validation/robustness evidence gate is intentionally conservative. Its purpose is to make a missed required replay much less likely, even at the cost of occasionally rerunning the full study when the change would not actually affect the numerical result.

## Replay-input fingerprint

`tools/ou_replay_fingerprint.py` computes a replay-input fingerprint from:

- **every tracked file anywhere under `tests/`, regardless of filename, extension, or executable bit**;
- tracked source, script, workflow, and build files elsewhere in the repository, using a deliberately broad set of executable/source/build extensions and standard build-system filenames;
- every tracked file whose Git mode is executable, regardless of its name or extension; and
- the SHA-256 of the exact `sim-data-files.zip` archive used by the study.

The unconditional `tests/**` rule is important. Test and study configuration can be hidden in files such as `.txt`, `.csv`, `.json`, `.dat`, fixtures, reference data, or extensionless files. The gate therefore does not try to infer which test files are executable or reachable. A content change to any tracked file under `tests/` invalidates the replay-input fingerprint.

This Markdown file is itself under `tests/`, so editing it also invalidates the replay-input fingerprint. That is an accepted false positive and is preferable to creating exceptions inside the directory that carries the study protocol and parameters.

## Results-tree fingerprint

The same tool also computes a second SHA-256 fingerprint over **every file anywhere under `reports/results/`**, including each relative path and file content. Adding, deleting, renaming, or modifying any result file changes this fingerprint, regardless of filename or extension.

This gives CI a direct indicator that the committed evidence tree changed independently of the executable/test inputs. A results-tree mismatch is treated conservatively as stale evidence and requires regeneration.

The fingerprint record itself is stored at `reports/ou_evidence_fingerprint.json`, outside `reports/results/`. This is deliberate: it lets the results-tree fingerprint cover every file under `reports/results/` without a self-reference exception.

The record contains both hashes:

- `replay_fingerprint`: executable/test/workflow/build content plus the simulation ZIP;
- `results_fingerprint`: the complete `reports/results/` tree.

## Why the Git commit SHA is not the replay key

The enclosing Git commit SHA is not the replay key. Generated evidence can be committed after a successful replay, producing a new commit SHA without changing the numerical computation. Keying directly on the commit SHA would therefore make the evidence commit invalidate itself.

Generated evidence and ordinary prose outside `tests/` do not change the replay-input hash unless they otherwise match one of the broad source/build/workflow/executable rules. Changes under `reports/results/` do change the separate results-tree hash.

## CI behavior

For a full OU publication/evidence run:

1. CI checks out the target source and downloads the versioned simulation archive.
2. It computes the current replay-input fingerprint and complete results-tree fingerprint.
3. If the committed fingerprint record is missing, or either fingerprint differs, CI runs the complete validation and robustness simulator studies, regenerates the evidence, and writes both new fingerprints.
4. If both fingerprints are unchanged, CI may reuse the committed full-study rows, but it still runs the validation/evidence contract tests.
5. After regeneration, the generated validation bundle, robustness bundle, mirrored manuscript inputs, and `reports/ou_evidence_fingerprint.json` are committed together.
6. If a write retry rebases onto a newer branch tip, both fingerprints are checked again before regenerated evidence can be pushed.

The policy is intentionally asymmetric:

> An unnecessary full replay is acceptable. Reusing stale evidence after a potentially relevant repository, test-data, simulation-data, or results-tree change is not.
