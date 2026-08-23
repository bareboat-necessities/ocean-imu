# OU full-study replay gate

The OU validation/robustness evidence gate is intentionally conservative. Its purpose is to make a missed required replay much less likely, even at the cost of occasionally rerunning the full study when the change would not actually affect the numerical result.

## What is fingerprinted

`tools/ou_replay_fingerprint.py` computes a content fingerprint from:

- **every tracked file anywhere under `tests/`, regardless of filename, extension, or executable bit**;
- tracked source, script, workflow, and build files elsewhere in the repository, using a deliberately broad set of executable/source/build extensions and standard build-system filenames;
- every tracked file whose Git mode is executable, regardless of its name or extension; and
- the SHA-256 of the exact `sim-data-files.zip` archive used by the study.

The unconditional `tests/**` rule is important. Test and study configuration can be hidden in files such as `.txt`, `.csv`, `.json`, `.dat`, fixtures, reference data, or extensionless files. The gate therefore does not try to infer which test files are executable or reachable. A content change to any tracked file under `tests/` invalidates the replay fingerprint.

This Markdown file is itself under `tests/`, so editing it also invalidates the fingerprint. That is an accepted false positive and is preferable to creating exceptions inside the directory that carries the study protocol and parameters.

## What is not used as the key

The enclosing Git commit SHA is not the replay key. Generated evidence can be committed after a successful replay, producing a new commit SHA without changing the numerical computation. Keying directly on the commit SHA would therefore make the evidence commit invalidate itself.

Generated evidence and ordinary prose outside `tests/` are not automatically included unless they otherwise match one of the broad source/build/workflow/executable rules.

## CI behavior

For a full OU publication/evidence run:

1. CI checks out the target source and downloads the versioned simulation archive.
2. It computes the current replay fingerprint.
3. If the committed fingerprint record is missing or differs, CI runs the complete validation and robustness simulator studies, regenerates the evidence, and writes the new fingerprint alongside that evidence.
4. If the fingerprint is unchanged, CI may reuse the committed full-study rows, but it still runs the validation/evidence contract tests.
5. If a write retry rebases onto a newer branch tip, the fingerprint is checked again before regenerated evidence can be pushed.

The policy is intentionally asymmetric:

> An unnecessary full replay is acceptable. Reusing stale evidence after a potentially relevant repository or test-data change is not.

The committed fingerprint record is `reports/results/ou_replay_fingerprint.json`.