# OU validation evidence provenance

These numerical rows come from a full simulator regeneration on the pinned
v1.2.1 28 ft vessel-response dataset. The immutable `replay_provenance` block
in `ou_validation_manifest.json` identifies the actual replay commit,
implementation/build dependencies, versioned inputs, and SHA-256 of the
normalized raw replay CSV. The manifest also records the compiler and Eigen
used to produce the bundle. Consult those fields for the revision of this
particular rerun; a later documentation commit does not change its source.

Statistical rows retain every replay and intentionally omit quality-gate status
and simulator return codes. Passing an evidence-consistency check does not mean
that the executable regression gates pass. Default-record failures are recorded
separately in `../rao_parameter_tuning/continuation/` and in CI simulation logs.

A statistical `--restat-from` records separate `restatement` provenance. It may
change derived statistics or presentation after verifying the implementation
closure and raw-row identity; it cannot replace the replay commit or replay
hashes. A replay-producing source change requires full regeneration.

`tools/ou_evidence_contract.py --check` verifies the replay hashes in a Git
checkout or a source archive without `.git`. `docs/ou-validation.md` describes
the policy. Historical bundles remain available through Git history.
