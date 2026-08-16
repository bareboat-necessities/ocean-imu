# OU validation evidence provenance

The numerical replay rows in this directory were regenerated from source commit
`5df1f3b42d3eb456961d12e598257b5859470451` by GitHub Actions run
`31943137398` and committed to `main` as
`35120f96d3f8c02872fd3e06fa94ebe547c3f5eb`.

That full replay used the isotropic OU-III implementation. The manifest's
`replay_provenance` block is immutable for these rows and pins the historical
replay commit, replay-producing implementation/build closure, versioned inputs,
and SHA-256 of the normalized raw replay CSV. The later removal of the legacy
`quality_gate_pass` and `simulator_return_code` columns is recorded as a verified
schema migration, not as a new simulator generation.

A statistical `--restat-from` records separate `restatement` provenance. It may
change derived statistics or presentation when analysis code changes, but it
first verifies the replay dependency closure and raw-row identity and never
changes the replay commit or replay hashes. Any replay-producing source change
requires a full simulator regeneration.

`tools/ou_evidence_contract.py --check` verifies the same replay hashes in a Git
checkout or in a source/release archive without `.git`. See
`ou_validation_manifest.json` for the machine-readable provenance and
`docs/ou-validation.md` for the policy.
