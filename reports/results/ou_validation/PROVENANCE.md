# OU validation evidence provenance

The numerical replay rows in this directory were regenerated from source commit
`538a6337fae3542b4249d4d250905c0630fdf2d9` by GitHub Actions run
`31953533100` and committed to `main` as
`bfe25ed2972d67ed806301129cc2a6547d4090ac`.

That full replay used the isotropic OU-III implementation. The manifest's
`replay_provenance` block is immutable for these rows and pins the historical
replay commit, replay-producing implementation/build closure, versioned inputs,
and SHA-256 of the normalized raw replay CSV. This full replay already uses the normalized statistical schema without
`quality_gate_pass` or `simulator_return_code`; the schema-v2 provenance migration
therefore changes manifest structure only and does not rewrite replay rows.

A statistical `--restat-from` records separate `restatement` provenance. It may
change derived statistics or presentation when analysis code changes, but it
first verifies the replay dependency closure and raw-row identity and never
changes the replay commit or replay hashes. Any replay-producing source change
requires a full simulator regeneration.

`tools/ou_evidence_contract.py --check` verifies the same replay hashes in a Git
checkout or in a source/release archive without `.git`. See
`ou_validation_manifest.json` for the machine-readable provenance and
`docs/ou-validation.md` for the policy.
