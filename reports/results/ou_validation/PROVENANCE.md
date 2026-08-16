# OU validation evidence provenance

The numerical replay rows in this directory were regenerated from source commit
`5df1f3b42d3eb456961d12e598257b5859470451` by GitHub Actions run
`31943137398` and committed to `main` as
`35120f96d3f8c02872fd3e06fa94ebe547c3f5eb`.

That regeneration used the current isotropic OU-III implementation. The
subsequent evidence-schema normalization removes only the historical simulator
pass/fail fields and adds implementation/analysis provenance hashes; it does not
change the continuous replay metrics, summary rows, paired effects, or figures.

See `ou_validation_manifest.json` for the machine-readable protocol, source
commit, implementation dependency hashes, analysis-pipeline hashes, input
hashes, and result hashes.
