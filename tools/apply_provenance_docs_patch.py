#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace(path: Path, old: str, new: str, count: int = 1) -> None:
    text = path.read_text(encoding="utf-8")
    if text.count(old) < count:
        raise SystemExit(f"anchor not found in {path}: {old[:180]!r}")
    path.write_text(text.replace(old, new, count), encoding="utf-8")

readme = ROOT / "README.md"
replace(
    readme,
    "- OU_III is using higher-order integral drift correction. It is 21 dimensional state Kalman filter. It seems works better on high precision IMUs.\n",
    "- OU_III uses higher-order integral drift correction. It is a 21-dimensional-state Kalman filter and is the OU variant used for the main 3-D navigation study.\n",
)
replace(
    readme,
    "- OU_II is using more direct integral drift correction and more responsive to sea state changes. It is 18 dimentional state Kalman filter.\n",
    "- OU_II uses more direct integral drift correction and is more responsive to sea-state changes. It is an 18-dimensional-state Kalman filter.\n",
)
replace(
    readme,
    """The OU-II/OU-III statistical validation runner keeps the executable regression\ngates, and scores both them and the primary full experiment over the final\n15 minutes of each 20-minute run. It pairs wave realization, IMU noise, and\ninitialization across both filters and the Adaptive, FixedNominal, and\nFixedOracle ablations.\n""",
    """The OU-II/OU-III statistical validation runner scores continuous metrics over\nthe final 15 minutes of each 20-minute run. Deterministic simulator regression\ngates remain executable diagnostics, but they are not Monte Carlo inclusion\ncriteria and are not exported as pass/fail evidence. Wave realization, IMU\nnoise, and initialization are paired across both filters and matched ablations.\n""",
)
replace(
    readme,
    """Full mode produces raw and summary CSV, JSON, LaTeX, paired-effect, manifest,\nand SVG plot artifacts under `reports/results/ou_validation/`. The completed\nten-seed, 300-run evidence used by the article is versioned in that directory;\nit shows a consistent stationary vertical OU-III benefit but mixed 3D,\ntransition, and adaptation outcomes rather than uniform superiority. See\n[`docs/ou-validation.md`](docs/ou-validation.md) for the protocol, seed controls,\nand interpretation.\n""",
    """Full mode produces raw and summary CSV, JSON, LaTeX, paired-effect, manifest,\nand SVG plot artifacts under `reports/results/ou_validation/`. The versioned\nten-seed study contains **840 simulator replays** across the declared stationary,\ntransition, comparison, covariance-policy, and channel-ablation configurations;\npaired seed-level aggregates, not 840 independent samples, are used for the\nprimary inference. See [`docs/ou-validation.md`](docs/ou-validation.md) for the\nprotocol, provenance contract, seed controls, and interpretation.\n""",
)

validation = ROOT / "docs/ou-validation.md"
replace(
    validation,
    """The current full bundle was regenerated from the current isotropic OU-III\nimplementation by GitHub Actions run `31943137398`. The replay source commit is\n`5df1f3b42d3eb456961d12e598257b5859470451`; the generated evidence was committed\nas `35120f96d3f8c02872fd3e06fa94ebe547c3f5eb`. The manifest records that source\ncommit and, through the evidence contract, hashes the transitive repository-local\nC/C++ implementation dependency closure of the OU simulator targets together\nwith the analysis pipeline. A filter/header change therefore makes the committed\nevidence stale until a full regeneration is performed.\n""",
    """The current full bundle was regenerated from the current isotropic OU-III\nimplementation by GitHub Actions run `31943137398`. Its immutable replay source\ncommit is `5df1f3b42d3eb456961d12e598257b5859470451`; the generated evidence was\ncommitted as `35120f96d3f8c02872fd3e06fa94ebe547c3f5eb`. The manifest now separates\nthat replay provenance from any later statistical restatement. Replay provenance\npins the simulator/filter dependency closure, build files, versioned inputs, and\nnormalized raw-row hash; a later restatement records its analysis context in a\nseparate block and cannot replace those replay pins.\n""",
)
start = validation.read_text(encoding="utf-8")
old = """## Provenance contract\n\n`ou_validation_manifest.json` records the command, versions, source commit,\nprotocol, seed triplets, input hashes, result hashes, analysis-pipeline hashes,\nand the transitive implementation hashes reachable from the simulator sources.\nThe implementation closure includes the simulator translation units and their\nrepository-local headers, including the actual OU-II/OU-III filter headers and\nshared simulator code.\n\n`tools/ou_evidence_contract.py --check` verifies that the recorded\nimplementation and analysis hashes still match the checkout. The normal\n`tests/validation` target runs this contract automatically. A source change in\nthe estimator path therefore fails the evidence check instead of silently\nleaving publication numbers attached to a different implementation.\n\nDuring a full regeneration, the same contract is allowed to normalize and stamp\nthe newly generated bundle only when the manifest's recorded replay commit is\nthe checked-out source commit. It cannot restamp old rows after implementation\ncode changes.\n"""
new = """## Provenance contract\n\n`ou_validation_manifest.json` has two provenance layers. `replay_provenance` is\nimmutable for a set of simulator rows: it records the replay source commit, the\ntransitive repository-local simulator/filter dependency closure, the simulator\nbuild Makefiles, input hashes, and a SHA-256 of the normalized raw replay CSV.\nThe build environment is also recorded when a new full replay is generated\n(compiler identity/version, Eigen identity, Python, NumPy, and platform), but\nthose environment fields are informational rather than cross-platform validity\nkeys.\n\nA later `--restat-from` may recompute summaries, paired effects, tables, or\npublication text from those same rows. Before writing anything it compares the\ncurrent replay dependency closure with the immutable replay closure and verifies\nthe raw-row hash. If a filter header, simulator translation unit, shared replay\ndependency, or recorded build file changed, restatement fails and requires a\nfull simulator replay. If only approved analysis/presentation code changed, the\nrestatement records a separate `restatement` block with its analysis hashes and\nsource-bundle hash; it never changes the replay commit or replay dependency\nhashes.\n\n`tools/ou_evidence_contract.py --check` enforces the replay pins and result\ninventory. In a Git checkout it can additionally validate the recorded replay\ncommit; in a GitHub source ZIP, release archive, Zenodo artifact, or copied tree\nwithout `.git`, it falls back to read-only manifest/file-hash verification and\nretains the same dependency-hash checks. `--auto` cannot promote a legacy\nrestatement or old rows to a new estimator revision.\n\nThe existing isotropic bundle was migrated to this schema without replay only\nafter Git history proved that its replay-producing dependency closure had not\nchanged since the recorded full replay and that the current normalized raw CSV\nwas exactly the historical replay CSV with only the retired regression-gate\ncolumns removed. No numerical replay metric or paired statistic was changed by\nthat provenance migration.\n"""
if old not in start:
    raise SystemExit("validation provenance section anchor missing")
validation.write_text(start.replace(old, new, 1), encoding="utf-8")

robustness = ROOT / "docs/ou-robustness.md"
text = robustness.read_text(encoding="utf-8")
old = """## Provenance\n\nThe robustness manifest now records two distinct source classes in addition to\ninput/result hashes:\n\n- `implementation_files`: the transitive repository-local C/C++ dependency\n  closure reachable from the OU-III simulator and shared simulation code;\n- `analysis_pipeline_files`: the validation/robustness Python analysis sources.\n\nThis is stronger than pinning only the simulator translation unit. A change in\n`SeaStateFusionFilter_OU_III.h`, `Kalman3D_Wave_OU_III.h`, or another local\nheader in the simulator include graph invalidates the committed evidence.\n\n`tools/ou_evidence_contract.py --check` verifies those hashes and the clean\nstatistical schema. The normal `tests/validation` target runs the same contract.\nThus a code change cannot leave the old Monte Carlo bundle looking current just\nbecause the wave-input hashes still match.\n\nThe protocol also records that all completed replays with machine-readable\nmetrics are included and that simulator regression gates are not exported. This\nkeeps regression testing and inferential evidence separate.\n"""
new = """## Provenance\n\nThe robustness manifest uses the same two-layer contract as the primary study.\nImmutable `replay_provenance` pins the full OU-III simulator/filter dependency\nclosure, the simulator build file, replay inputs, and the normalized raw-row\nSHA-256. A later `restatement` block may identify different Python analysis code\nand derived outputs, but it cannot replace the replay source commit or replay\nhashes.\n\nBefore `--restat-from` writes anything, the current replay closure is compared\nwith the recorded one. A change in `SeaStateFusionFilter_OU_III.h`,\n`Kalman3D_Wave_OU_III.h`, the simulator, shared replay code, or the recorded\nMakefile is therefore a hard replay requirement rather than a warning. Analysis\nor presentation changes may be restated only while the immutable replay and raw\nrows still verify.\n\n`tools/ou_evidence_contract.py --check` works both in a Git checkout and in a\nsource archive without `.git`; archive mode uses the committed manifest and\nfile hashes and does not weaken dependency verification. The protocol separately\nrecords that all completed machine-readable replays are included and that\ndeterministic simulator regression gates are not Monte Carlo acceptance\ncriteria.\n"""
if old not in text:
    raise SystemExit("robustness provenance section anchor missing")
robustness.write_text(text.replace(old, new, 1), encoding="utf-8")

method = ROOT / "doc/kalman_ou_iii/w3d-sim-charts.tex-part"
replace(
    method,
    """Descriptive secondary comparisons retain paired intervals\nwithout promoting each one to an additional confirmatory hypothesis.  The full\nprotocol, seed triplets, implementation-source hashes, and generated results are\nversioned in the committed evidence bundle.\n""",
    """Descriptive secondary comparisons retain paired intervals\nwithout promoting each one to an additional confirmatory hypothesis.  A full\nsimulator generation establishes immutable replay provenance: replay source\ncommit, replay-producing implementation/build hashes, input hashes, and the\nnormalized raw-row hash.  A later statistical restatement is recorded\nseparately and is permitted only while that replay dependency closure and the\nraw rows still match; it cannot substitute for a replay after implementation\nchanges.  The evidence contract enforces the same file-hash guarantees in a\nsource archive even when Git metadata is unavailable.\n""",
)

prov = ROOT / "reports/results/ou_validation/PROVENANCE.md"
prov.write_text("""# OU validation evidence provenance\n\nThe numerical replay rows in this directory were regenerated from source commit\n`5df1f3b42d3eb456961d12e598257b5859470451` by GitHub Actions run\n`31943137398` and committed to `main` as\n`35120f96d3f8c02872fd3e06fa94ebe547c3f5eb`.\n\nThat full replay used the isotropic OU-III implementation. The manifest's\n`replay_provenance` block is immutable for these rows and pins the historical\nreplay commit, replay-producing implementation/build closure, versioned inputs,\nand SHA-256 of the normalized raw replay CSV. The later removal of the legacy\n`quality_gate_pass` and `simulator_return_code` columns is recorded as a verified\nschema migration, not as a new simulator generation.\n\nA statistical `--restat-from` records separate `restatement` provenance. It may\nchange derived statistics or presentation when analysis code changes, but it\nfirst verifies the replay dependency closure and raw-row identity and never\nchanges the replay commit or replay hashes. Any replay-producing source change\nrequires a full simulator regeneration.\n\n`tools/ou_evidence_contract.py --check` verifies the same replay hashes in a Git\ncheckout or in a source/release archive without `.git`. See\n`ou_validation_manifest.json` for the machine-readable provenance and\n`docs/ou-validation.md` for the policy.\n""", encoding="utf-8")

workflow = ROOT / ".github/workflows/ou-validation.yml"
text = workflow.read_text(encoding="utf-8")
for anchor in ('      - "tools/ou_validation.py"\n',):
    occurrences = text.count(anchor)
    if occurrences != 2:
        raise SystemExit(f"expected two workflow path anchors, found {occurrences}")
text = text.replace(
    '      - "tools/ou_validation.py"\n',
    '      - "tools/ou_validation.py"\n      - "tools/ou_evidence_provenance.py"\n      - "tools/ou_evidence_contract.py"\n',
)
workflow.write_text(text, encoding="utf-8")

print("updated README, provenance docs, paper methodology, and workflow paths")
