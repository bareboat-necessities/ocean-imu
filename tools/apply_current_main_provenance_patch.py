#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace(path: Path, old: str, new: str, count: int = 1) -> None:
    text = path.read_text(encoding="utf-8")
    if text.count(old) < count:
        raise SystemExit(f"anchor not found in {path}: {old[:180]!r}")
    path.write_text(text.replace(old, new, count), encoding="utf-8")

prov = ROOT / "tools/ou_evidence_provenance.py"
text = prov.read_text(encoding="utf-8")
old = '''    old_analysis = {\n        name: normalize_record(record)\n        for name, record in current.get("analysis_pipeline_files", {}).items()\n    }\n    replay_environment = {\n'''
new = '''    old_analysis = {\n        name: normalize_record(record)\n        for name, record in current.get("analysis_pipeline_files", {}).items()\n    }\n    current_command = [str(item) for item in current.get("command", [])]\n    current_is_restatement = bool(\n        current.get("restated_from")\n        or current.get("protocol", {}).get("restated_from")\n        or "--restat-from" in current_command\n    )\n    historical_fields = csv.DictReader(\n        io.StringIO(historical_raw.decode("utf-8"), newline="")\n    ).fieldnames or []\n    removed_gate_columns = [name for name in GATE_FIELDS if name in historical_fields]\n    replay_environment = {\n'''
if old not in text:
    raise SystemExit("migration analysis anchor missing")
text = text.replace(old, new, 1)
text = text.replace(
    '''            "normalization_only_columns_removed": list(GATE_FIELDS),\n''',
    '''            "normalization_only_columns_removed": removed_gate_columns,\n''',
    1,
)
old = '''    if old_analysis:\n        current["restatement"] = {\n            "git_commit": current.get("git_commit"),\n            "analysis_pipeline_files": old_analysis,\n            "source_bundle_sha256": None,\n            "source_manifest_sha256": None,\n            "restated_at": None,\n            "environment": {\n                "python": current.get("python", "unrecorded"),\n                "numpy": current.get("numpy", "unrecorded"),\n                "platform": current.get("platform", "unrecorded"),\n            },\n            "provenance_migration_note": (\n                "restatement metadata recovered from the legacy manifest; "\n                "source bundle/manifest hashes and timestamp were not recorded "\n                "and are intentionally not fabricated"\n            ),\n        }\n    _replay_aliases(current)\n'''
new = '''    if old_analysis and current_is_restatement:\n        current["restatement"] = {\n            "git_commit": current.get("git_commit"),\n            "analysis_pipeline_files": old_analysis,\n            "source_bundle_sha256": None,\n            "source_manifest_sha256": None,\n            "restated_at": None,\n            "environment": {\n                "python": current.get("python", "unrecorded"),\n                "numpy": current.get("numpy", "unrecorded"),\n                "platform": current.get("platform", "unrecorded"),\n            },\n            "provenance_migration_note": (\n                "restatement metadata recovered from the legacy manifest; "\n                "source bundle/manifest hashes and timestamp were not recorded "\n                "and are intentionally not fabricated"\n            ),\n        }\n    else:\n        current.pop("restatement", None)\n    _replay_aliases(current)\n'''
if old not in text:
    raise SystemExit("migration restatement anchor missing")
prov.write_text(text.replace(old, new, 1), encoding="utf-8")

# Committed robustness evidence may be either a fresh full generation (analysis
# hashes live at top level for backward compatibility) or a later restatement
# (analysis hashes live in the restatement block). Tests accept both without
# weakening replay dependency checks.
robust_test = ROOT / "tests/validation/test_ou_robustness.py"
rtext = robust_test.read_text(encoding="utf-8")
old = '''        analysis = set(manifest["restatement"]["analysis_pipeline_files"])\n'''
new = '''        analysis = set(\n            manifest.get("restatement", {}).get(\n                "analysis_pipeline_files", manifest.get("analysis_pipeline_files", {})\n            )\n        )\n'''
if rtext.count(old) != 1:
    raise SystemExit(f"expected one set(restatement analysis) anchor, found {rtext.count(old)}")
rtext = rtext.replace(old, new, 1)
old = '''        analysis = manifest["restatement"]["analysis_pipeline_files"]\n'''
new = '''        analysis = manifest.get("restatement", {}).get(\n            "analysis_pipeline_files", manifest.get("analysis_pipeline_files", {})\n        )\n'''
if rtext.count(old) != 1:
    raise SystemExit(f"expected one restatement analysis mapping anchor, found {rtext.count(old)}")
rtext = rtext.replace(old, new, 1)
robust_test.write_text(rtext, encoding="utf-8")

# Update current-evidence lineage in reproducibility documentation. This is a
# real full replay produced from the merged article commit, not a replay run by
# this feature branch.
for relative in ("docs/ou-validation.md", "reports/results/ou_validation/PROVENANCE.md"):
    path = ROOT / relative
    value = path.read_text(encoding="utf-8")
    value = value.replace("31943137398", "31953533100")
    value = value.replace("5df1f3b42d3eb456961d12e598257b5859470451", "538a6337fae3542b4249d4d250905c0630fdf2d9")
    value = value.replace("35120f96d3f8c02872fd3e06fa94ebe547c3f5eb", "bfe25ed2972d67ed806301129cc2a6547d4090ac")
    value = value.replace(
        "The later removal of the legacy\n`quality_gate_pass` and `simulator_return_code` columns is recorded as a verified\nschema migration, not as a new simulator generation.\n",
        "This full replay already uses the normalized statistical schema without\n`quality_gate_pass` or `simulator_return_code`; the schema-v2 provenance migration\ntherefore changes manifest structure only and does not rewrite replay rows.\n",
    )
    value = value.replace(
        "that its replay-producing dependency closure had not\nchanged since the recorded full replay and that the current normalized raw CSV\nwas exactly the historical replay CSV with only the retired regression-gate\ncolumns removed. No numerical replay metric or paired statistic was changed by\nthat provenance migration.",
        "that its replay-producing dependency closure had not\nchanged since the recorded full replay and that the current normalized raw CSV\nwas byte-identical to the archived full-replay CSV. No numerical replay metric\nor paired statistic was changed by that provenance migration.",
    )
    path.write_text(value, encoding="utf-8")

print("adjusted provenance migration for current-main full evidence")
