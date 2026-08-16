#!/usr/bin/env python3
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
path = TOOLS / "ou_evidence_provenance.py"
text = path.read_text(encoding="utf-8")

old = '''def _input_records_from_legacy(manifest: Mapping[str, Any]) -> dict[str, dict[str, object]]:\n    return {name: normalize_record(record) for name, record in manifest.get("source_files", {}).items()}\n'''
new = '''def _input_records_from_legacy(study: str, manifest: Mapping[str, Any]) -> dict[str, dict[str, object]]:\n    # Legacy robustness manifests mixed replay inputs, simulator sources, and\n    # Python analysis files under one source_files key. Preserve only actual\n    # replay inputs here; implementation/build and analysis provenance have\n    # their own explicit namespaces in schema v2.\n    excluded = set(implementation_records(study)) | set(analysis_records(study))\n    return {\n        name: normalize_record(record)\n        for name, record in manifest.get("source_files", {}).items()\n        if name not in excluded\n    }\n'''
if old not in text:
    raise SystemExit("legacy input-record helper anchor not found")
text = text.replace(old, new, 1)
text = text.replace('_input_records_from_legacy(manifest)', '_input_records_from_legacy(study, manifest)')
text = text.replace('_input_records_from_legacy(historical)', '_input_records_from_legacy(study, historical)')

old = '''    if git_available():\n        commit = replay.get("git_commit")\n        if not commit:\n            errors.append(f"{study}: replay source commit is missing")\n        elif subprocess.run(\n            ("git", "cat-file", "-e", f"{commit}^{{commit}}"), cwd=REPO_ROOT,\n            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,\n        ).returncode != 0:\n            errors.append(f"{study}: replay source commit is unavailable in this checkout: {commit}")\n    return errors\n'''
new = '''    commit = str(replay.get("git_commit", ""))\n    if not re.fullmatch(r"[0-9a-f]{40}", commit):\n        errors.append(f"{study}: replay source commit is missing or malformed")\n    # A shallow Git checkout may not contain the historical replay object.\n    # Hash validation above is authoritative for the committed source tree;\n    # full Git history is required only for the explicit legacy migration.\n    return errors\n'''
if old not in text:
    raise SystemExit("replay git-check anchor not found")
text = text.replace(old, new, 1)

anchor = '''def begin_restatement(study: str, source_bundle: Path) -> RestatementContext:\n'''
helper = '''def _assert_bundle_rows_match_raw_csv(source_bundle: Path, raw_csv: Path) -> None:\n    bundle = json.loads(source_bundle.read_text(encoding="utf-8"))\n    rows = bundle.get("raw_runs", [])\n    with raw_csv.open(newline="", encoding="utf-8") as handle:\n        reader = csv.DictReader(handle)\n        fields = reader.fieldnames or []\n        csv_rows = list(reader)\n    if len(rows) != len(csv_rows):\n        raise ProvenanceError("bundle raw_runs count differs from canonical raw replay CSV")\n    expected_fields = set(fields)\n    for index, (row, csv_row) in enumerate(zip(rows, csv_rows)):\n        if set(row) != expected_fields:\n            raise ProvenanceError(f"bundle raw_runs schema differs from canonical raw CSV at row {index}")\n        for field in fields:\n            if str(row[field]) != csv_row[field]:\n                raise ProvenanceError(\n                    f"bundle raw_runs differs from canonical raw replay CSV at row {index}, field {field}"\n                )\n\n\n'''
if anchor not in text:
    raise SystemExit("begin_restatement anchor missing")
text = text.replace(anchor, helper + anchor, 1)
old = '''    errors = replay_errors(study, manifest)\n    if errors:\n        raise ProvenanceError(\n'''
new = '''    errors = replay_errors(study, manifest)\n    if not source_raw.is_file():\n        errors.append(f"{study}: canonical raw replay CSV is missing: {source_raw}")\n    if errors:\n        raise ProvenanceError(\n'''
if old not in text:
    raise SystemExit("restatement replay-errors anchor missing")
text = text.replace(old, new, 1)
old = '''            + "\\n".join(errors)\n        )\n    return RestatementContext(\n'''
new = '''            + "\\n".join(errors)\n        )\n    _assert_bundle_rows_match_raw_csv(source_bundle, source_raw)\n    return RestatementContext(\n'''
if old not in text:
    raise SystemExit("restatement return anchor missing")
text = text.replace(old, new, 1)

old = '''            "source_bundle_sha256": sha256_file(out / cfg["json"]),\n            "source_manifest_sha256": None,\n            "restated_at": None,\n'''
new = '''            "source_bundle_sha256": None,\n            "source_manifest_sha256": None,\n            "restated_at": None,\n'''
if old not in text:
    raise SystemExit("legacy recovered restatement anchor missing")
text = text.replace(old, new, 1)

# Fresh replay provenance is initialized in the commit/check job, but the
# compiler and Eigen environment must describe the combine job that actually
# rebuilt the simulator. Full-generation manifests now carry that environment
# forward explicitly.
old = '''        "environment": environment_metadata(),\n    }\n\n\ndef _compare_record_maps'''
new = '''        "environment": manifest.get("build_environment") or environment_metadata(),\n    }\n\n\ndef _compare_record_maps'''
if old not in text:
    raise SystemExit("fresh replay environment anchor missing")
text = text.replace(old, new, 1)
path.write_text(text, encoding="utf-8")

validation = TOOLS / "ou_validation.py"
vtext = validation.read_text(encoding="utf-8")
old = '''        "numpy": np.__version__,\n        "command": [sys.executable, *sys.argv],\n        "protocol": protocol,\n'''
new = '''        "numpy": np.__version__,\n        "build_environment": evidence_provenance.environment_metadata(),\n        "command": [sys.executable, *sys.argv],\n        "protocol": protocol,\n'''
# The first occurrence is restat; the second is full generation. Only full
# generation should claim a replay build environment.
if vtext.count(old) < 2:
    raise SystemExit("validation manifest environment anchors missing")
first = vtext.find(old)
second = vtext.find(old, first + len(old))
vtext = vtext[:second] + vtext[second:].replace(old, new, 1)
validation.write_text(vtext, encoding="utf-8")

robustness = TOOLS / "ou_robustness.py"
rtext = robustness.read_text(encoding="utf-8")
old = '''        "matplotlib": matplotlib.__version__,\n        "command": [sys.executable, *sys.argv],\n        "protocol": protocol,\n'''
new = '''        "matplotlib": matplotlib.__version__,\n        "build_environment": evidence_provenance.environment_metadata(),\n        "command": [sys.executable, *sys.argv],\n        "protocol": protocol,\n'''
if old not in rtext:
    raise SystemExit("robustness full-generation environment anchor missing")
rtext = rtext.replace(old, new, 1)
robustness.write_text(rtext, encoding="utf-8")

print("strengthened replay/restatement provenance safety checks")
