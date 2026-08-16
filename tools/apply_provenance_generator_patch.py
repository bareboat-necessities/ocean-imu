#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace(path: Path, old: str, new: str, count: int = 1) -> None:
    text = path.read_text(encoding="utf-8")
    if text.count(old) < count:
        raise SystemExit(f"anchor not found in {path}: {old[:180]!r}")
    path.write_text(text.replace(old, new, count), encoding="utf-8")


validation = ROOT / "tools/ou_validation.py"
robustness = ROOT / "tools/ou_robustness.py"

replace(
    validation,
    """The runner reports the simulator's own executable quality gates alongside, but\nseparately from, its inferential score.  Validation uses a configurable long\n""",
    """The simulator retains deterministic executable regression gates, but they are\nnot Monte Carlo inclusion criteria and are not exported as statistical-row\nfields. Validation uses a configurable long\n""",
)
replace(
    validation,
    "import numpy as np\n\n\nREPO_ROOT",
    "import numpy as np\n\nimport ou_evidence_provenance as evidence_provenance\n\n\nREPO_ROOT",
)
replace(
    validation,
    '''    with source.open(encoding="utf-8") as stream:\n        bundle = json.load(stream)\n''',
    '''    restatement_context = evidence_provenance.begin_restatement(\n        "validation", source\n    )\n    with source.open(encoding="utf-8") as stream:\n        bundle = json.load(stream)\n''',
)
# Only the restat raw write is the first occurrence after restat_bundle.
text = validation.read_text(encoding="utf-8")
needle = "def restat_bundle("
pos = text.index(needle)
raw_write = "    write_csv(raw_path, raw_rows)\n"
raw_pos = text.index(raw_write, pos)
text = text[:raw_pos] + "    evidence_provenance.preserve_raw_rows(restatement_context, raw_path)\n" + text[raw_pos + len(raw_write):]
validation.write_text(text, encoding="utf-8")
replace(
    validation,
    '''    write_json(manifest_path, manifest)\n\n    for path in (*result_paths, manifest_path):\n''',
    '''    manifest = evidence_provenance.finalize_restatement_manifest(\n        "validation", manifest, restatement_context\n    )\n    write_json(manifest_path, manifest)\n\n    for path in (*result_paths, manifest_path):\n''',
)
replace(
    validation,
    '''                                "quality_gate_pass": int(gate_pass),\n                                "simulator_return_code": return_code,\n''',
    "",
)

replace(
    robustness,
    "import ou_validation as core  # noqa: E402\n",
    "import ou_validation as core  # noqa: E402\nimport ou_evidence_provenance as evidence_provenance  # noqa: E402\n",
)
replace(
    robustness,
    '''            "quality_gate_pass": int(gate_pass),\n            "simulator_return_code": return_code,\n''',
    "",
)
replace(
    robustness,
    '''    with source.open(encoding="utf-8") as stream:\n        bundle = json.load(stream)\n''',
    '''    restatement_context = evidence_provenance.begin_restatement(\n        "robustness", source\n    )\n    with source.open(encoding="utf-8") as stream:\n        bundle = json.load(stream)\n''',
)
text = robustness.read_text(encoding="utf-8")
pos = text.index("def restat_bundle(")
raw_write = "    core.write_csv(raw_path, rows)\n"
raw_pos = text.index(raw_write, pos)
text = text[:raw_pos] + "    evidence_provenance.preserve_raw_rows(restatement_context, raw_path)\n" + text[raw_pos + len(raw_write):]
robustness.write_text(text, encoding="utf-8")
# Replace the old re-pin/warn block with immutable replay + separate restatement provenance.
old = '''    source_files, moved = _restated_source_files(manifest.get("source_files", {}))\n    manifest.update(\n        {\n            "git_commit": core.git_output("rev-parse", "HEAD"),\n            "git_branch": core.git_output("branch", "--show-current"),\n            "git_diff_stat": core.git_output("diff", "--stat"),\n            "python": sys.version,\n            "platform": platform.platform(),\n            "numpy": np.__version__,\n            "command": [sys.executable, *sys.argv],\n            "protocol": protocol,\n            "restated_from": protocol["restated_from"],\n            "source_files": source_files,\n            "result_files": result_files,\n        }\n    )\n    # Only ever present after a restat, and only for sources that actually\n    # moved: the pin above is what the tables were computed under, and this is\n    # what the replays were run under.\n    if moved:\n        manifest["sources_moved_since_rows"] = moved\n    else:\n        manifest.pop("sources_moved_since_rows", None)\n\n    core.write_json(manifest_path, manifest)\n\n    for name in moved:\n        print(\n            f"Note: {name} has changed since the archived rows were produced. "\n            "The rows are carried, not re-run; regenerate the study if that "\n            "file decides what the simulator does."\n        )\n'''
new = '''    manifest.update(\n        {\n            "protocol": protocol,\n            "restated_from": protocol["restated_from"],\n            "result_files": result_files,\n        }\n    )\n    manifest = evidence_provenance.finalize_restatement_manifest(\n        "robustness", manifest, restatement_context\n    )\n    core.write_json(manifest_path, manifest)\n'''
replace(robustness, old, new)

print("patched OU evidence generators")
