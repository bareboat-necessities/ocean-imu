#!/usr/bin/env python3
"""Keep generated OU publication text aligned with the manuscript claim scope.

The statistical generator deliberately owns numbers and tables.  Editorial
claim scope is owned by the article.  This small post-generation step removes a
retired interpretation sentence from the generated 3-D table caption and then
refreshes the manifest record for the publication file.  It changes no replay
rows or statistics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


RETIRED_AXIS_CLAIM = (
    " The vertical gain of OU--III is paid for in the horizontal channels."
)
CURRENT_CAPTION_TAIL = r"positive $\Delta$ favors OU--II.}"
PUBLICATION_NAME = "ou_validation_publication.tex"
MANIFEST_NAME = "ou_validation_manifest.json"


def _file_record(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def sync_validation_publication(validation_dir: Path) -> bool:
    publication = validation_dir / PUBLICATION_NAME
    manifest_path = validation_dir / MANIFEST_NAME
    if not publication.is_file():
        raise FileNotFoundError(publication)
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)

    text = publication.read_text(encoding="utf-8")
    count = text.count(RETIRED_AXIS_CLAIM)
    if count > 1:
        raise RuntimeError(
            f"retired publication claim occurs {count} times; refusing ambiguous rewrite"
        )
    if count == 1:
        text = text.replace(RETIRED_AXIS_CLAIM, "", 1)
        publication.write_text(text, encoding="utf-8")
    if CURRENT_CAPTION_TAIL not in text:
        raise RuntimeError(
            "OU 3-D publication caption no longer matches the article wording contract"
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    result_files = manifest.get("result_files")
    if not isinstance(result_files, dict) or PUBLICATION_NAME not in result_files:
        raise RuntimeError(
            f"manifest result_files has no {PUBLICATION_NAME} record"
        )
    new_record = _file_record(publication)
    changed = result_files[PUBLICATION_NAME] != new_record or count == 1
    result_files[PUBLICATION_NAME] = new_record
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return changed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, required=True)
    args = parser.parse_args()
    changed = sync_validation_publication(args.validation_dir)
    print(
        "Aligned OU validation publication wording and manifest."
        if changed
        else "OU validation publication wording already aligned."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
