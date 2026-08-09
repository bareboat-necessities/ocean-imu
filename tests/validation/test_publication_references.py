"""Static reference audit for the OU-III publication source.

LaTeX normally treats an unresolved cross-reference or citation as a warning and
still produces a PDF.  For the publication manuscript that is too permissive: a
stripped table, appendix, or bibliography item can silently become ``??``.  This
test walks the main manuscript's reachable source files and requires every
source-level cross-reference to have a reachable label and every citation key to
exist in the manuscript bibliography databases.
"""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
MAIN = DOC / "kalman_ou-w3d.tex"

INPUT_RE = re.compile(r"\\input\{([^}]+)\}")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
REF_RE = re.compile(r"\\(?:eqref|ref|cref|Cref|autoref|pageref)\{([^}]+)\}")
CITE_RE = re.compile(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}")
BIBLIOGRAPHY_RE = re.compile(r"\\bibliography\{([^}]+)\}")
BIB_ENTRY_RE = re.compile(r"@\w+\s*\{\s*([^,\s]+)\s*,", re.I)
COMMENT_RE = re.compile(r"(?<!\\)%.*$")


def strip_comments(text: str) -> str:
    return "\n".join(COMMENT_RE.sub("", line) for line in text.splitlines())


def resolve_input(name: str) -> Path | None:
    # The standalone *-study sources are selected only when
    # OUStandaloneDetailedResults is defined, which the main paper never does.
    if name.endswith("-study.tex-part"):
        return None
    path = DOC / name
    if path.is_file():
        return path
    for suffix in (".tex", ".tex-part"):
        candidate = DOC / f"{name}{suffix}"
        if candidate.is_file():
            return candidate
    return None


def reachable_sources() -> dict[Path, str]:
    pending = [MAIN]
    sources: dict[Path, str] = {}
    while pending:
        path = pending.pop()
        if path in sources:
            continue
        text = strip_comments(path.read_text(encoding="utf-8"))
        sources[path] = text
        for name in INPUT_RE.findall(text):
            child = resolve_input(name)
            if child is not None and child not in sources:
                pending.append(child)
    return sources


def bibliography_keys(main_source: str) -> set[str]:
    keys: set[str] = set()
    for group in BIBLIOGRAPHY_RE.findall(main_source):
        for name in group.split(","):
            path = DOC / f"{name.strip()}.bib"
            if not path.is_file():
                raise AssertionError(f"missing bibliography database: {path.name}")
            keys.update(BIB_ENTRY_RE.findall(strip_comments(path.read_text(encoding="utf-8"))))
    return keys


class PublicationReferenceTests(unittest.TestCase):
    def test_every_reachable_cross_reference_has_a_reachable_label(self):
        sources = reachable_sources()
        labels = {
            label
            for text in sources.values()
            for label in LABEL_RE.findall(text)
        }
        missing: list[tuple[str, str]] = []
        for path, text in sources.items():
            for group in REF_RE.findall(text):
                for target in (item.strip() for item in group.split(",")):
                    if target and target not in labels:
                        missing.append((path.name, target))
        self.assertEqual(missing, [], f"unresolved publication references: {missing}")

    def test_every_reachable_citation_has_a_bibliography_entry(self):
        sources = reachable_sources()
        available = bibliography_keys(sources[MAIN])
        missing: list[tuple[str, str]] = []
        for path, text in sources.items():
            for group in CITE_RE.findall(text):
                for key in (item.strip() for item in group.split(",")):
                    if key and key not in available:
                        missing.append((path.name, key))
        self.assertEqual(missing, [], f"undefined publication citations: {missing}")

    def test_removed_unevaluated_appendices_are_not_referenced(self):
        sources = reachable_sources()
        combined = "\n".join(sources.values())
        for label in ("sec:heel", "sec:gps_fusion"):
            self.assertNotIn(rf"\ref{{{label}}}", combined)


if __name__ == "__main__":
    unittest.main()
