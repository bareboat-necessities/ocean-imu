"""Static reference audit for the OU-III publication source.

LaTeX normally treats an unresolved cross-reference or citation as a warning and
still produces a PDF.  For the publication manuscript that is too permissive: a
stripped table, appendix, or bibliography item can silently become ``??``.  This
test walks the main manuscript's reachable source files and requires every
source-level cross-reference to have a reachable label and every citation key to
exist in the manuscript bibliography databases.

Macros are stricter than references.  A cross-reference resolves from the .aux
file, so a label may be defined anywhere in the document, but a macro has to
exist when TeX reads the character that uses it: quoting a generated value in a
section that is typeset before the file defining it is input is a fatal
undefined-control-sequence error, not a warning.  The last test therefore
replays the input order of each compiled document and requires every generated
macro to be defined before it is quoted.
"""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
MAIN = DOC / "kalman_ou-w3d.tex"
RESULTS = DOC / "kalman_ou-w3d-results.tex"

INPUT_RE = re.compile(r"\\input\{([^}]+)\}")
DEFINITION_RE = re.compile(r"\\(?:provide|new|renew)command\s*\{\s*\\([A-Za-z]+)\s*\}")
MACRO_USE_RE = re.compile(r"\\((?:OUValidation|OURobustness)[A-Za-z]*)")
STANDALONE_FLAG = r"\def\OUStandaloneDetailedResults"
STANDALONE_GUARD = r"\ifdefined\OUStandaloneDetailedResults"
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
REF_RE = re.compile(r"\\(?:eqref|ref|cref|Cref|autoref|pageref)\{([^}]+)\}")
CITE_RE = re.compile(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}")
BIBLIOGRAPHY_RE = re.compile(r"\\bibliography\{([^}]+)\}")
BIB_ENTRY_RE = re.compile(r"@\w+\s*\{\s*([^,\s]+)\s*,", re.I)
COMMENT_RE = re.compile(r"(?<!\\)%.*$")


def strip_comments(text: str) -> str:
    return "\n".join(COMMENT_RE.sub("", line) for line in text.splitlines())


def resolve_input(name: str, standalone: bool = False) -> Path | None:
    # The standalone *-study sources are selected only when
    # OUStandaloneDetailedResults is defined, which the main paper never does.
    if name.endswith("-study.tex-part") and not standalone:
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


def typeset_order(root: Path) -> list[tuple[Path, int, str]]:
    """Flatten a document into the line order TeX actually reads.

    Only the two constructs the manuscript uses to select sources are modelled:
    ``\\input`` is expanded in place, and the ``\\ifdefined
    \\OUStandaloneDetailedResults`` header that opens a shared part redirects to
    its *-study source and ends the file when the standalone results document
    is being read.  ``\\IfFileExists`` needs no special handling because its
    guarded input comes before the fallback definitions in the source, which is
    the order TeX sees when the file is present.
    """

    def expand(
        path: Path, standalone: bool, stack: tuple[Path, ...]
    ) -> list[tuple[Path, int, str]]:
        if path in stack:
            return []
        lines = strip_comments(path.read_text(encoding="utf-8")).splitlines()
        flattened: list[tuple[Path, int, str]] = []
        guarded = lines and lines[0].strip().startswith(STANDALONE_GUARD)
        for number, line in enumerate(lines, start=1):
            if guarded:
                if not standalone:
                    # Skip the redirect block; the shared text follows it.
                    if line.strip().startswith(r"\fi"):
                        guarded = False
                    continue
                if line.strip().startswith(r"\expandafter\endinput"):
                    return flattened
            if line.lstrip().startswith(STANDALONE_FLAG):
                standalone = True
            flattened.append((path, number, line))
            for name in INPUT_RE.findall(line):
                child = resolve_input(name, standalone)
                if child is not None:
                    flattened.extend(expand(child, standalone, (*stack, path)))
        return flattened

    return expand(root, False, ())


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

    def test_generated_macros_are_defined_before_they_are_quoted(self):
        # A generated macro quoted ahead of its definition is a fatal
        # undefined-control-sequence error, and it is invisible to a set
        # comparison of definitions against uses: the name is defined, just too
        # late.  That is how a Riccati-section sentence quoting the robustness
        # spans -- defined with the robustness tables, several sections further
        # down -- took the paper build down on main.
        for root in (MAIN, RESULTS):
            with self.subTest(document=root.name):
                defined: set[str] = set()
                early: list[str] = []
                for path, number, line in typeset_order(root):
                    names = DEFINITION_RE.findall(line)
                    body = DEFINITION_RE.sub("", line)
                    for name in MACRO_USE_RE.findall(body):
                        if name not in defined:
                            early.append(f"{path.name}:{number}: \\{name}")
                    defined.update(names)
                self.assertEqual(
                    early,
                    [],
                    "generated macros quoted before they are defined: "
                    f"{early}",
                )

    def test_removed_unevaluated_appendices_are_not_referenced(self):
        sources = reachable_sources()
        combined = "\n".join(sources.values())
        for label in ("sec:heel", "sec:gps_fusion"):
            self.assertNotIn(rf"\ref{{{label}}}", combined)


if __name__ == "__main__":
    unittest.main()
