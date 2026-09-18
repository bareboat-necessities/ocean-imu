"""Protect the independent article/study source graphs after the proof split."""
from collections import Counter
from pathlib import Path
import re
import unittest

from test_publication_references import (
    DOC, MAIN, STUDY, INPUT_RE, LABEL_RE, reachable_sources, resolve_input,
    strip_comments,
)

# These are the retained detailed proof modules, not copies or rewritten proofs.
PROOF_PARTS = {
    'w3d-init.tex-part',
    'w3d-mems-bias-preconditions.tex-part',
    'w3d-brmm-physical-wave-condition.tex-part',
    'w3d-brmm-stability-theorem.tex-part',
    'w3d-brmm-period-bridge.tex-part',
    'w3d-iss-stability.tex-part',
    'w3d-analytical-stability-widening.tex-part',
    'w3d-stability-widening-phase-c.tex-part',
    'w3d-stability-widening-source-path.tex-part',
    'w3d-p5-staged-capture.tex-part',
    'w3d-semiglobal-stability.tex-part',
    'w3d-finite-live-capture.tex-part',
    'w3d-hybrid-stability.tex-part',
    'w3d-stability-widening-phase-d.tex-part',
    'w3d-stability-widening-phase-e.tex-part',
    'w3d-stability-widening-phase-f.tex-part',
}


def input_counts(root: Path) -> Counter:
    """Count inclusions, retaining multiplicity and rejecting cycles."""
    counts = Counter()

    def visit(path: Path, stack: tuple[Path, ...]) -> None:
        if path in stack:
            raise AssertionError(f'cyclic input: {path}')
        for name in INPUT_RE.findall(strip_comments(path.read_text())):
            child = resolve_input(name)
            if child is not None:
                counts[child.name] += 1
                visit(child, (*stack, path))

    visit(root, ())
    return counts


class StabilityDocumentSplitTests(unittest.TestCase):
    def test_proof_parts_are_in_study_exactly_once_and_not_in_article(self):
        article = input_counts(MAIN)
        study = input_counts(STUDY)
        for part in sorted(PROOF_PARTS):
            with self.subTest(part=part):
                self.assertEqual(article[part], 0)
                self.assertEqual(study[part], 1)

    def test_filter_design_and_results_stay_in_article(self):
        article = input_counts(MAIN)
        study = input_counts(STUDY)
        for part in (
            'w3d-state.tex-part', 'w3d-meas.tex-part',
            'w3d-kalm-up.tex-part', 'w3d-adaptation-motivation.tex-part',
            'w3d-results.tex-part', 'w3d-sim-charts.tex-part',
            'w3d-mag-hard-iron.tex-part', 'w3d-initialization-overview.tex-part',
        ):
            with self.subTest(part=part):
                self.assertEqual(article[part], 1)
                self.assertEqual(study[part], 0)

    def test_study_reader_guide_precedes_detailed_lemmas(self):
        source = STUDY.read_text()
        parts = (
            'w3d-stability-reader-guide.tex-part',
            'w3d-stability-proof-roadmap.tex-part',
            'w3d-stability-certificates.tex-part',
            'w3d-brmm-physical-wave-condition.tex-part',
            'w3d-iss-stability.tex-part',
            'w3d-init.tex-part',
            'w3d-stability-widening-phase-f.tex-part',
        )
        positions = [source.index(r'\input{' + part + '}') for part in parts]
        self.assertEqual(positions, sorted(positions))

    def test_study_sources_exist_and_do_not_depend_on_other_documents_aux(self):
        sources = reachable_sources(STUDY)
        for path, text in sources.items():
            self.assertNotIn(r'\externaldocument', text, path.name)
            for target in re.findall(r'\\studysource\{([^}]+)\}', text):
                self.assertTrue((DOC.parents[1] / target).is_file(), target)
        labels = Counter(label for text in sources.values()
                         for label in LABEL_RE.findall(text))
        self.assertEqual([label for label, count in labels.items() if count > 1], [])

    def test_study_uses_same_ieee_two_column_style_as_article(self):
        article = MAIN.read_text()
        study = STUDY.read_text()
        self.assertIn(r'\\documentclass[conference]{IEEEtran}', article)
        self.assertIn(r'\\documentclass[conference]{IEEEtran}', study)
        self.assertNotIn(r'\\usepackage[margin=27mm]{geometry}', study)
        self.assertIn(r'\\IEEEauthorblockN', study)
        self.assertIn(r'\\IEEEauthorblockA', study)

    def test_both_top_level_documents_are_covered_by_existing_latex_build(self):
        workflow = (DOC.parents[1] / '.github/workflows/build.yml').read_text()
        self.assertIn('doc/kalman_ou_iii', workflow)
        self.assertIn('root_file: "*.tex"', workflow)
        for root in (MAIN, STUDY):
            self.assertIn(r'\begin{document}', root.read_text())
            self.assertIn(r'\end{document}', root.read_text())


if __name__ == '__main__':
    unittest.main()
