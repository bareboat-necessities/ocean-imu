"""Experimental overlays reproduce exact headers and reject stale sources."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
import encounter_pseudo_study as study


class PseudoRecipeTests(unittest.TestCase):
    def setUp(self):
        snapshot = json.loads(study.RECIPES.read_text())
        self.recipes = snapshot['recipes']
        self.sources = snapshot['source_texts']

    def test_every_recorded_overlay_reconstructs_exact_bytes(self):
        for name, recipe in self.recipes.items():
            with self.subTest(recipe=name), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                study.apply_recipe(recipe, directory, self.sources)
                self.assertEqual({str(p.relative_to(directory)) for p in directory.rglob('*.h')},
                                 set(recipe['patches']) | set(recipe['added']))
                for path, patch in recipe['patches'].items():
                    self.assertEqual(study.digest(directory/path), patch['overlay_sha256'])

    def test_stale_source_or_corrupted_overlay_is_rejected(self):
        for field in ('base_sha256', 'overlay_sha256'):
            recipe = copy.deepcopy(self.recipes['wave_matched:ou3'])
            next(iter(recipe['patches'].values()))[field] = '0'*64
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                with self.assertRaises(ValueError):
                    study.apply_recipe(recipe, Path(temporary), self.sources)

    def test_overlapping_edits_are_rejected(self):
        recipe = copy.deepcopy(self.recipes['wave_matched:ou3'])
        edits = next(iter(recipe['patches'].values()))['edits']
        edits.append(copy.deepcopy(edits[0]))
        with tempfile.TemporaryDirectory() as temporary, self.assertRaises(ValueError):
            study.apply_recipe(recipe, Path(temporary), self.sources)


if __name__ == '__main__':
    unittest.main()
