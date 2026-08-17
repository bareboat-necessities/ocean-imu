from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
ARTICLE = REPO_ROOT / "doc" / "kalman_ou_iii" / "w3d-baseline-comparison.tex-part"
GENERATOR = REPO_ROOT / "plots" / "kalman_ou_iii" / "baseline-comparison.py"


class DeterministicBaselinePublicationContractTests(unittest.TestCase):
    def test_prose_uses_generator_owned_aggregate_macros(self):
        article = ARTICLE.read_text(encoding="utf-8")
        generator = GENERATOR.read_text(encoding="utf-8")
        macros = (
            "OUBaselinePIIMeanVerticalPercent",
            "OUBaselineTVGNLOMeanVerticalPercent",
            "OUBaselineOUIIMeanVerticalPercent",
            "OUBaselineOUIIIMeanVerticalPercent",
        )
        for macro in macros:
            self.assertIn(rf"\{macro}", article)
            self.assertIn(macro, generator)

        for stale in (r"\pct{4.5}", r"\pct{6.5}", r"\pct{11.9}", r"\pct{14.5}"):
            self.assertNotIn(stale, article)

    def test_macro_values_come_from_same_aggregate_as_table(self):
        generator = GENERATOR.read_text(encoding="utf-8")
        self.assertIn('value = tex_num(aggregates[method]["mean_z_pct"], 1)', generator)
        self.assertIn('tex_num(a[\'mean_z_pct\'],1)', generator)
        self.assertIn(r'rf"\renewcommand{{\{macro}}}{{{value}}}"', generator)


if __name__ == "__main__":
    unittest.main()
