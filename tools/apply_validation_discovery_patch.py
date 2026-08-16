#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
core = ROOT / "tests/validation/test_ou_validation.py"
compat = ROOT / "tests/validation/test_reorg_compat_overrides.py"
pub = ROOT / "tests/validation/test_zz_publication_contract.py"

text = core.read_text(encoding="utf-8")

def replace_method(name: str, module: str, helper: str) -> None:
    global text
    pattern = re.compile(
        rf"^    def {re.escape(name)}\(self\):\n.*?(?=^    def |^class |\Z)",
        re.MULTILINE | re.DOTALL,
    )
    replacement = (
        f"    def {name}(self):\n"
        f"        from {module} import {helper}\n"
        f"        {helper}(self)\n\n"
    )
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise SystemExit(f"expected exactly one method {name}, found {count}")

for name, module, helper in (
    ("test_the_deterministic_table_is_generated_not_typed", "test_reorg_compat_overrides", "_deterministic_table_remains_generated_evidence"),
    ("test_singer_relationship_and_contribution_wording", "test_reorg_compat_overrides", "_singer_relationship_and_contribution_boundary"),
    ("test_tracker_and_heel_material_are_scoped_as_unevaluated_appendices", "test_reorg_compat_overrides", "_unevaluated_tracker_and_heel_material_stays_out_of_main_claims"),
    ("test_abstract_reports_committed_stationary_aggregate", "test_zz_publication_contract", "_abstract_reports_committed_stationary_aggregate"),
    ("test_fixed_reference_and_transition_limits_are_stated", "test_zz_publication_contract", "_fixed_reference_and_transition_limits_are_stated"),
    ("test_inference_is_qualified_rather_than_asserted", "test_zz_publication_contract", "_inference_is_qualified_rather_than_asserted"),
    ("test_three_dimensional_and_channel_results_are_reported", "test_zz_publication_contract", "_three_dimensional_and_channel_results_are_reported"),
    ("test_transition_and_secondary_ensembles_are_rescored", "test_zz_publication_contract", "_transition_and_secondary_ensembles_are_rescored"),
    ("test_baseline_fairness_thresholds_and_hardware_limits_are_recorded", "test_zz_publication_contract", "_baseline_fairness_thresholds_and_hardware_limits_are_recorded"),
    ("test_the_contribution_is_framed_at_the_width_of_the_evidence", "test_zz_publication_contract", "_contribution_is_framed_at_the_width_of_the_evidence"),
):
    replace_method(name, module, helper)
core.write_text(text, encoding="utf-8")

# The helper modules remain ordinary semantic-check libraries (and the
# publication module keeps its real TestCase class); they no longer mutate the
# test classes at import time, so discovery order is irrelevant.
text = compat.read_text(encoding="utf-8")
text = re.sub(
    r"\ncore\.ManuscriptMethodologyTests\.test_the_deterministic_table_is_generated_not_typed.*\Z",
    "\n",
    text,
    flags=re.DOTALL,
)
compat.write_text(text, encoding="utf-8")

text = pub.read_text(encoding="utf-8")
start = text.find("\ncore.CommittedFullResultsTests.test_abstract_reports_committed_stationary_aggregate")
end = text.find("\n\nclass ReorganizedPublicationContractTests", start)
if start < 0 or end < 0:
    raise SystemExit("publication monkeypatch block not found")
text = text[:start] + text[end:]
pub.write_text(text, encoding="utf-8")

print("converted publication tests to vanilla-discovery semantics")
