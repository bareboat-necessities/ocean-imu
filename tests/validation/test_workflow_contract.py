from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-validation.yml"
BRANCH_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-full-evidence-branch.yml"
BUILD_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "build.yml"


def _mapping_child_keys(text, section):
    lines = text.splitlines()
    start = lines.index(f"{section}:")
    keys = set()
    for line in lines[start + 1 :]:
        if line and not line.startswith(" "):
            break
        if not line.startswith("  ") or line.startswith("    "):
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or ":" not in stripped:
            continue
        keys.add(stripped.split(":", 1)[0])
    return keys


def _job_block(workflow, job_name, next_job_name):
    start = workflow.index(f"  {job_name}:")
    end = workflow.index(f"  {next_job_name}:", start)
    return workflow[start:end]


def _inline_sequence(stage, key):
    prefix = f"{key}: ["
    for line in stage.splitlines():
        stripped = line.strip()
        if stripped.startswith(prefix) and stripped.endswith("]"):
            values = stripped[len(prefix) : -1]
            return [value.strip().strip("'\"") for value in values.split(",")]
    raise AssertionError(f"inline sequence {key!r} not found")


def _folded_scalar(stage, key):
    lines = stage.splitlines()
    marker = f"{key}: >-"
    for index, line in enumerate(lines):
        if line.strip() != marker:
            continue
        indent = len(line) - len(line.lstrip())
        parts = []
        for continuation in lines[index + 1 :]:
            if not continuation.strip():
                continue
            continuation_indent = len(continuation) - len(continuation.lstrip())
            if continuation_indent <= indent:
                break
            parts.append(continuation.strip())
        return " ".join(parts)
    raise AssertionError(f"folded scalar {key!r} not found")


def _compact(expression):
    return "".join(expression.split())


class WorkflowContractTests(unittest.TestCase):
    def test_full_regeneration_validates_before_commit(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        check_marker = "- name: Check the manuscript against the regenerated evidence"
        commit_marker = "- name: Commit the regenerated evidence"
        check = workflow.index(check_marker)
        commit = workflow.index(commit_marker)

        self.assertLess(check, commit)
        gate = workflow[check:commit]
        self.assertIn("make -C tests/validation test", gate)
        self.assertNotIn("continue-on-error", gate)

    def test_validation_publication_is_aligned_before_bundle_upload(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        combine = workflow.index("- name: Combine paired validation shards")
        align = workflow.index(
            "- name: Align validation publication wording with the article"
        )
        upload = workflow.index("- name: Upload regenerated bundle")

        self.assertLess(combine, align)
        self.assertLess(align, upload)
        stage = workflow[align:upload]
        self.assertIn("tools/ou_publication_sync.py", stage)
        self.assertIn("reports/results/ou_validation", stage)

    def test_full_replay_is_gated_by_replay_and_results_fingerprints(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        fingerprint = workflow.index("  fingerprint:")
        regenerate = workflow.index("  regenerate:")
        self.assertLess(fingerprint, regenerate)

        gate = workflow[fingerprint:regenerate]
        self.assertIn("tools/ou_replay_fingerprint.py", gate)
        self.assertIn("sim-data-files.zip", gate)
        self.assertIn("reports/ou_evidence_fingerprint.json", gate)
        self.assertIn("replay_required=false", gate)
        self.assertIn("replay_required=true", gate)
        self.assertIn("complete results tree", gate)
        self.assertIn("make -C tests/validation test", gate)

        regen_header = workflow[regenerate:workflow.index("    runs-on:", regenerate)]
        self.assertIn("needs: fingerprint", regen_header)
        self.assertIn(
            "needs.fingerprint.outputs.replay_required == 'true'", regen_header
        )

    def test_every_full_replay_uses_the_fingerprinted_zip_bytes(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        full = workflow[workflow.index("  fingerprint:"):]
        self.assertIn("- name: Preserve fingerprinted simulation archive", full)
        self.assertIn("name: ou-fingerprinted-simulation-data", full)
        self.assertGreaterEqual(
            full.count("name: ou-fingerprinted-simulation-data"),
            4,
        )
        regenerate = full[full.index("  regenerate:"):]
        self.assertNotIn("gh release download v1.1.3", regenerate)

    def test_regenerated_evidence_hashes_the_final_validated_tree(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        check = workflow.index("- name: Check the manuscript against the regenerated evidence")
        record = workflow.index("- name: Record replay and results fingerprints")
        verify = workflow.index("- name: Verify final evidence fingerprints")
        commit = workflow.index("- name: Commit the regenerated evidence")
        self.assertLess(check, record)
        self.assertLess(record, verify)
        self.assertLess(verify, commit)

        stage = workflow[record:commit]
        self.assertIn("tools/ou_replay_fingerprint.py", stage)
        self.assertIn("--write reports/ou_evidence_fingerprint.json", stage)
        self.assertIn("--check reports/ou_evidence_fingerprint.json", stage)

        commit_stage = workflow[commit:]
        self.assertIn("reports/ou_evidence_fingerprint.json", commit_stage)
        self.assertIn("reports/results/ou_validation", commit_stage)
        self.assertIn("reports/results/ou_robustness", commit_stage)

    def test_results_tree_changes_are_in_workflow_trigger(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn('- "reports/results/**"', workflow)

    def test_evidence_workflow_changes_trigger_smoke_validation(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn('- ".github/workflows/ou-validation.yml"', workflow)
        self.assertIn('- ".github/workflows/build.yml"', workflow)
        self.assertIn('- ".github/workflows/ou-full-evidence-branch.yml"', workflow)

    def test_push_retry_revalidates_after_rebase(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        start = workflow.index("for attempt in 1 2 3 4; do")
        end = workflow.index("\n          done", start)
        loop = workflow[start:end]

        rebase = loop.index("git pull --rebase")
        fingerprint = loop.index("tools/ou_replay_fingerprint.py")
        validate = loop.index("make -C tests/validation test")
        self.assertLess(rebase, fingerprint)
        self.assertLess(fingerprint, validate)

    def test_failure_message_cannot_run_after_a_push(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("The regenerated bundle is committed, but", workflow)
        self.assertNotIn("Fail if the manuscript no longer matches", workflow)

    def test_branch_full_evidence_is_manual_only(self):
        workflow = BRANCH_WORKFLOW.read_text(encoding="utf-8")
        self.assertEqual(_mapping_child_keys(workflow, "on"), {"workflow_dispatch"})
        self.assertIn("validation_mode: full", workflow)

    def test_main_build_is_the_authoritative_automatic_full_evidence_path(self):
        workflow = BUILD_WORKFLOW.read_text(encoding="utf-8")
        evidence = workflow.index("  ou-evidence:")
        build = workflow.index("  build:", evidence)
        stage = workflow[evidence:build]
        self.assertIn("github.ref == 'refs/heads/main'", stage)
        self.assertIn("uses: ./.github/workflows/ou-validation.yml", stage)
        self.assertIn("validation_mode: full", stage)

    def test_main_pdf_build_uses_post_evidence_head_and_compiles_ou_iii(self):
        workflow = BUILD_WORKFLOW.read_text(encoding="utf-8")
        stage = _job_block(workflow, "build", "ou-tuning")
        self.assertIn("needs: [ou-evidence, classify]", stage)
        self.assertIn("kalman_ou_iii", _inline_sequence(stage, "dir"))
        self.assertIn(
            "ref: ${{ github.ref == 'refs/heads/main' && 'refs/heads/main' || '' }}",
            stage,
        )
        self.assertIn("- name: Compile LaTeX document (${{ matrix.dir }})", stage)
        self.assertIn("working_directory: doc/${{ matrix.dir }}", stage)

    def test_main_pdf_build_requires_successful_evidence(self):
        workflow = BUILD_WORKFLOW.read_text(encoding="utf-8")
        stage = _job_block(workflow, "build", "ou-tuning")
        condition = _compact(_folded_scalar(stage, "if"))
        self.assertEqual(
            condition,
            "${{!cancelled()&&needs.classify.outputs.run_build=='true'&&"
            "(github.ref!='refs/heads/main'||needs['ou-evidence'].result=='success')}}",
        )


if __name__ == "__main__":
    unittest.main()
