from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-validation.yml"
BRANCH_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-full-evidence-branch.yml"
BUILD_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "build.yml"
PROOF_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou3-proof.yml"
VALIDATION_MAKEFILE = REPO_ROOT / "tests" / "validation" / "Makefile"


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
        self.assertIn("make -C tests/validation evidence-test", gate)
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
        self.assertIn("make -C tests/validation evidence-test", gate)

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
        validate = loop.index("make -C tests/validation evidence-test")
        self.assertLess(rebase, fingerprint)
        self.assertLess(fingerprint, validate)

    def test_evidence_gate_is_the_suite_without_the_proof_searches(self):
        """What the publication gate runs, and why it is not the whole suite.

        The staged OU-III proof searches are three quarters of an hour of
        interval arithmetic and read source and tooling only: no bundle this
        repository publishes can move their verdict. Running them inside the
        twenty-minute commit job is what killed it mid-suite on every main push
        from 2026-09-03 on, leaving the branch with unpublished evidence and
        skipping the document build behind it. ou3-proof.yml owns them and
        budgets hours; `test` still runs everything for a local full pass and
        for the pull-request smoke gate.
        """
        makefile = VALIDATION_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn(
            "PROOF_SEARCH_TESTS := $(wildcard test_ou3_p2_*.py test_ou3_p3_*.py "
            "test_ou3_p4_*.py test_ou3_p5_*.py)",
            makefile,
        )
        self.assertIn(
            "EVIDENCE_TESTS := $(filter-out $(PROOF_SEARCH_TESTS),"
            "$(wildcard test_*.py))",
            makefile,
        )
        self.assertIn("evidence-test: evidence-contract", makefile)
        self.assertIn("test: evidence-contract", makefile)
        self.assertIn("python3 -m unittest discover -v -p 'test_*.py'", makefile)

        # Most of the skipped modules are named in the proof workflow. The few
        # that are not stay covered because the pull-request smoke gate runs
        # the whole suite, so nothing here drops out of CI entirely.
        proof = PROOF_WORKFLOW.read_text(encoding="utf-8")
        directory = REPO_ROOT / "tests" / "validation"
        skipped = sorted(
            path.stem
            for prefix in ("p2", "p3", "p4", "p5")
            for path in directory.glob(f"test_ou3_{prefix}_*.py")
        )
        self.assertTrue(skipped)
        self.assertTrue(any(name in proof for name in skipped))

        workflow = WORKFLOW.read_text(encoding="utf-8")
        smoke = workflow[workflow.index("  validate:"):workflow.index("  fingerprint:")]
        self.assertIn("make -C tests/validation test", smoke)

    def test_commit_job_validates_against_the_regenerated_commit(self):
        """The bundles are made at github.sha, so the gate must see that tree.

        `regenerate` and `combine` check out this run's commit, and the
        freshness gate in tools/ou_evidence_provenance.py stamps replay
        provenance only when the bundle's recorded git_commit equals HEAD.
        Checking out the moving branch name here made the gate compare the
        bundle against whatever had landed since and refuse, both when the
        branch advanced mid-run and on any re-run of this job.
        """
        workflow = WORKFLOW.read_text(encoding="utf-8")
        commit = workflow.index("  commit:")
        checkout = workflow.index("- name: Checkout", commit)
        skip = workflow.index("- name: Skip when this commit's evidence", commit)
        stage = workflow[checkout:skip]

        self.assertIn("ref: ${{ github.sha }}", stage)
        self.assertNotIn("ref: ${{ github.ref_name }}", stage)
        # The push retry rebases onto the branch tip, which needs real history.
        self.assertIn("fetch-depth: 0", stage)

    def test_republishing_an_already_committed_evidence_commit_is_a_no_op(self):
        """A re-run of `commit` must not try to publish the same evidence twice."""
        workflow = WORKFLOW.read_text(encoding="utf-8")
        commit = workflow.index("  commit:")
        skip = workflow.index("- name: Skip when this commit's evidence", commit)
        install = workflow.index("- name: Install dependencies", skip)
        guard = workflow[skip:install]

        self.assertIn("id: published", guard)
        self.assertIn("replay_provenance", guard)
        self.assertIn('echo "already=true" >> "$GITHUB_OUTPUT"', guard)

        # Every step that validates, stamps, or publishes must be gated on it,
        # so the skip cannot leave a half-applied bundle behind.
        gated = (
            "Install dependencies",
            "Reuse fingerprinted simulation data",
            "Unpack fingerprinted simulation data",
            "Download regenerated bundles",
            "Place bundles and mirror the manuscript copies",
            "Check the manuscript against the regenerated evidence",
            "Record replay and results fingerprints",
            "Verify final evidence fingerprints",
            "Commit the regenerated evidence",
        )
        rest = workflow[install:]
        for name in gated:
            with self.subTest(step=name):
                start = rest.index(f"- name: {name}")
                head = rest[start : start + 400]
                self.assertIn("if: steps.published.outputs.already != 'true'", head)

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
        stage = _job_block(workflow, "build", "release")
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
        stage = _job_block(workflow, "build", "release")
        condition = _compact(_folded_scalar(stage, "if"))
        self.assertEqual(
            condition,
            "${{!cancelled()&&needs.classify.outputs.run_build=='true'&&"
            "(github.ref!='refs/heads/main'||needs['ou-evidence'].result=='success')}}",
        )


if __name__ == "__main__":
    unittest.main()
