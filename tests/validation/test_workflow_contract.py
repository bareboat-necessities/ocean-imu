#!/usr/bin/env python3
from pathlib import Path
import subprocess
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-validation.yml"
BUILD_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "build.yml"
TFG_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "tfg-validation.yml"
PROOF_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou3-proof.yml"
VALIDATION_MAKEFILE = REPO_ROOT / "tests" / "validation" / "Makefile"
ROOT_MAKEFILE = REPO_ROOT / "Makefile"


class WorkflowContractTests(unittest.TestCase):
    def test_validation_workflow_is_smoke_only_and_does_not_commit(self):
        """Standalone validation stays smoke-only; reusable main owns evidence."""
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("python3-matplotlib", text)
        self.assertIn("libeigen3-dev", text)
        self.assertIn("build-for-validation", text)
        self.assertIn("python3 tools/ou_validation.py --repetitions 1", text)
        self.assertIn("python3 tools/ou_robustness.py --repetitions 1", text)
        self.assertIn("--output /tmp/ou_study.json", text)
        self.assertIn("--generated-dir /tmp/ou_generated", text)
        self.assertIn("python3 tools/ou_rs_law_ablation.py", text)
        self.assertIn("--output /tmp/ou_rs_law_ablation.json", text)
        self.assertIn("--generated-dir /tmp/ou_rs_law_generated", text)
        self.assertIn("--refactor-study-json /tmp/ou_refactor.json", text)
        self.assertIn("--bootstrap-resamples 200", text)
        self.assertIn("--snapshots-dir /tmp/ou_snapshots", text)
        self.assertIn("--output-dir /tmp/ou_detrending", text)
        self.assertIn("make -C tests/validation test", text)
        self.assertIn("name: ou-validation-smoke", text)
        self.assertNotIn("git pull --rebase", text)
        self.assertNotIn("push:\n", text)
        self.assertIn("workflow_call:", text)
        self.assertIn("validation_mode:", text)
        self.assertIn("default: smoke", text)
        self.assertIn("workflow_dispatch:", text)
        self.assertIn("if: inputs.validation_mode == 'full'", text)
        self.assertIn("permissions:\n  contents: read", text)
        self.assertIn("contents: write", text)

    def test_main_build_owns_full_validation_artifacts(self):
        build = BUILD_WORKFLOW.read_text(encoding="utf-8")
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("ou-evidence:", build)
        self.assertIn("uses: ./.github/workflows/ou-validation.yml", build)
        self.assertIn("validation_mode: ${{ github.ref == 'refs/heads/main' && 'full' || 'smoke' }}", build)
        self.assertIn("stage an equivalent fresh bundle before the test prerequisite", build)
        self.assertIn("python3 tools/ou_validation.py --repetitions 10", workflow)
        self.assertIn("python3 tools/ou_robustness.py --repetitions 10", workflow)
        self.assertIn("--bootstrap-resamples 10000", workflow)
        self.assertIn("--no-figures", workflow)
        self.assertIn("--benchmark-only", workflow)
        self.assertIn("--refresh-only", workflow)
        self.assertIn("--from-json", workflow)
        self.assertIn("--combine-shards", workflow)
        self.assertIn("--shard-index ${{ matrix.shard }}", workflow)
        self.assertIn("shard: [0, 1, 2]", workflow)
        self.assertIn("--snapshot-prefix matrix_", workflow)
        self.assertIn("--snapshot-refactor-study-json", workflow)
        self.assertIn("name: ou-validation-full", workflow)
        self.assertIn("name: ou-robustness-full", workflow)
        self.assertIn("gh artifact", workflow)
        self.assertIn("git push origin HEAD:${GITHUB_REF_NAME}", workflow)
        self.assertIn("generated-evidence-${{ github.ref }}", workflow)
        self.assertIn("--from-json reports/results/ou_validation/ou_validation.json", workflow)
        self.assertIn("--from-json reports/results/ou_robustness/ou_robustness.json", workflow)
        self.assertIn("matrix_ou_refactor.json", workflow)
        self.assertIn("matrix_ou_detrending.json", workflow)
        self.assertIn("matrix_sigma_ablation.csv", workflow)
        self.assertIn("Record replay and results fingerprints", workflow)
        self.assertIn("--write reports/ou_evidence_fingerprint.json", workflow)
        self.assertIn("--check reports/ou_evidence_fingerprint.json", workflow)
        self.assertIn("git add reports/results/ou_validation reports/results/ou_robustness reports/ou_evidence_fingerprint.json", workflow)

    def test_full_evidence_jobs_reuse_the_fingerprinted_simulation_archive(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        # Only the fingerprint job may resolve a release archive. Replay shards,
        # combine and commit must consume those exact uploaded archive bytes.
        self.assertEqual(workflow.count("- name: Download versioned simulation data"), 1)
        self.assertEqual(workflow.count("/releases/tags/simulation"), 1)
        self.assertIn("name: ou-fingerprinted-simulation-data", workflow)
        self.assertGreaterEqual(workflow.count("- name: Reuse fingerprinted simulation data"), 3)
        self.assertGreaterEqual(workflow.count("gh artifact download --repo"), 5)

    def test_validation_workflow_combines_each_study_on_its_own_runner(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        combine_job = workflow.split("\n  combine:\n", 1)[1].split("\n  commit:\n", 1)[0]
        self.assertIn("study: [validation, robustness]", combine_job)
        self.assertIn('"ou-${{ matrix.study }}-shard-*"', combine_job)
        self.assertNotIn('--name "ou-*-shard-*"', combine_job)
        self.assertIn("name: ou-${{ matrix.study }}-full", combine_job)

    def test_evidence_commit_materializes_provenance_before_running_tests(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        commit_job = text.split("\n  commit:\n", 1)[1].split("\n  validate:\n", 1)[0]
        self.assertIn("make -C tests/validation evidence-test", commit_job)
        self.assertIn("Skip when this commit's evidence is already published", commit_job)
        self.assertIn("if git log --format=%s", commit_job)
        self.assertIn("replay_commit", commit_job)
        makefile = VALIDATION_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn("evidence-test: build evidence-contract", makefile)
        self.assertIn("evidence-contract: publication-sync", makefile)
        self.assertIn("ou_evidence_contract.py --auto", makefile)
        self.assertIn("ou_publication_sync.py --sync", makefile)
        self.assertIn("ou_evidence_contract.py --check", makefile)

    def test_evidence_gate_is_partitioned_from_expensive_proof_search(self):
        """What the publication gate runs, and why it is not the whole suite.

        Full make test remains the development/CI target, but publishing an
        unchanged evidence bundle must not also rerun half-hour interval word
        searches. Dedicated proof CI owns those tests. Keep this contract near
        the workflow rather than scattering ad-hoc -k/-p lists across YAML.
        """
        makefile = VALIDATION_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn(
            "PROOF_SEARCH_TESTS := $(wildcard test_ou3_p2_*.py test_ou3_p3_*.py "
            "test_ou3_p4_*.py test_ou3_p5_*.py test_ou3_sea3_*.py)",
            makefile,
        )
        self.assertIn(
            "EVIDENCE_TESTS := $(filter-out $(PROOF_SEARCH_TESTS),$(wildcard test_*.py))",
            makefile,
        )
        self.assertIn("python3 -m unittest -v $(EVIDENCE_TESTS:.py=)", makefile)
        self.assertIn("test: build evidence-contract", makefile)
        self.assertIn("python3 -m unittest discover -v -p 'test_*.py'", makefile)

        # Most of the skipped modules are named in the proof workflow. Some
        # low-level primitives also run in source-foundation; none is skipped
        # from the ordinary full test target.
        proof = PROOF_WORKFLOW.read_text(encoding="utf-8")
        skipped = [
            "test_ou3_p2_correlation_path_memory",
            "test_ou3_p3_correlated_translation_covariance_upper",
            "test_ou3_p3_correlated_translation_segment",
            "test_ou3_p3_ltv_postmeasurement_certificate",
            "test_ou3_p3_ltv_translation_ucc_probe",
            "test_ou3_p3_source_uniform_translation_covariance",
        ]
        # This is intentionally a best-effort declaration check. The full suite
        # continues to run every module, while dedicated proof CI names the
        # expensive certificate obligations it owns.
        self.assertTrue(any(module in proof for module in skipped))
        workflow = WORKFLOW.read_text(encoding="utf-8")
        evidence_test_job = workflow.split("\n  commit:\n", 1)[1].split("\n  validate:\n", 1)[0]
        self.assertNotIn("make -C tests/validation test", evidence_test_job)
        self.assertIn("make -C tests/validation evidence-test", evidence_test_job)

    def test_main_build_keeps_compile_and_run_shards(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("Build test executables only", text)
        self.assertIn("${{ matrix.test_dir }} build", text)
        self.assertIn("Run precompiled test shard", text)
        self.assertIn("SKIP_EVIDENCE_CONTRACT=1", text)
        self.assertIn("Enforce evidence provenance before validation shard tests", text)
        self.assertIn("Regenerate OU-III motion diagnostics from this commit", text)
        self.assertIn("python3 tools/ou3_motion_evidence.py", text)
        self.assertIn("--repo . --generated-dir doc/kalman_ou_iii --rooted-relative-paths", text)
        self.assertIn("motion diagnostics stayed fresh but OU replay dependencies changed", text)
        self.assertIn("--rebind-cached-snapshots", text)
        self.assertNotIn("make -C tests/validation evidence-contract", text)
        self.assertIn("collect-results:", text)
        self.assertIn("Download test results", text)
        self.assertIn("Download the exact test source checkout", text)
        self.assertIn("Require verified source archive provenance", text)
        self.assertIn("source-metadata", text)
        self.assertIn("--fail-fast", text)
        self.assertIn("source_archive_provenance.py verify-extracted", text)
        self.assertIn("source_archive_provenance.py create", text)
        self.assertIn("source_archive_provenance.py require-clean", text)
        self.assertIn("Required generated result is absent from the source worktree", text)
        self.assertNotIn('cp -f "results/$f"', text)
        self.assertIn("Collect OU-III multi-seed artifacts produced by test shards", text)
        self.assertNotIn("Rebuild OU-III article ablation panels from validated sources", text)
        self.assertNotIn("Rebuild OU-III matched-law R_S panels from validated sources", text)
        self.assertNotIn("Rebuild OU-III multi-seed evidence from validated sources", text)

    def test_main_build_preserves_matching_materialized_evidence_for_papers(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        collect = text.split("\n  collect-results:\n", 1)[1].split("\n  build-MCU:\n", 1)[0]
        artifact_restore = collect.index('for snapshot in /tmp/ou-final-evidence')
        materialized_check = collect.index('python3 tools/ou_evidence_contract.py --check')
        historical_refresh = collect.index('MIGRATE=(--auto --migrate-from-history')
        self.assertLess(artifact_restore, materialized_check)
        self.assertLess(materialized_check, historical_refresh)
        self.assertIn('publication already emitted by this run', collect)
        self.assertIn('if (( replay_materialized == 0 )); then', collect)

    def test_main_build_downloads_test_results_without_evidence_bundle_collisions(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        collect = text.split("\n  collect-results:\n", 1)[1].split("\n  build-MCU:\n", 1)[0]
        self.assertIn('mapfile -t RESULT_ARTIFACTS', collect)
        self.assertIn('select(startswith("result-"))', collect)
        self.assertIn('artifact_args+=(--name "$artifact")', collect)
        self.assertIn('gh artifact download', collect)
        self.assertIn('"${artifact_args[@]}" --dir /tmp/results', collect)

    def test_compile_check_verifies_the_fingerprint_used_by_papers(self):
        makefile = VALIDATION_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn("ou_evidence_contract.py --auto", makefile)
        self.assertIn("ou_publication_sync.py --sync", makefile)
        self.assertIn("ou_evidence_contract.py --check", makefile)
        text = TFG_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("command: make -C tests/validation test", text)
        self.assertIn('"doc/kalman_ou_iii/ou3_*.svg"', text)
        self.assertIn('"doc/kalman_ou_iii/ou3_*generated.tex-part"', text)
        self.assertIn('"doc/kalman_ou_iii/*tuning*generated.tex-part"', text)

    def test_main_only_jobs_remain_main_only(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("if: github.ref == 'refs/heads/main'", text)
        self.assertIn("Publish Release", text)
        self.assertIn("make-release:", text)

    def test_tfg_validation_avoids_old_fingerprint_stamping(self):
        text = TFG_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("make -C tests/validation test", text)
        self.assertNotIn("generated/", text)
        self.assertNotIn("ou_replay_fingerprint.py --write", text)
        self.assertNotIn("kalman_ou2.hpp", text)
        self.assertNotIn("kalman_ou3.hpp", text)

    def test_main_only_math_typesetting_uses_main(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("make-papers:", text)
        papers = text.split("\n  make-papers:\n", 1)[1].split("\n  make-release:\n", 1)[0]
        self.assertIn("github.ref == 'refs/heads/main'", papers)
        self.assertIn("[build, run-tests, ou-evidence, collect-results]", papers)
        self.assertIn("Download the exact test source checkout", papers)
        self.assertIn("Copy in the recorded source checkout", papers)
        self.assertIn("source_archive_provenance.py verify-extracted", papers)
        self.assertIn("type: Float64x2", papers)

    def test_branches_still_receive_smoke_validation(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("uses: ./.github/workflows/ou-validation.yml", text)
        self.assertIn("validation_mode: ${{ github.ref == 'refs/heads/main' && 'full' || 'smoke' }}", text)
        self.assertIn("name: classify latest push", text)
        self.assertIn("needs.classify.outputs.replay_inputs_changed == 'true'", text)
        self.assertIn("no broad build needed for proof/docs-only push", text)

    def test_smoke_validation_skips_simulations_for_proof_and_doc_pushes(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("name: classify latest push", text)
        self.assertIn("git diff --name-only \"$BASE\" \"$HEAD\"", text)
        self.assertIn("src/*|tests/*|Makefile|.github/workflows/build.yml|.github/workflows/ou-validation.yml", text)
        self.assertIn("tools/ou3_*|tools/stability/*", text)
        self.assertIn("name: no simulation smoke needed", text)
        self.assertIn("proof-only/doc-only changes do not launch simulator smoke jobs", text)
        self.assertIn("concurrency:\n  group: ou-validation-${{ github.ref }}", text)
        self.assertIn("cancel-in-progress: true", text)

    def test_root_make_all_reuses_single_sim_archive_download(self):
        text = ROOT_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn("if [ -n \"$$missing\" ]; then", text)
        self.assertIn("for d in $$missing; do", text)
        self.assertIn("Found existing", text)


    def test_resolved_partition_retains_every_proof_and_evidence_module(self):
        directory = REPO_ROOT / "tests" / "validation"
        result = subprocess.run(
            ["make", "--no-print-directory", "-s", "print-test-partition"],
            cwd=directory, check=True, capture_output=True, text=True, timeout=30,
        )
        rows = dict(line.split("=", 1) for line in result.stdout.splitlines())
        self.assertEqual(set(rows), {"EVIDENCE", "PROOF"})
        groups = {key: value.split() for key, value in rows.items()}
        all_modules = {p.stem for p in directory.glob("test_*.py")}
        proof = {p.stem for prefix in ("p2", "p3", "p4", "p5", "sea3")
                 for p in directory.glob(f"test_ou3_{prefix}_*.py")}
        self.assertTrue(proof)
        self.assertEqual(set(groups["PROOF"]), proof)
        self.assertEqual(set(groups["EVIDENCE"]), all_modules - proof)
        for group in groups.values():
            self.assertEqual(len(group), len(set(group)))


if __name__ == "__main__":
    unittest.main()
