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


def _job(text: str, name: str, next_name: str | None = None) -> str:
    start = text.index(f"\n  {name}:\n")
    if next_name is not None:
        end = text.index(f"\n  {next_name}:\n", start + 1)
        return text[start:end]
    return text[start:]


class WorkflowContractTests(unittest.TestCase):
    def test_validation_workflow_separates_smoke_from_full_publication(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("workflow_call:", text)
        self.assertIn("workflow_dispatch:", text)
        self.assertIn("default: smoke", text)
        self.assertIn("permissions:\n  contents: read", text)
        validate = _job(text, "validate", "fingerprint")
        self.assertIn("if: inputs.validation_mode != 'full'", validate)
        self.assertIn("make -C tests/validation test", validate)
        self.assertIn("--mode smoke", validate)
        self.assertNotIn("git push", validate)
        fingerprint = _job(text, "fingerprint", "regenerate")
        self.assertIn("if: inputs.validation_mode == 'full'", fingerprint)
        commit = _job(text, "commit")
        self.assertIn("permissions:\n      contents: write", commit)
        self.assertIn("git push origin HEAD:refs/heads/${{ github.ref_name }}", commit)

    def test_main_build_requires_full_evidence_before_main_build(self):
        build = BUILD_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("ou-evidence:", build)
        evidence = _job(build, "ou-evidence", "build")
        self.assertIn("if: github.ref == 'refs/heads/main'", evidence)
        self.assertIn("uses: ./.github/workflows/ou-validation.yml", evidence)
        self.assertIn("validation_mode: full", evidence)
        build_job = _job(build, "build", "release")
        self.assertIn("needs: [ou-evidence, classify]", build_job)
        self.assertIn("needs['ou-evidence'].result == 'success'", build_job)
        self.assertIn("github.ref != 'refs/heads/main'", build_job)
        self.assertIn("refs/heads/main", build_job)

    def test_full_evidence_reuses_one_fingerprinted_archive(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        fingerprint = _job(text, "fingerprint", "regenerate")
        regenerate = _job(text, "regenerate", "combine")
        combine = _job(text, "combine", "commit")
        commit = _job(text, "commit")

        # Smoke has its own download, but the full path resolves the release
        # only in fingerprint and reuses those exact uploaded bytes thereafter.
        self.assertEqual(fingerprint.count("Download versioned simulation data"), 1)
        self.assertIn("name: ou-fingerprinted-simulation-data", fingerprint)
        for section in (regenerate, combine, commit):
            self.assertIn("Reuse fingerprinted simulation data", section)
            self.assertNotIn("gh release download", section)
        self.assertIn("--simulation-zip sim-data-files.zip", fingerprint)
        self.assertIn("--simulation-zip sim-data-files.zip", commit)

    def test_full_evidence_shards_and_combines_each_study_without_cross_mix(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        regenerate = _job(text, "regenerate", "combine")
        combine = _job(text, "combine", "commit")
        self.assertIn("study: [validation, robustness]", regenerate)
        self.assertIn("shard: [0, 1, 2]", regenerate)
        self.assertIn("--shard-count ${{ env.SHARD_COUNT }}", regenerate)
        self.assertIn("--shard-index ${{ matrix.shard }}", regenerate)
        self.assertIn("name: ou-${{ matrix.study }}-shard-${{ matrix.shard }}", regenerate)
        self.assertIn("study: [validation, robustness]", combine)
        self.assertIn("pattern: ou-${{ matrix.study }}-shard-*", combine)
        self.assertIn("name: ou-${{ matrix.study }}-full", combine)
        self.assertNotIn("pattern: ou-*-shard-*", combine)

    def test_publication_is_verified_before_fingerprint_commit_and_push(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        commit = _job(text, "commit")
        self.assertIn("ref: ${{ github.sha }}", commit)
        self.assertIn("Skip when this commit's evidence is already published", commit)
        self.assertIn("replay_provenance", commit)
        self.assertIn("Place bundles and mirror the manuscript copies", commit)
        self.assertIn("make -C tests/validation evidence-test", commit)
        self.assertIn("--write reports/ou_evidence_fingerprint.json", commit)
        self.assertIn("--check reports/ou_evidence_fingerprint.json", commit)
        self.assertIn("git add reports/results/ou_validation reports/results/ou_robustness", commit)
        self.assertIn("git push origin HEAD:refs/heads/${{ github.ref_name }}", commit)
        self.assertLess(commit.index("make -C tests/validation evidence-test"),
                        commit.index("--write reports/ou_evidence_fingerprint.json"))
        self.assertLess(commit.index("--check reports/ou_evidence_fingerprint.json"),
                        commit.index("git add reports/results/ou_validation"))

    def test_evidence_gate_and_proof_gate_are_disjoint_and_complete(self):
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
        self.assertIn("evidence-test: evidence-contract", makefile)
        self.assertIn("proof-test: build", makefile)
        self.assertIn("test: evidence-contract", makefile)
        self.assertIn("python3 -m unittest discover -v -p 'test_*.py'", makefile)

        result = subprocess.run(
            ["make", "--no-print-directory", "-s", "print-test-partition"],
            cwd=REPO_ROOT / "tests" / "validation",
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        rows = dict(line.split("=", 1) for line in result.stdout.splitlines())
        self.assertEqual(set(rows), {"EVIDENCE", "PROOF"})
        groups = {key: value.split() for key, value in rows.items()}
        all_modules = {p.stem for p in (REPO_ROOT / "tests" / "validation").glob("test_*.py")}
        expected_proof = {
            p.stem
            for prefix in ("p2", "p3", "p4", "p5", "sea3")
            for p in (REPO_ROOT / "tests" / "validation").glob(f"test_ou3_{prefix}_*.py")
        }
        self.assertTrue(expected_proof)
        self.assertEqual(set(groups["PROOF"]), expected_proof)
        self.assertEqual(set(groups["EVIDENCE"]), all_modules - expected_proof)
        self.assertFalse(set(groups["PROOF"]) & set(groups["EVIDENCE"]))
        for modules in groups.values():
            self.assertEqual(len(modules), len(set(modules)))

    def test_canonical_proof_workflow_owns_all_expensive_proof_modules(self):
        proof = PROOF_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("make -C tests/validation proof-test", proof)
        self.assertIn("PROOF_TEST_SHARD_INDEX=${{ matrix.shard }}", proof)
        self.assertIn("PROOF_TEST_SHARD_COUNT=6", proof)
        self.assertIn("shard: [0, 1, 2, 3, 4, 5]", proof)
        self.assertIn('"tests/validation/Makefile"', proof)

    def test_build_keeps_current_matrix_and_runs_real_tests_and_documents(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        build = _job(text, "build", "release")
        self.assertIn(
            "dir: [wave_sim, freq, spectrum, ahrs, pii_observer, kalman_ou_ii, "
            "kalman_ou_iii, kalman_tfg, imu_calibrate, detrend, nlo]",
            build,
        )
        self.assertIn("make build", build)
        self.assertIn("./run_tests.sh", build)
        self.assertIn("tools/ou_sim_table.py", build)
        self.assertIn("tools/tfg_sim_table.py", build)
        self.assertIn("Compile LaTeX document", build)
        self.assertIn("Upload PDF artifacts", build)

    def test_build_proof_only_classifier_keeps_complete_sea3_paths_focused(self):
        text = BUILD_WORKFLOW.read_text(encoding="utf-8")
        classify = _job(text, "classify", "ou-evidence")
        self.assertIn("git diff --name-only \"$BEFORE\" \"$AFTER\"", classify)
        self.assertIn(".github/workflows/ou3-complete-sea3.yml", classify)
        self.assertIn("tools/stability/ou3_*.py", classify)
        self.assertIn("tests/validation/test_ou3_sea3_*.py", classify)
        self.assertIn('echo "run_build=false"', classify)
        self.assertIn("refs/heads/main", classify)

    def test_validation_classifier_skips_proof_only_simulation_smoke(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        classify = _job(text, "classify", "proof_only_skip")
        self.assertIn("git diff --name-only \"$BEFORE\" \"$AFTER\"", classify)
        self.assertIn("tools/stability/ou3_*.py", classify)
        self.assertIn('echo "run_smoke=false"', classify)
        skip = _job(text, "proof_only_skip", "validate")
        self.assertIn("no simulation smoke needed", skip)
        self.assertIn("focused proof tests cover it", skip)
        self.assertIn("cancel-in-progress: true", text)

    def test_tfg_validation_does_not_mutate_ou_evidence_or_restamp_fingerprints(self):
        text = TFG_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("Build TFG tests", text)
        self.assertIn("Run TFG mathematical unit tests", text)
        self.assertIn("Checkout exact pre-PR baseline", text)
        self.assertIn("Run SpectralMSE TFG arm", text)
        self.assertIn("Run pre-PR main baseline arm", text)
        self.assertNotIn("ou_replay_fingerprint.py --write", text)
        self.assertNotIn("git push", text)
        self.assertNotIn("reports/results/ou_validation", text)
        self.assertNotIn("reports/results/ou_robustness", text)

    def test_root_make_all_fetches_sim_archive_only_when_missing(self):
        text = ROOT_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn("all: build test", text)
        self.assertIn("test: ensure-sim-data", text)
        self.assertIn("fetch-sim-data:", text)
        self.assertIn('curl -fL "$(SIM_DATA_URL)" -o "$(SIM_DATA_ZIP)"', text)
        self.assertIn("ensure-sim-data:", text)
        self.assertIn('if [ ! -f "$(SIM_DATA_CHECK_FILE)" ]; then', text)
        self.assertIn('$(MAKE) -C "$(REPO_ROOT)" fetch-sim-data', text)
        self.assertEqual(text.count('curl -fL "$(SIM_DATA_URL)"'), 1)

    def test_validation_makefile_materializes_and_checks_evidence_without_hash_editing(self):
        makefile = VALIDATION_MAKEFILE.read_text(encoding="utf-8")
        self.assertIn("publication-sync: build", makefile)
        self.assertIn("ou_publication_sync.py", makefile)
        self.assertIn("evidence-contract: publication-sync", makefile)
        self.assertIn("ou_evidence_contract.py --auto", makefile)
        self.assertNotIn("sed -i", makefile)
        self.assertNotIn("reports/ou_evidence_fingerprint.json", makefile)


if __name__ == "__main__":
    unittest.main()
