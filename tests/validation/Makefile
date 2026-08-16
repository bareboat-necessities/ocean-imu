.PHONY: build evidence-contract test clean

build:
	python3 -m py_compile \
		../../tools/ou_validation.py \
		../../tools/ou_robustness.py \
		../../tools/ou_evidence_contract.py \
		../../tools/rotate_record_heading.py \
		test_ou_validation.py \
		test_ou_robustness.py \
		test_ou_evidence_contract.py \
		test_record_conventions.py \
		test_publication_references.py \
		test_workflow_contract.py \
		test_reorg_compat_overrides.py \
		test_zz_publication_contract.py

evidence-contract: build
	python3 ../../tools/ou_evidence_contract.py --auto

test: evidence-contract
	# The publication modules intentionally replace obsolete prose/layout-pinned
	# manuscript assertions in test_ou_validation. Import them before discovery
	# so unittest constructs the suite from the semantic replacements.
	python3 -c 'import unittest, test_zz_publication_contract, test_reorg_compat_overrides; suite=unittest.defaultTestLoader.discover(".", pattern="test_*.py"); result=unittest.TextTestRunner(verbosity=2).run(suite); raise SystemExit(0 if result.wasSuccessful() else 1)'

clean:
	rm -rf __pycache__ ../../tools/__pycache__
