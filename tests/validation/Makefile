.PHONY: build publication-sync evidence-contract test evidence-test proof-test print-test-partition clean

# The staged OU-III proof searches. Each of these modules is minutes of
# interval arithmetic and the set is around three quarters of an hour, which is
# why ou3-proof.yml runs them in jobs that budget one to two hours each. They
# read source and tooling only: no bundle this repository publishes can change
# their verdict. `evidence-test` therefore leaves them to the workflow that
# owns them, and runs every contract that does depend on the evidence tree.
PROOF_SEARCH_TESTS := $(wildcard test_ou3_p2_*.py test_ou3_p3_*.py test_ou3_p4_*.py test_ou3_p5_*.py test_ou3_sea3_*.py)
PROOF_TEST_MODULES := $(patsubst %.py,%,$(sort $(PROOF_SEARCH_TESTS)))
PROOF_TEST_SHARD_INDEX ?= 0
PROOF_TEST_SHARD_COUNT ?= 1
EVIDENCE_TESTS := $(filter-out $(PROOF_SEARCH_TESTS),$(wildcard test_*.py))
EVIDENCE_TEST_MODULES := $(EVIDENCE_TESTS:.py=)

build:
	python3 -m py_compile ../../tools/*.py ../../tools/stability/*.py test_*.py

publication-sync: build
	python3 ../../tools/ou_publication_sync.py \
		--validation-dir ../../reports/results/ou_validation \
		--direction-results ../../doc/kalman_ou_iii/w3d-wave-direction-results.tex-part \
		--deterministic-results ../../doc/kalman_ou_iii/w3d-sim-results-generated.tex-part

evidence-contract: publication-sync
	python3 ../../tools/ou_evidence_contract.py --auto

test: evidence-contract
	python3 -m unittest discover -v -p 'test_*.py'

evidence-test: evidence-contract
	python3 -m unittest -v $(EVIDENCE_TESTS:.py=)

# The canonical proof workflow runs every source/proof test, in disjoint
# shards. No proof threshold or test is removed from the complete `test` target.
proof-test: build
	PROOF_TEST_SHARD_INDEX=$(PROOF_TEST_SHARD_INDEX) PROOF_TEST_SHARD_COUNT=$(PROOF_TEST_SHARD_COUNT) python3 -c 'import os, sys, unittest; i = int(os.environ["PROOF_TEST_SHARD_INDEX"]); n = int(os.environ["PROOF_TEST_SHARD_COUNT"]); modules = sys.argv[1:]; assert n > 0 and 0 <= i < n, "invalid proof shard"; selected = modules[i::n]; assert selected, "empty proof shard"; print("PROOF_TEST_MODULES", selected, flush=True); result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromNames(selected)); sys.exit(not result.wasSuccessful())' $(PROOF_TEST_MODULES)

print-test-partition:
	@echo EVIDENCE=$(EVIDENCE_TEST_MODULES)
	@echo PROOF=$(PROOF_TEST_MODULES)

clean:
	rm -rf __pycache__ ../../tools/__pycache__ ../../tools/stability/__pycache__
