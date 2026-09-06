.PHONY: build publication-sync evidence-contract test evidence-test clean

# The staged OU-III proof searches. Each of these modules is minutes of
# interval arithmetic and the set is around three quarters of an hour, which is
# why ou3-proof.yml runs them in jobs that budget one to two hours each. They
# read source and tooling only: no bundle this repository publishes can change
# their verdict. `evidence-test` therefore leaves them to the workflow that
# owns them, and runs every contract that does depend on the evidence tree.
PROOF_SEARCH_TESTS := $(wildcard test_ou3_p2_*.py test_ou3_p3_*.py test_ou3_p4_*.py test_ou3_p5_*.py)
EVIDENCE_TESTS := $(filter-out $(PROOF_SEARCH_TESTS),$(wildcard test_*.py))

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

clean:
	rm -rf __pycache__ ../../tools/__pycache__ ../../tools/stability/__pycache__
