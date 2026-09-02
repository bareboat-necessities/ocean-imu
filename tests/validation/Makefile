.PHONY: build publication-sync evidence-contract test clean

build:
	python3 -m py_compile ../../tools/*.py test_*.py

publication-sync: build
	python3 ../../tools/ou_publication_sync.py \
		--validation-dir ../../reports/results/ou_validation \
		--direction-results ../../doc/kalman_ou_iii/w3d-wave-direction-results.tex-part \
		--deterministic-results ../../doc/kalman_ou_iii/w3d-sim-results-generated.tex-part

evidence-contract: publication-sync
	python3 ../../tools/ou_evidence_contract.py --auto

test: evidence-contract
	python3 -m unittest discover -v -p 'test_*.py'

clean:
	rm -rf __pycache__ ../../tools/__pycache__
