.PHONY: build test clean

build:
	python3 -m py_compile \
		../../tools/model_mismatch_ablation.py \
		test_model_mismatch_ablation.py

test: build
	python3 -m unittest -v test_model_mismatch_ablation.py

clean:
	rm -rf __pycache__ ../../tools/__pycache__
