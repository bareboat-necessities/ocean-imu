.PHONY: build test clean

build:
	python3 -m py_compile ../../tools/ou_validation.py test_ou_validation.py

test: build
	python3 -m unittest -v test_ou_validation.py

clean:
	rm -rf __pycache__ ../../tools/__pycache__
