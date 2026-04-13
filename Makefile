PROJECT := earthkit-geo
CONDA := conda
CONDAFLAGS :=
COV_REPORT := html

setup:
	pre-commit install

default: qa unit-tests

qa:
	pre-commit run --all-files

unit-tests:
	python -m pytest -vv --cov=. --cov-report=$(COV_REPORT)
