PYTHON ?= python3
VENV := .venv
BIN := $(VENV)/bin

.PHONY: env verify test fetch reproduce

env:
	$(PYTHON) -m venv $(VENV)
	$(BIN)/pip install -r requirements.lock
	$(BIN)/pip install --no-deps -e .

verify:
	$(BIN)/python scripts/verify_release.py

test:
	$(BIN)/pytest -q

fetch:
	$(BIN)/python scripts/fetch_st002081.py
	$(BIN)/factorized-fetch-human-gem
	$(BIN)/factorized-fetch-ccle

reproduce: fetch
	$(BIN)/factorized-ccle-catalog --config config/ccle-expanded.yaml --output-dir artifacts/ccle-expanded --skip-fetch
	$(BIN)/factorized-st002081
	$(BIN)/factorized-st002081-structural
	$(BIN)/factorized-st002081-null-sensitivity
	$(BIN)/factorized-st000818-replication
	$(BIN)/factorized-ccle
	$(BIN)/factorized-human-sensitivity
	$(BIN)/factorized-human-graph-mixing
	$(BIN)/factorized-ccle-sensitivity
	$(BIN)/factorized-ccle-null-ensemble
	$(BIN)/factorized-ccle-null-ensemble --config config/pathway-score-ccle-property-matched-null.yaml --output artifacts/pathway-score-ccle-property-matched-null
	$(BIN)/factorized-ccle-reaction-cluster
	$(BIN)/factorized-mapping-supplement
	$(BIN)/factorized-publication-figures
	$(BIN)/factorized-publication-synthesis
	$(BIN)/python scripts/verify_release.py

