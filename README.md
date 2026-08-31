# Matched nulls for metabolic pathway scores

This repository is the reproducible code and evidence release for **“Matched nulls reveal
specificity limits of metabolic pathway scores.”** It tests a simple qualification rule:
a biological representation must predict outside the people, populations, or cancer lineages
used for fitting and outperform a random representation matched for effective complexity.

The release contains no private or participant-level source data. It contains frozen protocols,
analysis code, exact configurations, checksummed aggregate outputs, complete sensitivity grids,
accepted/rejected mapping ledgers, figures, and the manuscript package.

- **Release:** [`v1.0.1`](https://github.com/Metastatebio/factorized-pathway-scores/releases/tag/v1.0.1)
- **Archive DOI:** [10.5281/zenodo.22207315](https://doi.org/10.5281/zenodo.22207315)

## Main results

- In 1,539 repeated samples from 112 people, coarse lipid families reconstructed hidden markers
  but failed against size-matched random groups.
- A 95-descriptor lipid-structure representation beat 20 fixed-degree graph nulls and replicated
  in a separately locked cohort of 450 people from 15 population categories.
- In 913 CCLE cell lines, direct HumanGEM/GPR features predicted 60 held-out metabolites across
  unseen lineages and beat 20 dimension-matched and 20 network-degree-and-coverage-matched random
  feature draws.
- Transcript interactions did not improve the general direct-network model.

These results qualify representations for hidden-measurement reconstruction. They do **not**
establish physiological flux, causality, treatment response, clinical utility, or same-person
genomic personalization.

## Start here

- [Manuscript](docs/submission/factorized-pathway-scores-nature-communications.pdf)
- [Supplement](docs/submission/factorized-pathway-scores-supplement.pdf)
- [Frozen benchmark protocol](docs/protocol-factorized-pathway-score-benchmark.md)
- [Machine-readable claim ledger](artifacts/factorized-pathway-publication/claim-ledger.csv)
- [Source-data and license audit](docs/source-data.md)
- [Reproduction instructions](docs/reproducibility.md)

## Verify the released evidence

Python 3.12.3 is the reference interpreter. The exact package set is frozen in
`requirements.lock`.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install --no-deps -e .
.venv/bin/python scripts/verify_release.py
.venv/bin/pytest -q
```

The verification command checks every declared output checksum and the principal result
invariants without downloading source data. See `make reproduce` for a full rerun from public
inputs.

## Repository structure

```text
artifacts/   Checksummed aggregate outputs, sensitivity grids, mappings, and figures
config/      Frozen analysis configurations and random seeds
docs/        Protocols, manuscript, supplement, and source-data audit
scripts/     Input retrieval and release-integrity checks
src/         Analysis implementation
tests/       Unit and pipeline-contract tests
```

## License and citation

Code is Apache-2.0. Manuscripts, protocols, figures, and aggregate result tables are CC BY 4.0.
Third-party raw data are excluded and remain under their source terms. Cite the archived release
using [DOI 10.5281/zenodo.22207315](https://doi.org/10.5281/zenodo.22207315) and
`CITATION.cff`; cite each upstream dataset as listed in `docs/source-data.md`.
