# Source-data and redistribution audit

No participant-level or third-party raw data are distributed in this repository. Retrieval scripts
write to the ignored `data/raw/` directory and fail closed on the checksums used by the analyses.
The repository distributes only original code, protocols, aggregate analysis outputs, mapping
ledgers, and figures.

| Source | Locked object | Primary location | Integrity | Redistribution decision |
| --- | --- | --- | --- | --- |
| Metabolomics Workbench ST002081 / AN003790 | mwTab analysis file | `https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN003790/mwtab/txt` | Trailing-newline-normalized SHA-256 `e5d68bea4f6adf113f9bafa0e8bd333b949f289f0dc645a18435d542c3ae9c9b` | Raw file excluded; users retrieve from the repository |
| Metabolomics Workbench ST000818 / AN001299 | Factors and processed measurements | Workbench REST endpoints under `https://www.metabolomicsworkbench.org/rest` | SHA-256 recorded in the released analysis manifest | Raw responses excluded; users retrieve from the repository |
| CCLE 2019 | Metabolomics, RNA expression, cell-line annotations | `https://data.broadinstitute.org/ccle` | Three SHA-256 values hard-coded in `datasets/ccle.py` and checked before use | Raw files excluded; cite Li et al. 2019 and comply with Broad/source terms |
| Human-GEM v2.0.0 | Model YAML, genes and metabolites | `https://github.com/SysBioChalmers/Human-GEM/tree/v2.0.0/model` | Three SHA-256 values hard-coded in `fetch_human_gem.py` | Raw files excluded; upstream model is CC BY 4.0 and must be attributed |

Metabolomics Workbench exposes the cited studies for programmatic download through its documented
REST service. The absence of raw redistribution here is deliberately conservative: public access
is not treated as permission to relicense deposited participant-level data. Human-GEM identifies
its model repository as CC BY 4.0; this release cites the project and retains source identifiers.

Derived tables can contain assay labels, stable metabolite identifiers, reaction identifiers, and
aggregate performance statistics. They do not contain direct participant identifiers or source
abundance matrices. The project licenses do not override source-dataset terms, privacy rights,
trademark rights, or citation requirements.

Primary citations and dataset DOIs are listed in the manuscript Data Availability statement.

