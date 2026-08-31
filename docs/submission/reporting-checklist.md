# Submission and reporting checklist

Target: *Nature Communications*, Article format. This is a working audit, not a journal-issued
checklist. The journal's current formatting instructions should be rechecked on the submission
date: <https://www.nature.com/ncomms/submit/article>.

## Manuscript structure

| Item | Status | Evidence or action |
| --- | --- | --- |
| Concise title | Complete | 9 words; no punctuation-dependent subtitle |
| Unstructured abstract, no references | Complete | 191 words |
| Introduction, Results, Discussion, Methods | Complete | Present in submission manuscript |
| Main text within target length | Complete | Under 5,000 words excluding Methods, references and end matter |
| Data Availability | Complete | Public accession/DOI statements included |
| Code Availability | Complete | Public repository, `v1.0.0`, licenses and Zenodo DOI stated |
| References | Complete for draft | 25 primary/source references; bibliographic manager validation still recommended |
| Acknowledgements | Complete | No external funding; Metastate compute resources disclosed |
| Author Contributions | Complete | Sole-author contribution and CRediT roles stated |
| Competing Interests | Complete | Founder, affiliation and commercial-benefit interest disclosed |
| Ethics statement | Not applicable to reanalysis as drafted | Confirm source-study consent and repository-use terms before submission |

## Reproducibility and statistics

| Item | Status | Evidence or action |
| --- | --- | --- |
| Analysis chronology disclosed | Complete | Frozen, adaptive, locked and post-result roles named |
| Biological holdout units disclosed | Complete | Participant, population category and lineage |
| Leakage controls disclosed | Complete | Training-only imputation, scaling and fitting |
| Null-model construction disclosed | Complete | Family-size, fixed-degree and dimension-matched nulls |
| Sample sizes stated | Complete | Main text and Supplementary Table 1 |
| Statistical units stated | Complete | Supplementary Methods |
| Confidence-interval methods stated | Complete | Paired group or target bootstrap |
| Multiplicity stated for selected interaction | Complete | FDR screen described |
| Complete sensitivity results | Complete | 18 human rows and 10 CCLE settings supplied |
| Missingness/mapping exclusions | Complete | Five machine-readable mapping and feature-set tables |
| Random seeds/configurations available | Complete in repository | Frozen YAML and deterministic implementations |
| Source/output checksums | Complete | Artifact manifests and publication synthesis |
| Causal/clinical/flux boundary | Complete | Abstract, Methods, Discussion and supplement |

## Data and code release tasks

| Task | Status |
| --- | --- |
| Create clean public source repository | Complete |
| Select and add open-source license | Complete: Apache-2.0 code; CC BY 4.0 research outputs |
| Remove private paths, credentials and non-redistributable raw data | Complete |
| Tag immutable release matching manuscript | Complete: `v1.0.0` |
| Archive release and obtain DOI | Complete: 10.5281/zenodo.22207315 |
| Insert repository URL, release tag and DOI into Code Availability | Complete |
| Validate all source URLs and checksums from a clean environment | Complete: 83/83 declared outputs byte-identical |
| Deposit Supplementary Data files in journal-supported archive if required | Pending submission workflow |

## Author and editorial declarations

| Item | Status |
| --- | --- |
| Final author list, affiliations and correspondence | Complete |
| All-author approval | Sole author; final submission approval remains with O.Ü. |
| Funding and grant numbers | Complete: no external funding declared |
| Competing-interest and Metastate/IP statement | Complete |
| Related manuscripts/preprints disclosed | Pending author confirmation |
| Suggested/opposed reviewers with conflict checks | Shortlist prepared; author approval pending |
| Data-use and source-license review | Complete for public release; journal-level confirmation still advised |
| Cover-letter originality statement | Drafted; author confirmation required |

## Figures and files

| Item | Status | Action |
| --- | --- | --- |
| Main manuscript DOCX/PDF | Complete | Final exports inspected at A4 page size; title and figure-caption rendering corrected |
| Supplement DOCX/PDF | Complete | Final landscape export inspected; wide tables and symbols are readable |
| Main figure set | Complete | Study design, human validation and CCLE validation in 300-dpi PNG and vector PDF |
| Editable/source figures | Complete | Deterministic matplotlib source plus machine-readable source index |
| Figure legends | Complete | Three main legends and one supplementary legend present |
| Accessibility/color review | Pending | Check contrast, grayscale and readable labels |

## Submission blockers

The scientific draft, declarations and reproducibility release are complete. Remaining editorial
gates are independent domain review, accessibility review, final reference-manager validation and
the author's explicit approval of the journal submission.
