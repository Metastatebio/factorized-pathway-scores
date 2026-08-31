# Factorized pathway-score publication synthesis

**Decision:** `ROBUST_HUMAN_STRUCTURAL_AND_REACTION_TOPOLOGY_SIGNAL`

## Evidence domains

| domain | dataset | role | result | status |
| --- | --- | --- | --- | --- |
| Repeated human metabolomics | ST002081 | Frozen coarse benchmark + adaptive structural follow-up | Coarse R2=0.588 failed its null; adaptive structural R2=0.907 passed across 20 graph nulls | ADAPTIVE_SIGNAL |
| External human structural replication | ST000818 | Locked population-held-out replication | 255 lipids across 15 groups; passed 20 graph nulls and 9 settings; minimum edge replacement=0.771 | SUPPORTED |
| Reaction topology prediction | CCLE + HumanGEM | New lineage-held-out metabolite benchmark | Mechanistic RMSE=0.797; improvement vs random=0.124; 10 settings, 20 dimension-matched and 20 property-matched random seeds pass; cluster lower bound=0.086 | SUPPORTED |
| Mapping and feature provenance | ST002081 + ST000818 + CCLE | Accepted/rejected feature and model-input audit | 5 machine-readable tables; all source integrity verified | SUPPORTED |
| Transcript-conditioned coupling | CCLE | Prior reaction-resolved inference | 1 FDR-significant state among 284 tests; no general interaction increment | CANDIDATE_SPECIFIC |
| Human genetic topology | METSIM + HumanGEM | Prior degree-matched topology benchmark | 26/71 direct recoveries; P=1e-05 | SUPPORTIVE_NOT_INDEPENDENT |
| Genotype-by-context reaction coordinate | UKB/BBJ summaries | Prior reciprocal summary-statistic scan | 1 strict replicated row; 0 novel mutable-context rows | LIMITED |
| Isotope-resolved external challenge | NSCLC | Prior falsification | 0 expression and 0 oncogenotype FDR signals | NULL |
| Constraint sensitivity | HumanGEM | Prior 1,000-draw model ensemble | 3 stable and 2 unstable configured candidates | MODEL_RELATIVE |
| Directional longitudinal prediction | ST002081 | Prior temporal falsification | forward R2=0.452; reverse R2=0.459 | NOT_SUPPORTED |

## Claim ledger

| claim | decision | basis |
| --- | --- | --- |
| Compact lipid-family scores reconstruct hidden human markers | SUPPORTED_VS_POPULATION | Participant-bootstrap RMSE interval above zero versus population mean |
| Broad lipid-family scores are pathway-specific | NOT_SUPPORTED | Random-group RMSE gate failed |
| Structure-resolved lipid scores beat a degree-preserving null | SUPPORTED_ADAPTIVE_SAME_COHORT | Adaptive participant-paired RMSE and precision intervals above zero |
| Structure-resolved lipid scores replicate in a separate human cohort | SUPPORTED_EXTERNAL_POPULATION_HELD_OUT | ST000818 population-held-out validation across 20 graph nulls |
| Direct HumanGEM/GPR features outperform size-matched random features | SUPPORTED | Lineage-isolated target-bootstrap comparison across CCLE metabolites |
| Transcript interactions improve the general direct-network model | NOT_SUPPORTED | Interaction model did not beat network-metabolite or additive models |
| Transcript-conditioned coupling can exist for selected reactions | SUPPORTED_IN_MODEL_SYSTEM | FDR-significant CDA–cytidine–uridine interaction with robustness checks |
| The evidence validates same-person genomic personalization | NOT_TESTED | No available cohort contains the required participant-level modalities |
| The model estimates physiological flux or treatment response | NOT_SUPPORTED | Constraint results are model-relative and external drug/isotope tests are null |
| The human structural result is robust across declared analysis settings | SUPPORTED_POST_RESULT | All nine settings pass both endpoints in both human cohorts |
| Human fixed-degree null graphs are moved and diverse | SUPPORTED_POST_RESULT | Minimum edge replacement=0.771; maximum pairwise Jaccard=0.131 |
| The CCLE reaction-topology result is robust across declared settings | SUPPORTED_POST_RESULT | Matched random-control lower bounds pass in all ten settings |
| Transcript interactions become generally useful under alternative settings | NOT_SUPPORTED | Interactions fail to beat the additive model in all ten settings |
| The CCLE topology result is stable across random feature draws | SUPPORTED_POST_RESULT | All 20 null seeds pass; expected-null lower bound=0.093 SD |
| CCLE topology beats network-degree and coverage-matched features | SUPPORTED_POST_RESULT | All 20 hard-null seeds pass; expected-null lower bound=0.097 SD |
| The CCLE topology effect is not concentrated in one subsystem | SUPPORTED_POST_RESULT | 20-cluster bootstrap lower bound=0.086 SD; leave-one-cluster gate passes |
| Analysis-facing mappings and exclusions are auditable | SUPPORTED_RESOURCE | Accepted/rejected human and CCLE mapping tables with verified source integrity |

## Interpretation

Direct metabolic-network neighborhoods carry reproducible predictive information in lineage-held-out cancer-cell-line data, but adding transcript interactions does not improve the general predictor. Broad lipid-family compression recovers substantial human profile information but does not outperform size-matched random compression. An explicitly adaptive structure-resolved lipid representation does beat degree-preserving random graphs in the same cohort, and this result replicates in a separate human cohort with complete population categories held out. Together, the results support empirical calibration of pathway-score resolution and candidate-specific interaction testing, not universal multi-omic fusion.

## Boundary

This synthesis benchmarks evidence components across different public cohorts and model systems. It is not same-participant multi-omics, does not infer missing DNA or RNA as measured, and does not establish physiological flux, prospective clinical prediction, diagnosis, treatment selection or drug efficacy.
