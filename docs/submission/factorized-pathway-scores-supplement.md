# Supplementary information

## Matched nulls reveal specificity limits of metabolic pathway scores

This supplement accompanies the main manuscript. It distinguishes frozen primary analyses,
adaptive analyses, locked external replication, post-result robustness analyses and supporting
evidence. Machine-readable tables remain authoritative where rounded values are shown here.

## Supplementary Methods

### Evidence roles and analysis chronology

The ST002081 coarse-family benchmark was frozen before execution. Its structural-descriptor
follow-up was specified after the coarse benchmark failed its matched-random specificity gate and
is therefore adaptive. ST000818 was selected from a pre-existing public-study registry using only
metadata and feature names; its population-held-out protocol was frozen before the abundance
matrix was retrieved. The human and CCLE parameter grids were declared after their respective
primary results and are robustness analyses, not prospective confirmation. Earlier genetics,
constraint, isotope, temporal and interaction analyses are supporting domains and were not pooled
with the primary observations.

### Matched-null rationale

The no-information population-mean baseline tests whether a representation predicts at all. It
does not test whether the declared biological labels matter. The coarse benchmark therefore used
random groups preserving family sizes. The lipid-structure benchmark randomized the binary
feature-by-descriptor graph while preserving every feature degree and every descriptor degree.
The CCLE benchmark used random metabolite, GPR and interaction feature sets matched to each
target's feature counts and evaluated with the same lineage folds. These nulls address effective
representation complexity; they do not prove causal mechanism.

### Group-isolated validation

All samples from a ST002081 participant were held out together. All samples from a ST000818
population category were held out together. All cell lines from a CCLE lineage were held out
together in two repeats of five deterministic size-balanced folds. Imputation, scaling, feature
summaries and ridge models were fitted using training observations only.

### Statistical units

Human confidence intervals used paired resampling of the outer biological groups: participants in
ST002081 and population categories in ST000818. The primary analyses used 5,000 bootstrap draws;
the post-result grid used 2,000 draws per setting. CCLE comparisons used 5,000 paired bootstrap
draws over prediction targets. The stochastic CCLE null ensemble used 2,000 target-bootstrap draws
per random seed and 5,000 for the expected-null contrast. These intervals quantify variation across the declared resampling
units. They do not make correlated metabolites or reaction neighborhoods independent.

### Mapping and source integrity

The feature maps retain every accepted and rejected row and a reason for exclusion. CCLE mappings
record assay names, stable identifiers, HumanGEM metabolites and compartments, GPR eligibility,
and the exact feature set used for each target. Source and output manifests contain SHA-256
checksums. The publication synthesis verified all 15 contributing artifact manifests.

## Supplementary Table 1 | Cohorts and benchmark roles

| Resource | Analysis unit | Retained data | Holdout unit | Evidential role |
| --- | --- | ---: | --- | --- |
| ST002081 | Longitudinal human plasma lipidomics | 1,539 samples; 112 people; 493 eligible lipids | Participant | Frozen coarse benchmark; adaptive structural analysis |
| ST000818 / AN001299 | Human blood lipidomics | 450 samples; 15 population categories; 255 eligible lipids | Population category | Protocol-locked external structural replication |
| CCLE | Cancer-cell-line metabolomics and RNA | 913 cell lines; 76 mapped assay metabolites; 60 primary targets | Lineage | Reaction-neighborhood prediction and transcript ablation |
| METSIM summary results | Human mQTL states | 71 eligible gene–metabolite states | Not re-fitted as individual data | Positive-control topology audit |
| NSCLC public data | Cell-line abundance, isotope and molecular state | Source-specific eligible observations | Source-defined groups | External falsification of stronger mechanistic claims |
| HumanGEM / Human1 | Curated human metabolic network | Reaction, metabolite and GPR records | Not applicable | Topology and model-relative constraint evidence |

## Supplementary Table 2 | Primary benchmark contrasts

| Benchmark | Biological representation | Matched null | Out-of-group unit | Primary inference |
| --- | --- | --- | --- | --- |
| ST002081 coarse | Lipid-family median scores | Random groups preserving family sizes | Participant | Predictive compression, but no family specificity |
| ST002081 structural | Lipid class, carbon, unsaturation and chain descriptors | Fixed feature and descriptor degree graph | Participant | Adaptive structural specificity |
| ST000818 structural | Same parser and representation contract | Fixed feature and descriptor degree graph | Population category | External structural replication |
| CCLE topology | Direct HumanGEM neighbors plus reaction-linked GPR features | Dimension-matched random metabolites, GPRs and products | Cancer lineage | Reaction topology has predictive specificity |
| CCLE transcript ablation | Direct neighbors plus additive or interaction GPR terms | Nested direct-neighbor/additive models | Cancer lineage | Additive context helps compact prediction; interactions do not improve it generally |

## Supplementary Table 3 | Human structural sensitivity grid

All entries are structural-minus-degree-null improvements. Parentheses are paired 95% bootstrap
intervals. All randomized graphs preserved exact row and column degrees.

| Cohort | Setting | Descriptors | Structural RMSE | Null RMSE | RMSE improvement (95% CI) | Precision improvement (95% CI) |
| --- | --- | ---: | ---: | ---: | --- | --- |
| ST002081 | Primary | 95 | 0.309 | 0.411 | 0.104 (0.098, 0.109) | 0.133 (0.119, 0.148) |
| ST002081 | Ridge alpha 1 | 95 | 0.311 | 0.408 | 0.102 (0.097, 0.108) | 0.142 (0.129, 0.155) |
| ST002081 | Ridge alpha 100 | 95 | 0.326 | 0.432 | 0.105 (0.099, 0.112) | 0.139 (0.125, 0.153) |
| ST002081 | Descriptor minimum 3 | 112 | 0.290 | 0.375 | 0.089 (0.084, 0.093) | 0.108 (0.095, 0.120) |
| ST002081 | Descriptor minimum 10 | 69 | 0.343 | 0.448 | 0.109 (0.103, 0.116) | 0.179 (0.161, 0.195) |
| ST002081 | Maximum prevalence 0.80 | 95 | 0.309 | 0.407 | 0.101 (0.095, 0.106) | 0.133 (0.119, 0.146) |
| ST002081 | Maximum prevalence 0.95 | 95 | 0.309 | 0.410 | 0.103 (0.097, 0.110) | 0.140 (0.127, 0.153) |
| ST002081 | Five swaps per edge | 95 | 0.309 | 0.404 | 0.097 (0.092, 0.102) | 0.127 (0.113, 0.140) |
| ST002081 | Twenty-five swaps per edge | 95 | 0.309 | 0.399 | 0.093 (0.088, 0.099) | 0.126 (0.114, 0.139) |
| ST000818 | Primary | 81 | 1.237 | 1.790 | 0.563 (0.057, 1.006) | 0.103 (0.086, 0.121) |
| ST000818 | Ridge alpha 1 | 81 | 1.288 | 1.788 | 0.511 (0.062, 0.903) | 0.100 (0.084, 0.117) |
| ST000818 | Ridge alpha 100 | 81 | 1.253 | 1.786 | 0.543 (0.041, 0.969) | 0.074 (0.058, 0.097) |
| ST000818 | Descriptor minimum 3 | 112 | 1.192 | 1.682 | 0.497 (0.050, 0.882) | 0.073 (0.054, 0.092) |
| ST000818 | Descriptor minimum 10 | 41 | 1.764 | 1.821 | 0.056 (0.042, 0.120) | 0.097 (0.081, 0.114) |
| ST000818 | Maximum prevalence 0.80 | 81 | 1.237 | 1.767 | 0.541 (0.062, 0.960) | 0.102 (0.085, 0.121) |
| ST000818 | Maximum prevalence 0.95 | 81 | 1.237 | 1.787 | 0.561 (0.068, 0.993) | 0.102 (0.086, 0.119) |
| ST000818 | Five swaps per edge | 81 | 1.237 | 1.794 | 0.569 (0.065, 1.006) | 0.102 (0.082, 0.125) |
| ST000818 | Twenty-five swaps per edge | 81 | 1.237 | 1.773 | 0.547 (0.057, 0.975) | 0.093 (0.071, 0.113) |

Decision: all nine settings passed both endpoints in each cohort. Because this grid was designed
after the primary results, the decision is `HUMAN_STRUCTURAL_SENSITIVITY_ROBUST`, not independent
confirmation.

The graph-quality diagnostic regenerated 200 panel-specific nulls across the two cohorts. All row
and column degrees were exact, every panel contained 20 unique nulls, minimum observed-edge
replacement was 0.771 and maximum pairwise null Jaccard was 0.131. This demonstrates graph movement
and diversity but not perfectly uniform sampling from every graph with the same degrees.

## Supplementary Table 4 | CCLE reaction-topology sensitivity grid

The improvement is random-factorized RMSE minus reaction-factorized RMSE, in training-fold outcome
standard-deviation units. Positive values favor topology. The final column is interaction-model
RMSE improvement over the additive GPR model; negative values favor the additive model.

| Setting | Targets | Signatures | Factorized RMSE | Median target R² | Improvement vs random (95% CI) | Improvement vs additive (95% CI) |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Primary | 60 | 169 | 0.797 | 0.241 | 0.124 (0.084, 0.166) | −0.009 (−0.015, −0.004) |
| Maximum 4 metabolites/reaction | 49 | 75 | 0.864 | 0.193 | 0.084 (0.028, 0.135) | −0.008 (−0.014, −0.002) |
| Maximum 6 metabolites/reaction | 59 | 154 | 0.820 | 0.236 | 0.099 (0.060, 0.140) | −0.009 (−0.015, −0.004) |
| Maximum 12 metabolites/reaction | 60 | 198 | 0.795 | 0.253 | 0.131 (0.088, 0.176) | −0.011 (−0.016, −0.005) |
| Unbounded reaction size | 60 | 202 | 0.795 | 0.253 | 0.125 (0.080, 0.171) | −0.011 (−0.017, −0.006) |
| Include transport | 63 | 192 | 0.709 | 0.469 | 0.219 (0.169, 0.273) | −0.014 (−0.020, −0.009) |
| Ridge alpha 1 | 60 | 169 | 0.799 | 0.237 | 0.124 (0.083, 0.166) | −0.010 (−0.016, −0.005) |
| Ridge alpha 100 | 60 | 169 | 0.794 | 0.252 | 0.117 (0.079, 0.157) | −0.006 (−0.010, −0.002) |
| GPR mean | 60 | 169 | 0.797 | 0.241 | 0.123 (0.086, 0.166) | −0.009 (−0.014, −0.004) |
| GPR sum | 60 | 169 | 0.797 | 0.241 | 0.123 (0.083, 0.167) | −0.009 (−0.014, −0.004) |

Decision: reaction topology beat its matched random representation in all ten settings. Transcript
interactions failed to beat additive GPR terms in all ten settings. The grid is post-result
robustness evidence.

## Supplementary Table 5 | CCLE stochastic random-feature null ensemble

The biological feature set and lineage folds were fixed. Each row independently regenerated the
dimension-matched random features. Effects are random-minus-topology equal-lineage RMSE.

| Random seed | Random RMSE | Improvement (95% target-bootstrap CI) |
| ---: | ---: | --- |
| 20262000 | 0.927 | 0.130 (0.082, 0.181) |
| 20262001 | 0.924 | 0.127 (0.085, 0.171) |
| 20262002 | 0.902 | 0.105 (0.055, 0.157) |
| 20262003 | 0.925 | 0.128 (0.080, 0.172) |
| 20262004 | 0.936 | 0.139 (0.098, 0.185) |
| 20262005 | 0.921 | 0.124 (0.071, 0.176) |
| 20262006 | 0.939 | 0.142 (0.091, 0.195) |
| 20262007 | 0.947 | 0.150 (0.106, 0.195) |
| 20262008 | 0.922 | 0.126 (0.074, 0.173) |
| 20262009 | 0.919 | 0.122 (0.080, 0.168) |
| 20262010 | 0.930 | 0.134 (0.087, 0.184) |
| 20262011 | 0.924 | 0.127 (0.085, 0.177) |
| 20262012 | 0.939 | 0.143 (0.102, 0.189) |
| 20262013 | 0.920 | 0.123 (0.081, 0.168) |
| 20262014 | 0.926 | 0.129 (0.082, 0.179) |
| 20262015 | 0.931 | 0.134 (0.091, 0.182) |
| 20262016 | 0.951 | 0.154 (0.092, 0.232) |
| 20262017 | 0.970 | 0.173 (0.111, 0.240) |
| 20262018 | 0.925 | 0.128 (0.076, 0.183) |
| 20262019 | 0.941 | 0.144 (0.099, 0.196) |
| Expected null over 20 seeds | 0.931 | 0.134 (0.093, 0.177) |

All 20 per-seed intervals excluded zero. The biological-factorized RMSE remained exactly 0.797 SD
across seeds. This post-result ensemble addresses stochastic null selection, not external cohort
generalization.

## Supplementary Table 6 | Property-matched null and subsystem dependence audit

The hard null matched random metabolites on direct-network degree and assay coverage and random
GPR signatures on candidate-network degree. It did not use target values or target correlations.

| Audit | Units | Effect or balance | 95% interval / gate |
| --- | ---: | ---: | --- |
| Property-matched expected null | 20 seeds; 60 targets | 0.137 SD | 0.097–0.181 |
| Property-matched seed range | 20 seeds | 0.116–0.149 SD | Every per-seed lower bound >0 |
| Metabolite log-degree imbalance | 20 seeds | 0.055 mean absolute difference | Threshold ≤0.35 |
| Metabolite assay-coverage imbalance | 20 seeds | 0.000 mean absolute difference | Threshold ≤0.02 |
| GPR log-degree imbalance | 20 seeds | 0.009 mean absolute difference | Threshold ≤0.35 |
| HumanGEM subsystem-cluster bootstrap | 20 clusters; 60 targets | 0.137 SD | 0.086–0.194 |
| Leave-one-subsystem-out | 20 exclusions | Minimum 0.120 SD | Every effect >0 |

The subsystem audit assigns each target to its most frequent eligible HumanGEM subsystem and is a
dependence sensitivity, not proof that pathways or targets are independent.

## Supplementary Table 7 | Mapping and eligibility audit

| Domain | Total | Accepted | Rejected | Acceptance fraction |
| --- | ---: | ---: | ---: | ---: |
| ST000818 lipid features | 260 | 255 | 5 | 98.1% |
| ST002081 lipid features | 845 | 493 | 352 | 58.3% |
| CCLE assay metabolites | 225 | 76 | 149 | 33.8% |
| CCLE GPR signatures | 195 | 169 | 26 | 86.7% |
| CCLE primary prediction targets | 60 | 60 | 0 | 100.0% |

The CCLE conclusion is explicitly limited to the mapped, reaction-accessible assay subset. Rejected
rows are retained in `ccle-metabolite-map.csv`; no inferred identity is assigned to an unmapped
assay feature.

## Supplementary Table 8 | Machine-readable claim ledger

| Claim | Decision | Basis |
| --- | --- | --- |
| Compact lipid-family scores reconstruct hidden human markers | Supported versus population mean | Participant-bootstrap RMSE interval above zero |
| Broad lipid-family scores are pathway-specific | Not supported | Size-matched random-group gate failed |
| Structure-resolved lipid scores beat a fixed-degree null | Adaptive same-cohort support | Participant-paired RMSE and precision intervals above zero |
| Structure-resolved scores replicate in a separate human cohort | External population-held-out support | ST000818 passed all 20 graph nulls |
| Direct HumanGEM/GPR features beat size-matched random features | Supported in CCLE | Lineage-isolated target-bootstrap comparison |
| Transcript interactions improve the general direct-network model | Not supported | Interaction model did not beat direct-metabolite or additive models |
| Selected transcript-conditioned coupling can exist | Supported in model system | CDA–cytidine–uridine state passed the declared screen and checks |
| The evidence validates same-person genomic personalization | Not tested | No cohort contains the required participant-level modalities |
| The model estimates physiological flux or treatment response | Not supported | Constraint evidence is model-relative; external challenges are null |
| Human structural results survive declared settings | Post-result support | Nine of nine settings pass in each cohort |
| CCLE reaction topology survives declared settings | Post-result support | Ten of ten settings pass the matched-random gate |
| CCLE topology is stable across random feature draws | Post-result support | Twenty of twenty seed intervals pass; expected-null lower bound 0.093 SD |
| CCLE topology beats degree-and-coverage-matched features | Post-result support | Twenty of twenty hard-null intervals pass; expected-null lower bound 0.097 SD |
| CCLE topology is not concentrated in one subsystem | Post-result support | Twenty-cluster lower bound 0.086 SD; all leave-one-subsystem effects pass |
| Transcript interactions become generally useful in sensitivity analysis | Not supported | Zero of ten settings beat the additive model |
| Analysis-facing mappings and exclusions are auditable | Resource supported | Complete accepted/rejected tables and verified checksums |

## Supplementary Table 9 | Evidence domains not pooled with primary observations

| Domain | Result | Adjudication |
| --- | --- | --- |
| METSIM–HumanGEM topology | 26 of 71 direct recoveries; plus-one empirical P=9.9999×10⁻⁶ | Supportive positive-control audit; nominations are knowledge based |
| CCLE reaction-resolved interaction | One FDR-significant state among 284 tests | Candidate-specific model-system evidence |
| UKB/BBJ summary-statistic scan | One strict replicated row; no novel mutable-context rows | Limited |
| NSCLC isotope challenge | No expression or oncogenotype FDR signal | Null external challenge |
| HumanGEM constraint ensemble | Three stable and two unstable configured candidates in 1,000 draws | Model-relative only |
| ST002081 temporal direction | Forward R²=0.452; reverse R²=0.459 | Directional prediction not supported |

## Supplementary Figure 1 | Factorized benchmark summary

The publication-synthesis figure compares the principal matched-null effects across evidence
domains. Effect definitions differ by benchmark and are labeled in the source data; they must not
be interpreted as a pooled meta-analytic effect. Source:
`artifacts/factorized-pathway-publication/factorized-pathway-benchmarks.png`.

## Supplementary Data index

| File | Contents |
| --- | --- |
| `artifacts/pathway-score-human-sensitivity/sensitivity-grid.csv` | Complete 18-row human parameter grid with paired intervals and degree checks |
| `artifacts/pathway-score-human-sensitivity/sensitivity-summary.csv` | Cohort-level grid adjudication |
| `artifacts/pathway-score-human-graph-mixing/null-distance-audit.csv` | Complete 200-row graph movement and degree audit |
| `artifacts/pathway-score-human-graph-mixing/mask-diversity-summary.csv` | Per-panel uniqueness and pairwise null Jaccard summaries |
| `artifacts/pathway-score-ccle-sensitivity/setting-summary.csv` | Complete ten-setting CCLE adjudication |
| `artifacts/pathway-score-ccle-sensitivity/model-metrics-grid.csv` | Metrics for all six models in every CCLE setting |
| `artifacts/pathway-score-ccle-null-ensemble/seed-summary.csv` | Complete 20-seed stochastic-null results |
| `artifacts/pathway-score-ccle-null-ensemble/target-null-ensemble.csv` | Target-level random and biological RMSE for every seed |
| `artifacts/pathway-score-ccle-property-matched-null/seed-summary.csv` | Degree- and coverage-matched null results and balance metrics |
| `artifacts/pathway-score-ccle-property-matched-null/target-null-ensemble.csv` | Hard-null target-level RMSE values |
| `artifacts/pathway-score-ccle-reaction-cluster/target-effects-with-subsystem.csv` | Deterministic target-to-subsystem assignments and effects |
| `artifacts/pathway-score-ccle-reaction-cluster/leave-one-subsystem-out.csv` | Complete leave-one-subsystem-out effects |
| `artifacts/pathway-score-publication-figures/figure-source-index.csv` | Source path and checksum for every main-figure input |
| `artifacts/pathway-score-publication-figures/manifest.json` | Figure implementation, configuration and output checksums |
| `artifacts/pathway-score-mapping-supplement/human-lipid-feature-map.csv` | Human feature eligibility, descriptors and hidden panels |
| `artifacts/pathway-score-mapping-supplement/ccle-metabolite-map.csv` | Accepted and rejected CCLE assay-to-HumanGEM mappings |
| `artifacts/pathway-score-mapping-supplement/ccle-gpr-map.csv` | GPR signatures and eligibility |
| `artifacts/pathway-score-mapping-supplement/ccle-target-feature-sets.csv` | Exact target-level feature sets |
| `artifacts/pathway-score-mapping-supplement/mapping-qc-summary.csv` | Mapping totals |
| `artifacts/factorized-pathway-publication/evidence-domain-summary.csv` | Evidence roles and adjudications |
| `artifacts/factorized-pathway-publication/claim-ledger.csv` | Machine-readable claim decisions and bases |
| `artifacts/factorized-pathway-publication/manifest.json` | Synthesis decision, source manifests and output checksums |

## Supplementary boundary statement

The analyses qualify representations for hidden-measurement reconstruction and candidate
prioritization. They do not establish pathway activity, reaction direction, physiological flux,
causal genotype effects, diagnosis, treatment response, drug efficacy or same-person genomic
personalization.
