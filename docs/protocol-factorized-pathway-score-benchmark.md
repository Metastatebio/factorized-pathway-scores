# Frozen protocol: factorized pathway-score benchmark

**Protocol version:** 1.0  
**Freeze date:** 30 August 2026  
**Data policy:** local, checksum-pinned public data only; no new controlled-access data  
**Primary role:** retrospective methods benchmark with prospectively frozen new prediction tasks

## Scientific question

Can a compact, mechanism-aware representation of a metabolic profile recover held-out
measurements better than population baselines and size-matched random feature groupings, while
retaining explicit uncertainty and refusing to treat absent omics modalities as measured?

The factorized pathway state is represented as

`P = (M, G, T, R, E, U)`

where `M` is measured metabolomic state, `G` is genetic/GPR support, `T` is transcript-state
support, `R` is constraint-derived reserve, `E` is external context or perturbation evidence and
`U` records uncertainty and unavailable modalities. No scalar is allowed to imply that all six
components were observed in the same person.

## New benchmark 1: ST002081 missing-panel reconstruction

### Cohort and grouping

- Use the locally pinned ST002081 mwTab source.
- Exclude repository subject identifier `NA`.
- Retain features complete across eligible samples.
- Define lipid families deterministically from the assay name: `TAG*` is `TAG`; otherwise the
  prefix before the first opening parenthesis is the family.
- Retain families with at least five complete features.
- Keep every sample, but isolate all samples from a participant in one outer fold.

### Mask contract

Within every retained family, seeded permutation followed by round-robin allocation creates five
disjoint masks. Every eligible feature is hidden in exactly one mask. The target feature is never
used to construct its own score.

### Models

1. `population_mean`: training-fold target mean.
2. `random_group_score_ridge`: ridge regression from median scores of a size-preserving random
   reassignment of visible features to families.
3. `family_score_ridge`: ridge regression from median standardized scores of the declared lipid
   families.
4. `all_visible_ridge`: ridge regression from every visible marker; this is a high-dimensional
   predictive reference, not the mechanistic baseline.

All imputation, centering and scaling are fitted inside the training fold. Primary evaluation uses
five participant-isolated folds, a fixed ridge alpha of 10 and one prediction for every
sample-feature pair.

### Endpoints and gates

Primary endpoints are row-weighted RMSE, equal-participant-weighted RMSE and pooled out-of-fold
R². Secondary endpoints are precision@3 and NDCG for recovering markers in the sample-specific
top 10% of absolute standardized hidden deviations.

The pathway representation passes its bounded missing-panel claim only if:

- participant-bootstrap RMSE improvement for `family_score_ridge` versus `population_mean` has
  a 95% interval wholly above zero;
- the same improvement versus `random_group_score_ridge` has a 95% interval wholly above zero;
- participant-bootstrap precision@3 improvement versus the random grouping has a 95% interval
  wholly above zero; and
- at least eight biochemical families and 300 features enter the benchmark.

Failure does not invalidate the longitudinal individuality result; it blocks a pathway-specific
missing-panel claim.

## New benchmark 2: CCLE pathway-constrained metabolite reconstruction

### Cohort and candidate family

- Use pinned CCLE 2019 metabolomics, RNA expression and lineage annotations.
- Use the complete stable-ID mapping and direct non-transport HumanGEM v2.0.0 catalogue produced
  by `config/ccle-expanded.yaml`.
- A target metabolite must have at least one assayed direct-reaction neighbor and one available
  reaction-linked expression signature.
- All rows with the target measured are eligible; all prediction missingness is imputed using
  training-fold medians.
- Validation uses deterministic, size-balanced lineage-isolated folds.

### Models

1. `population_mean`.
2. `random_factorized_ridge`: the same numbers of metabolite, expression and interaction features
   as the mechanistic model, sampled without outcome access.
3. `network_metabolites_ridge`: directly connected assayed metabolites only.
4. `network_additive_ridge`: direct metabolites plus their HumanGEM GPR expression signatures.
5. `factorized_interaction_ridge`: direct metabolites, GPR signatures and prespecified
   metabolite-by-signature products.
6. `all_metabolites_ridge`: every other mapped metabolite, as a high-dimensional reference.

The fixed primary alpha is 10. Feature construction, median imputation and scaling occur within
each training fold. Models are evaluated twice across five lineage-isolated folds.

### Endpoints and gates

Primary endpoints are equal-lineage-weighted RMSE and row-weighted R², first per target and then
aggregated with target-bootstrap uncertainty. The factorized model passes the bounded CCLE claim
only if:

- its target-bootstrap RMSE improvement over `random_factorized_ridge` has a 95% interval wholly
  above zero;
- its improvement over `network_metabolites_ridge` has a 95% interval wholly above zero;
- its median target R² is positive; and
- at least 30 target metabolites and 10 lineages are evaluated.

Comparison with `all_metabolites_ridge` is descriptive. The compact model need not beat an
all-marker predictor to demonstrate pathway specificity.

## Prior evidence incorporated after the new runs

The following already-opened results may be synthesized but are not independent confirmation of
the new tasks:

- METSIM 71-state HumanGEM topology benchmark;
- CCLE reaction-coupling inference and CDA result;
- UKB/BBJ G×environment reaction-coordinate summary scans;
- NSCLC isotope-resolved falsification;
- BCAA public perturbation reanalyses;
- HumanGEM tissue, task-safety and constraint ensembles; and
- temporal target-entry analysis.

Each evidence domain remains separate in the claim ledger. Cross-domain convergence is not
treated as same-participant multi-omics.

## Publication decisions

`PATHWAY_SCORE_BENCHMARK_SIGNAL` requires both new benchmark gates to pass. If only one passes,
the paper becomes a domain-specific benchmark with the other domain reported as a falsification.
If neither passes, the output remains a public benchmark resource and the Metastate pathway-score
claim is not promoted.

No outcome supports diagnosis, treatment selection, physiological flux, germline
personalization, clinical utility or a Nature-level same-person multi-omic claim.
