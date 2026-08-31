# Locked external replication protocol: structure-resolved pathway scores in ST000818

**Protocol version:** 1.0  
**Freeze date:** 30 August 2026, before retrieval of the abundance matrix  
**Role:** separate public-cohort human replication

## Study and eligibility

Metabolomics Workbench ST000818 was selected from the already generated public intervention
registry using metadata and metabolite names only. It contains 450 blood samples from 450 reported
participants, evenly distributed across 15 African population categories. Analysis AN001299 is
locked because the registry reports the largest positive-mode structurally named lipid panel in
this study.

The abundance matrix was not retrieved or inspected before this protocol was written. All 450
samples with a nonmissing `Categorization` factor are eligible. A feature must:

- occur in AN001299;
- have a name recognized by the frozen lipid parser;
- be observed in every eligible sample; and
- belong to a lipid family with at least five complete features.

No feature is selected by variance, association, reconstruction performance or group difference.
The analysis proceeds only if at least 100 features and 10 population categories remain.

## Representation and hidden panels

The frozen structural parser represents headgroup, total carbon, total unsaturation,
headgroup-by-unsaturation, acyl identity, chain carbon and chain unsaturation. Descriptors must
occur in at least five features and no more than 90% of eligible features.

Seeded family-stratified assignment creates five disjoint hidden panels; every feature is an outcome
exactly once. Structural scores are medians of visible, training-standardized members. The hidden
feature never contributes to a predictor score.

## Population-held-out validation

Five deterministic outer folds isolate complete population categories: no category represented in
test data appears in training data. All centering, scaling and ridge fitting use training folds
only. Ridge alpha is fixed at 10.

Models are:

1. training-fold population mean;
2. coarse lipid-family scores;
3. randomly reassigned coarse groups of equal size;
4. all visible markers;
5. declared structural descriptors; and
6. degree-preserving randomized structural descriptors.

The structural null is repeated over 20 prespecified graph seeds. Every random graph must preserve
the exact feature and descriptor degree sequences. No null realization may be selected or omitted.

## Endpoints and promotion gates

The primary endpoint is equal-population-weighted RMSE in training-fold standard-deviation units.
Precision@3 for sample-specific top-decile hidden deviations is coprimary for prioritization.
Uncertainty is paired by population category with 5,000 bootstrap draws for every graph seed.

External structural replication requires:

- all eligibility gates;
- exact graph-degree preservation for every seed and hidden panel;
- an RMSE improvement interval wholly above zero for structural versus randomized descriptors in
  every graph realization; and
- a precision@3 improvement interval wholly above zero in every graph realization.

If RMSE passes but precision does not, only external reconstruction replicates. Failure is retained
and reported. This study contains metabolomics and population metadata, not participant genomics or
transcriptomics. It cannot validate genotype-conditioned scores, temporal direction, flux,
diagnosis, treatment response or clinical utility.
