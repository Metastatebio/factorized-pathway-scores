# Blinded external-review packet

## Purpose

This packet is designed for independent pre-submission review. It does not claim that external
review has occurred. Reviewers should assess the manuscript without relying on Metastate product
claims or commercial context.

## Files to provide

1. An author-free review export generated from
   `factorized-pathway-scores-nature-communications.md` — do not send the identified submission PDF.
2. `factorized-pathway-scores-supplement.pdf` — supplementary methods and tables.
3. `artifacts/factorized-pathway-publication/claim-ledger.csv` — claim-by-claim adjudication.
4. `artifacts/pathway-score-human-sensitivity/sensitivity-grid.csv` — complete human grid.
5. `artifacts/pathway-score-ccle-sensitivity/setting-summary.csv` — complete CCLE grid.
6. `artifacts/pathway-score-mapping-supplement/` — accepted/rejected mapping ledgers.
7. A clean, read-only repository snapshot and exact environment instructions.

Remove author names, correspondence, acknowledgements, contributions, competing interests,
repository owner identifiers and document metadata before a genuinely blinded review. The current
submission PDF is identified and must not be represented as anonymized.

## Reviewer brief

Please return a recommendation—advance, major revision or stop—followed by numbered major and minor
comments. Evaluate:

- whether out-of-group hidden-marker reconstruction is a scientifically useful complement to
  enrichment and single-sample pathway scoring;
- whether each matched null controls the relevant representation complexity;
- whether participant-, population- and lineage-level validation prevents material leakage;
- whether bootstrap units and uncertainty statements match the estimands;
- whether lipid-name descriptors risk reducing the result to chemical interpolation or assay
  nomenclature;
- whether CCLE mapping coverage and target selection limit generalization;
- whether HumanGEM adjacency plus GPR expression is described without implying flux or causal
  direction;
- whether the locked external cohort is sufficiently independent despite incomplete source-paper
  documentation;
- whether negative results are visible and interpreted consistently;
- whether novelty and general importance reach the target journal's threshold.

## Mandatory claim audit

For each statement below, classify it as supported, overstated, not tested or contradicted:

1. Broad lipid families predict hidden markers.
2. Broad lipid families carry pathway-specific information.
3. Lipid structural descriptors beat fixed-degree graph nulls.
4. The structural result replicates in an independent population-held-out cohort.
5. Direct reaction neighborhoods beat dimension-matched random features across held-out lineages.
6. Transcript interactions improve the general cross-omic predictor.
7. The framework validates same-person genomic personalization.
8. The framework estimates physiological flux or treatment response.

## Reproducibility checks

At minimum, verify that:

- the reported sample, feature and target counts match the mapping ledgers;
- the sensitivity summaries can be regenerated from the complete grids;
- exact row and column degrees are preserved in the human graph nulls;
- random CCLE models have the same declared feature dimensions and folds;
- no outcome-derived operation is fitted outside training data;
- the publication manifest verifies every contributing artifact checksum;
- the manuscript's strongest sentences match the machine-readable claim ledger.

## Expertise lanes

Seek at least three independent reviewers, with no recent collaboration, employment, financial or
supervisory conflict:

- metabolomics/pathway-analysis methodology;
- biostatistics and predictive-validation design;
- genome-scale metabolic modeling and GPR interpretation.

A fourth reviewer in translational metabolomics or precision health should evaluate whether the
personalized-metabolomics implications remain properly bounded.

## Confidentiality and independence statement

Reviewers should disclose any relationship to the authors, Metastate, cited datasets or competing
products. AI-assisted editorial review may supplement but must not be labeled independent human
peer review.
