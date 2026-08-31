# Post-result CCLE property-matched null ensemble protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the dimension-matched null ensemble  
**Role:** hard-null robustness analysis, not untouched confirmation

## Objective

Challenge the CCLE topology result with random feature sets matched not only in dimension but also
in outcome-independent feature propensity. The biological features, target set, folds, model,
reaction policy and GPR aggregation remain fixed.

## Matching contract

For each target, every direct-network metabolite is assigned a unique non-neighbor assay
metabolite by minimum-cost bipartite assignment on log direct-network degree and global assay
coverage. Every reaction-linked GPR signature is assigned a unique non-network signature on log
candidate-network degree. A small seeded jitter varies near-equivalent assignments without using
target abundance, target correlation or model performance. Random interaction products preserve
the biological interaction count.

Twenty deterministic seeds, 20262100 through 20262119, are fitted using the unchanged two repeats
of five complete-lineage folds. Complete target-level metrics are retained. Balance is reported as
the mean across targets of the absolute difference between biological and random feature-set means.

## Inference and gate

Per-seed target bootstraps use 2,000 draws. The expected-null contrast averages random RMSE within
target and uses 5,000 draws. Robustness requires invariant target and biological results, positive
point improvements for every seed, positive lower intervals for at least 90% of seeds, and a
positive expected-null lower interval. Mean absolute log-degree imbalance must be at most 0.35 for
metabolites and GPR signatures; mean assay-coverage imbalance must be at most 0.02.

This hard null reduces feature-selection confounding. It does not match target correlation because
doing so would use the outcome to define the null. It does not establish causal reaction direction,
physiological flux, human personalization, drug response or clinical utility.
