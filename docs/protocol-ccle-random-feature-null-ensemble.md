# Post-result CCLE random-feature null ensemble protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the primary and parameter-sensitivity results  
**Role:** stochastic-null robustness analysis, not untouched confirmation

## Objective

Test whether the lineage-held-out advantage of direct HumanGEM and GPR features depends on the
single random feature set used by the primary benchmark. The cohort, mappings, candidate reactions,
target gates, biological features, folds, ridge penalty and GPR aggregation remain identical to the
primary setting. Only the random-feature seed changes.

## Frozen ensemble

Twenty deterministic seeds, 20262000 through 20262019, generate target-specific random metabolite,
GPR and interaction features. Every random representation preserves the corresponding biological
representation's counts of metabolite, signature and interaction terms. Declared network features
are excluded from each target's random pool where the remaining pool is large enough.

Each realization is evaluated using two repeats of five complete-lineage holdout folds. The
factorized biological result must be numerically invariant across seeds. Complete seed- and
target-level results are retained.

## Inference and gate

For every seed, the primary contrast is random-factorized RMSE minus biological-factorized RMSE,
averaged over target-level equal-lineage RMSE values. A paired target bootstrap with 2,000 draws
produces a 95% interval. The ensemble also averages each target's random RMSE across all 20 seeds
and applies a 5,000-draw paired target bootstrap to that expected-null contrast.

Robustness requires:

1. identical target sets and biological-factorized metrics across seeds;
2. positive point improvement for all 20 seeds;
3. positive lower intervals for at least 90% of seeds; and
4. a positive lower interval for the expected-null contrast.

The analysis tests stochastic random-feature selection only. It does not establish causal
reaction direction, human germline effects, physiological flux, drug response or clinical utility.
