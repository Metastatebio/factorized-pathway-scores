# Post-result CCLE reaction-topology sensitivity protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the primary CCLE benchmark  
**Role:** robustness analysis, not untouched confirmation

## Objective

Test whether lineage-held-out pathway prediction depends on the maximum reaction size, exclusion of
transport reactions, GPR expression aggregation or ridge penalty. Metabolite mapping, target
eligibility, lineage folds, model families and random-control construction otherwise remain fixed.

## Frozen settings

Ten settings are evaluated:

1. primary: non-transport reactions with at most eight metabolites, limiting-subunit GPR, alpha 10;
2. maximum four metabolites per reaction;
3. maximum six metabolites per reaction;
4. maximum 12 metabolites per reaction;
5. no reaction-size ceiling;
6. direct transport reactions included at the primary size ceiling;
7. alpha 1;
8. alpha 100;
9. mean expression for multi-subunit GPRs; and
10. summed expression for multi-subunit GPRs.

Only the named parameter changes from primary. Random controls remain target-specific and preserve
the counts of metabolite, GPR and interaction features.

## Endpoints

The primary sensitivity endpoint is target-paired equal-lineage RMSE improvement of the complete
factorized feature set over its dimension-matched random control, with 2,000 bootstrap draws.
Target count, median target R² and comparisons with direct-metabolite and additive-GPR models are
reported. Robustness requires positive point improvement in every evaluable setting and a wholly
positive 95% interval in at least 80% of settings. Interaction terms are descriptive and are not
required to improve the additive model.

All settings and failures are retained. This post-result grid does not establish human physiology,
germline genetics, flux, drug response or clinical utility.
