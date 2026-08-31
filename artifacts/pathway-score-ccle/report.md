# CCLE factorized pathway held-out-metabolite benchmark

**Decision:** `RESOURCE_ONLY`  
**Human personalized prediction supported:** false  
**Measured flux supported:** false

## Result

The locked benchmark evaluated 60 target metabolites across 18 lineage groups using two repeats of five lineage-isolated folds.

| model | targets | mean_equal_lineage_rmse_sd | median_target_r2 | positive_r2_targets |
| --- | --- | --- | --- | --- |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 |
| factorized_interaction_ridge | 60 | 0.7968 | 0.241 | 52 |
| network_additive_ridge | 60 | 0.7877 | 0.2698 | 53 |
| network_metabolites_ridge | 60 | 0.7985 | 0.2406 | 53 |
| population_mean | 60 | 0.9853 | 0.0 | 0 |
| random_factorized_ridge | 60 | 0.9203 | 0.0447 | 39 |

## Paired target bootstrap

- `factorized_interaction_ridge` versus `population_mean`: mean equal-lineage RMSE improvement 0.1885 SD (95% interval 0.1429 to 0.2356).
- `factorized_interaction_ridge` versus `random_factorized_ridge`: mean equal-lineage RMSE improvement 0.1235 SD (95% interval 0.0826 to 0.1668).
- `factorized_interaction_ridge` versus `network_metabolites_ridge`: mean equal-lineage RMSE improvement 0.0017 SD (95% interval -0.0099 to 0.0139).
- `factorized_interaction_ridge` versus `network_additive_ridge`: mean equal-lineage RMSE improvement -0.0091 SD (95% interval -0.0144 to -0.0042).
- `factorized_interaction_ridge` versus `all_metabolites_ridge`: mean equal-lineage RMSE improvement -0.1588 SD (95% interval -0.2080 to -0.1146).

## Interpretation boundary

Lineage-isolated prediction in cancer cell lines can test whether compact direct-reaction and GPR features carry pathway-specific information. It does not establish human germline effects, personalized prediction, metabolite conversion rates, physiological flux, drug response or clinical utility.
