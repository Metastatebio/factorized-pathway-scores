# CCLE reaction-topology sensitivity

**Decision:** `CCLE_REACTION_TOPOLOGY_SENSITIVITY_ROBUST`

## Setting-level results

| setting | maximum_metabolites_per_reaction | include_transport | ridge_alpha | gpr_aggregation | candidate_rows | targets | signatures | mean_factorized_rmse_sd | median_factorized_target_r2 | positive_r2_targets | rmse_improvement_vs_random_sd | random_ci_lower | random_ci_upper | rmse_improvement_vs_additive_sd | additive_ci_lower | additive_ci_upper |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary | 8.0 | False | 10.0 | limiting_subunit | 317 | 60 | 169 | 0.7968 | 0.241 | 52 | 0.1235 | 0.0839 | 0.1662 | -0.0091 | -0.0147 | -0.0044 |
| reaction_max_4 | 4.0 | False | 10.0 | limiting_subunit | 139 | 49 | 75 | 0.8638 | 0.1933 | 40 | 0.084 | 0.0279 | 0.1347 | -0.0079 | -0.0145 | -0.0017 |
| reaction_max_6 | 6.0 | False | 10.0 | limiting_subunit | 263 | 59 | 154 | 0.8203 | 0.2365 | 52 | 0.0988 | 0.0603 | 0.1402 | -0.0094 | -0.0149 | -0.0043 |
| reaction_max_12 | 12.0 | False | 10.0 | limiting_subunit | 369 | 60 | 198 | 0.7949 | 0.2527 | 52 | 0.131 | 0.0882 | 0.1755 | -0.0107 | -0.0164 | -0.0054 |
| reaction_unbounded |  | False | 10.0 | limiting_subunit | 373 | 60 | 202 | 0.7951 | 0.2527 | 52 | 0.125 | 0.0797 | 0.1712 | -0.011 | -0.0167 | -0.0056 |
| include_transport | 8.0 | True | 10.0 | limiting_subunit | 691 | 63 | 192 | 0.7086 | 0.4692 | 58 | 0.2194 | 0.1694 | 0.2733 | -0.0143 | -0.0198 | -0.0092 |
| alpha_1 | 8.0 | False | 1.0 | limiting_subunit | 317 | 60 | 169 | 0.7985 | 0.2371 | 52 | 0.1238 | 0.0832 | 0.1658 | -0.0104 | -0.0157 | -0.0052 |
| alpha_100 | 8.0 | False | 100.0 | limiting_subunit | 317 | 60 | 169 | 0.7937 | 0.2517 | 53 | 0.117 | 0.0788 | 0.1569 | -0.0056 | -0.0096 | -0.002 |
| gpr_mean | 8.0 | False | 10.0 | mean | 317 | 60 | 169 | 0.7968 | 0.241 | 52 | 0.1232 | 0.0859 | 0.1664 | -0.009 | -0.0144 | -0.0044 |
| gpr_sum | 8.0 | False | 10.0 | sum | 317 | 60 | 169 | 0.7968 | 0.241 | 52 | 0.1232 | 0.0826 | 0.1669 | -0.009 | -0.0141 | -0.0041 |

## Complete model metrics

| model | targets | mean_equal_lineage_rmse_sd | median_target_r2 | positive_r2_targets | setting |
| --- | --- | --- | --- | --- | --- |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 | primary |
| factorized_interaction_ridge | 60 | 0.7968 | 0.241 | 52 | primary |
| network_additive_ridge | 60 | 0.7877 | 0.2698 | 53 | primary |
| network_metabolites_ridge | 60 | 0.7985 | 0.2406 | 53 | primary |
| population_mean | 60 | 0.9853 | 0.0 | 0 | primary |
| random_factorized_ridge | 60 | 0.9203 | 0.0447 | 39 | primary |
| all_metabolites_ridge | 49 | 0.6341 | 0.5701 | 48 | reaction_max_4 |
| factorized_interaction_ridge | 49 | 0.8638 | 0.1933 | 40 | reaction_max_4 |
| network_additive_ridge | 49 | 0.8559 | 0.1954 | 40 | reaction_max_4 |
| network_metabolites_ridge | 49 | 0.8728 | 0.1215 | 42 | reaction_max_4 |
| population_mean | 49 | 0.9904 | 0.0 | 0 | reaction_max_4 |
| random_factorized_ridge | 49 | 0.9478 | 0.0531 | 33 | reaction_max_4 |
| all_metabolites_ridge | 59 | 0.6442 | 0.5646 | 57 | reaction_max_6 |
| factorized_interaction_ridge | 59 | 0.8203 | 0.2365 | 52 | reaction_max_6 |
| network_additive_ridge | 59 | 0.8108 | 0.2428 | 53 | reaction_max_6 |
| network_metabolites_ridge | 59 | 0.8285 | 0.2163 | 53 | reaction_max_6 |
| population_mean | 59 | 0.9854 | 0.0 | 0 | reaction_max_6 |
| random_factorized_ridge | 59 | 0.9191 | 0.0382 | 39 | reaction_max_6 |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 | reaction_max_12 |
| factorized_interaction_ridge | 60 | 0.7949 | 0.2527 | 52 | reaction_max_12 |
| network_additive_ridge | 60 | 0.7842 | 0.2781 | 53 | reaction_max_12 |
| network_metabolites_ridge | 60 | 0.7958 | 0.263 | 53 | reaction_max_12 |
| population_mean | 60 | 0.9853 | 0.0 | 0 | reaction_max_12 |
| random_factorized_ridge | 60 | 0.9259 | 0.0427 | 41 | reaction_max_12 |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 | reaction_unbounded |
| factorized_interaction_ridge | 60 | 0.7951 | 0.2527 | 52 | reaction_unbounded |
| network_additive_ridge | 60 | 0.7841 | 0.2781 | 53 | reaction_unbounded |
| network_metabolites_ridge | 60 | 0.7958 | 0.263 | 53 | reaction_unbounded |
| population_mean | 60 | 0.9853 | 0.0 | 0 | reaction_unbounded |
| random_factorized_ridge | 60 | 0.9201 | 0.0419 | 44 | reaction_unbounded |
| all_metabolites_ridge | 63 | 0.6343 | 0.5701 | 61 | include_transport |
| factorized_interaction_ridge | 63 | 0.7086 | 0.4692 | 58 | include_transport |
| network_additive_ridge | 63 | 0.6944 | 0.4762 | 59 | include_transport |
| network_metabolites_ridge | 63 | 0.7069 | 0.4044 | 61 | include_transport |
| population_mean | 63 | 0.9866 | 0.0 | 0 | include_transport |
| random_factorized_ridge | 63 | 0.928 | 0.0818 | 48 | include_transport |
| all_metabolites_ridge | 60 | 0.6449 | 0.5618 | 58 | alpha_1 |
| factorized_interaction_ridge | 60 | 0.7985 | 0.2371 | 52 | alpha_1 |
| network_additive_ridge | 60 | 0.7881 | 0.2668 | 53 | alpha_1 |
| network_metabolites_ridge | 60 | 0.7986 | 0.2388 | 53 | alpha_1 |
| population_mean | 60 | 0.9853 | 0.0 | 0 | alpha_1 |
| random_factorized_ridge | 60 | 0.9223 | 0.0434 | 37 | alpha_1 |
| all_metabolites_ridge | 60 | 0.6319 | 0.5634 | 59 | alpha_100 |
| factorized_interaction_ridge | 60 | 0.7937 | 0.2517 | 53 | alpha_100 |
| network_additive_ridge | 60 | 0.7881 | 0.276 | 54 | alpha_100 |
| network_metabolites_ridge | 60 | 0.8004 | 0.2491 | 53 | alpha_100 |
| population_mean | 60 | 0.9853 | 0.0 | 0 | alpha_100 |
| random_factorized_ridge | 60 | 0.9107 | 0.0447 | 46 | alpha_100 |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 | gpr_mean |
| factorized_interaction_ridge | 60 | 0.7968 | 0.241 | 52 | gpr_mean |
| network_additive_ridge | 60 | 0.7878 | 0.2698 | 53 | gpr_mean |
| network_metabolites_ridge | 60 | 0.7985 | 0.2406 | 53 | gpr_mean |
| population_mean | 60 | 0.9853 | 0.0 | 0 | gpr_mean |
| random_factorized_ridge | 60 | 0.92 | 0.0447 | 39 | gpr_mean |
| all_metabolites_ridge | 60 | 0.638 | 0.5673 | 58 | gpr_sum |
| factorized_interaction_ridge | 60 | 0.7968 | 0.241 | 52 | gpr_sum |
| network_additive_ridge | 60 | 0.7878 | 0.2698 | 53 | gpr_sum |
| network_metabolites_ridge | 60 | 0.7985 | 0.2406 | 53 | gpr_sum |
| population_mean | 60 | 0.9853 | 0.0 | 0 | gpr_sum |
| random_factorized_ridge | 60 | 0.92 | 0.0447 | 39 | gpr_sum |

## Boundary

This post-result grid tests robustness of lineage-held-out CCLE prediction to reaction and model definitions. It does not establish human germline effects, personalized prediction, physiological flux, drug response or clinical utility.
