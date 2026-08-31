# Human structural pathway-score sensitivity

**Decision:** `HUMAN_STRUCTURAL_SENSITIVITY_ROBUST`

## Adjudication

| dataset | settings | all_degrees_preserved | all_rmse_points_positive | all_precision_points_positive | rmse_ci_pass_rate | precision_ci_pass_rate |
| --- | --- | --- | --- | --- | --- | --- |
| ST000818 | 9 | True | True | True | 1.0 | 1.0 |
| ST002081 | 9 | True | True | True | 1.0 | 1.0 |

## Complete grid

| dataset | setting | ridge_alpha | minimum_features | maximum_feature_fraction | swaps_per_edge | descriptors | structural_rmse_sd | random_rmse_sd | structural_precision_at_k | random_precision_at_k | reference | challenger | subjects | draws | rmse_improvement_sd | rmse_ci_lower | rmse_ci_upper | precision_improvement | precision_ci_lower | precision_ci_upper | row_degrees_preserved | column_degrees_preserved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ST000818 | alpha_1 | 1.0 | 5 | 0.9 | 10 | 81 | 1.2878 | 1.7878 | 0.5542 | 0.4539 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.511 | 0.0619 | 0.9027 | 0.1003 | 0.084 | 0.1167 | True | True |
| ST000818 | alpha_100 | 100.0 | 5 | 0.9 | 10 | 81 | 1.2529 | 1.7861 | 0.5673 | 0.4938 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5428 | 0.0412 | 0.9692 | 0.0735 | 0.0576 | 0.0966 | True | True |
| ST000818 | descriptor_max_080 | 10.0 | 5 | 0.8 | 10 | 81 | 1.2374 | 1.767 | 0.5739 | 0.4717 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5406 | 0.0622 | 0.9602 | 0.1022 | 0.0847 | 0.1206 | True | True |
| ST000818 | descriptor_max_095 | 10.0 | 5 | 0.95 | 10 | 81 | 1.2374 | 1.7867 | 0.5739 | 0.4716 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5609 | 0.0681 | 0.9932 | 0.1024 | 0.0864 | 0.119 | True | True |
| ST000818 | descriptor_min_10 | 10.0 | 10 | 0.9 | 10 | 41 | 1.7641 | 1.8205 | 0.4904 | 0.393 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.056 | 0.0417 | 0.1204 | 0.0973 | 0.0813 | 0.1141 | True | True |
| ST000818 | descriptor_min_3 | 10.0 | 3 | 0.9 | 10 | 112 | 1.1915 | 1.6824 | 0.5973 | 0.5244 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.4971 | 0.0502 | 0.8823 | 0.0729 | 0.0542 | 0.0923 | True | True |
| ST000818 | primary | 10.0 | 5 | 0.9 | 10 | 81 | 1.2374 | 1.7902 | 0.5739 | 0.4713 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5632 | 0.0566 | 1.0064 | 0.1027 | 0.0855 | 0.1207 | True | True |
| ST000818 | swaps_25 | 10.0 | 5 | 0.9 | 25 | 81 | 1.2374 | 1.7729 | 0.5739 | 0.4813 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5472 | 0.0569 | 0.9747 | 0.0926 | 0.0714 | 0.1126 | True | True |
| ST000818 | swaps_5 | 10.0 | 5 | 0.9 | 5 | 81 | 1.2374 | 1.7942 | 0.5739 | 0.4719 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 15 | 2000 | 0.5686 | 0.0645 | 1.0055 | 0.1021 | 0.0819 | 0.1246 | True | True |
| ST002081 | alpha_1 | 1.0 | 5 | 0.9 | 10 | 95 | 0.3106 | 0.4081 | 0.8852 | 0.7394 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1022 | 0.0966 | 0.1078 | 0.1418 | 0.129 | 0.1554 | True | True |
| ST002081 | alpha_100 | 100.0 | 5 | 0.9 | 10 | 95 | 0.3264 | 0.4318 | 0.878 | 0.7334 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1054 | 0.0988 | 0.1121 | 0.1392 | 0.1248 | 0.1531 | True | True |
| ST002081 | descriptor_max_080 | 10.0 | 5 | 0.8 | 10 | 95 | 0.3091 | 0.4074 | 0.8877 | 0.7504 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1007 | 0.0954 | 0.1064 | 0.133 | 0.1193 | 0.1463 | True | True |
| ST002081 | descriptor_max_095 | 10.0 | 5 | 0.95 | 10 | 95 | 0.3091 | 0.4098 | 0.8877 | 0.7476 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1032 | 0.0971 | 0.1097 | 0.1398 | 0.1267 | 0.1534 | True | True |
| ST002081 | descriptor_min_10 | 10.0 | 10 | 0.9 | 10 | 69 | 0.3431 | 0.4484 | 0.8601 | 0.6906 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1089 | 0.1027 | 0.1156 | 0.1785 | 0.1614 | 0.1953 | True | True |
| ST002081 | descriptor_min_3 | 10.0 | 3 | 0.9 | 10 | 112 | 0.2899 | 0.3753 | 0.9043 | 0.7943 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.0887 | 0.0844 | 0.093 | 0.1078 | 0.095 | 0.12 | True | True |
| ST002081 | primary | 10.0 | 5 | 0.9 | 10 | 95 | 0.3091 | 0.4105 | 0.8877 | 0.7475 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.1037 | 0.0983 | 0.1091 | 0.1326 | 0.1191 | 0.1478 | True | True |
| ST002081 | swaps_25 | 10.0 | 5 | 0.9 | 25 | 95 | 0.3091 | 0.399 | 0.8877 | 0.7537 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.0933 | 0.0877 | 0.0993 | 0.1263 | 0.1144 | 0.139 | True | True |
| ST002081 | swaps_5 | 10.0 | 5 | 0.9 | 5 | 95 | 0.3091 | 0.4036 | 0.8877 | 0.7588 | degree_preserving_random_structural_ridge | structural_descriptor_ridge | 112 | 2000 | 0.097 | 0.0922 | 0.1021 | 0.1265 | 0.1128 | 0.1399 | True | True |

## Boundary

This post-result grid tests robustness of structural-score prediction to reasonable model and descriptor settings in two public cohorts. It is not prospective confirmation and does not establish genetics, physiological flux, diagnosis, treatment response or clinical utility.
